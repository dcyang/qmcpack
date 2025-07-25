//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2025 QMCPACK developers.
//
// File developed by: Peter W. Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File refactored from: EstimatorManagerBase.cpp
//////////////////////////////////////////////////////////////////////////////////////

#include <array>
#include <functional>
#include <numeric>
#include <algorithm>
#include <variant>
#include <cstdint>

#include "EstimatorManagerNew.h"
#include "EstimatorInputDelegates.h"
#include "SpinDensityNew.h"
#include "MomentumDistribution.h"
#include "OneBodyDensityMatrices.h"
#include "SelfHealingOverlap.h"
#include "MagnetizationDensity.h"
#include "PerParticleHamiltonianLogger.h"
#include "EnergyDensityEstimator.h"
#include "QMCHamiltonians/QMCHamiltonian.h"
#include "Message/Communicate.h"
#include "Message/CommOperators.h"
#include "Message/CommUtilities.h"
#include <Pools/PooledData.h>
#include "Estimators/StructureFactorEstimator.h"
#include "Estimators/LocalEnergyEstimator.h"
#include "Estimators/LocalEnergyOnlyEstimator.h"
#include "Estimators/RMCLocalEnergyEstimator.h"
#include "QMCDrivers/SimpleFixedNodeBranch.h"
#include "QMCDrivers/WalkerProperties.h"
#include "Utilities/IteratorUtility.h"
#include "hdf/hdf_archive.h"
#include "OhmmsData/AttributeSet.h"
#include "Estimators/CSEnergyEstimator.h"
#include "type_traits/variant_help.hpp"

//leave it for serialization debug
//#define DEBUG_ESTIMATOR_ARCHIVE

namespace qmcplusplus
{

using SPOMap = SPOSet::SPOMap;

EstimatorManagerNew::EstimatorManagerNew(const QMCHamiltonian& ham, Communicate* c)
    : RecordCount(0), my_comm_(c), max4ascii(8), FieldWidth(20)
{
  app_log() << " Legacy constructor adding a default LocalEnergyEstimator for the MainEstimator " << std::endl;
  max4ascii = ham.sizeOfObservables() + 3;
  addMainEstimator(std::make_unique<LocalEnergyEstimator>(ham, true));
}

bool EstimatorManagerNew::areThereListeners() const
{
  return std::any_of(operator_ests_.begin(), operator_ests_.end(),
                     [](auto& oper_est) { return oper_est->isListenerRequired(); });
}

template<class EstInputType, typename... Args>
bool EstimatorManagerNew::createEstimator(EstimatorInput& input, Args&&... args)
{
  if (has<EstInputType>(input))
  {
    operator_ests_.push_back(std::make_unique<typename EstInputType::Consumer>(std::move(std::get<EstInputType>(input)),
                                                                               std::forward<Args>(args)...));
    return true;
  }
  else
    return false;
}

template<class EstInputType, typename... Args>
bool EstimatorManagerNew::createScalarEstimator(ScalarEstimatorInput& input, Args&&... args)
{
  if (has<EstInputType>(input))
  {
    auto estimator = std::make_unique<typename EstInputType::Consumer>(std::move(std::get<EstInputType>(input)),
                                                                       std::forward<Args>(args)...);
    if (estimator->isMainEstimator())
      addMainEstimator(std::move(estimator));
    else
      scalar_ests_.push_back(std::move(estimator));
    return true;
  }
  else
    return false;
}

//initialize the name of the primary estimator
void EstimatorManagerNew::constructEstimators(EstimatorManagerInput&& emi,
                                              const ParticleSet& pset,
                                              const TrialWaveFunction& twf,
                                              const QMCHamiltonian& H,
                                              const PSPool& pset_pool)
{
  for (auto& est_input : emi.get_estimator_inputs())
    if (!(createEstimator<SpinDensityInput>(est_input, pset.getLattice(), pset.getSpeciesSet()) ||
          createEstimator<MomentumDistributionInput>(est_input, pset.getTotalNum(), pset.getTwist(),
                                                     pset.getLattice()) ||
          createEstimator<SelfHealingOverlapInput>(est_input, twf) ||
          createEstimator<OneBodyDensityMatricesInput>(est_input, pset.getLattice(), pset.getSpeciesSet(),
                                                       twf.getSPOMap(), pset) ||
          createEstimator<MagnetizationDensityInput>(est_input, pset.getLattice()) ||
          createEstimator<PerParticleHamiltonianLoggerInput>(est_input, my_comm_->rank()) ||
          createEstimator<EnergyDensityInput>(est_input, pset_pool) ||
          createEstimator<StructureFactorInput>(est_input, pset_pool)))
      throw UniformCommunicateError(std::string(error_tag_) +
                                    "cannot construct an estimator from estimator input object.");

  for (auto& scalar_input : emi.get_scalar_estimator_inputs())
    if (!(createScalarEstimator<LocalEnergyInput>(scalar_input, H) ||
          createScalarEstimator<CSLocalEnergyInput>(scalar_input, H) ||
          createScalarEstimator<RMCLocalEnergyInput>(scalar_input, H)))
      throw UniformCommunicateError(std::string(error_tag_) +
                                    "cannot construct a scalar estimator from scalar estimator input object.");

  if (main_estimator_ == nullptr)
  {
    app_log() << "  Adding a default LocalEnergyEstimator for the MainEstimator " << std::endl;
    max4ascii = H.sizeOfObservables() + 3;
    addMainEstimator(std::make_unique<LocalEnergyEstimator>(H, true));
  }

  makeConfigReport(app_log());
}

EstimatorManagerNew::~EstimatorManagerNew() = default;

/** reset names of the properties
 *
 * \todo this should be in in constructor object shouldn't be reused
 * Warning this is different from some "resets" in the code, it does not clear the object
 *
 * The number of estimators and their order can vary from the previous state.
 * reinitialized properties before setting up a new BlockAverage data list.
 *
 * The object is still not completely valid.
 *
 */
void EstimatorManagerNew::reset()
{
  //It's essential that this one is first, other code assumes weightInd always == 0
  weightInd = BlockProperties.add("BlockWeight");
  assert(weightInd == 0);
  cpuInd         = BlockProperties.add("BlockCPU");
  acceptRatioInd = BlockProperties.add("AcceptRatio");
  BlockAverages.clear(); //cleaup the records
  // Side effect of this is that BlockAverages becomes size to number of scalar values tracked by all the
  // scalar estimators
  main_estimator_->add2Record(BlockAverages);
  for (int i = 0; i < scalar_ests_.size(); i++)
    scalar_ests_[i]->add2Record(BlockAverages);
  // possibly redundant variable
  max4ascii += BlockAverages.size();
}

void EstimatorManagerNew::makeConfigReport(std::ostream& os) const
{
  os << "EstimatorManager setup for this section:\n"
     << "  Main Estimator:  " << main_estimator_->getSubTypeStr() << '\n';
  if (scalar_ests_.size() > 0)
  {
    os << "  ScalarEstimators:\n";
    for (auto& scalar_est : scalar_ests_)
      os << "    " << scalar_est->getSubTypeStr() << '\n';
  }
  if (operator_ests_.size() > 0)
  {
    os << "  General Estimators:\n";
    for (auto& est : operator_ests_)
      os << "    " << est->getMyType() << "  (" << est->getMyName() << ")\n";
  }
}

void EstimatorManagerNew::addHeader(std::ostream& o)
{
  o.setf(std::ios::scientific, std::ios::floatfield);
  o.setf(std::ios::left, std::ios::adjustfield);
  o.precision(10);
  for (int i = 0; i < BlockAverages.size(); i++)
    FieldWidth = std::max(FieldWidth, BlockAverages.Names[i].size() + 2);
  for (int i = 0; i < BlockProperties.size(); i++)
    FieldWidth = std::max(FieldWidth, BlockProperties.Names[i].size() + 2);
  int maxobjs = std::min(BlockAverages.size(), max4ascii);
  o << "#   index    ";
  for (int i = 0; i < maxobjs; i++)
    o << std::setw(FieldWidth) << BlockAverages.Names[i];
  for (int i = 0; i < BlockProperties.size(); i++)
    o << std::setw(FieldWidth) << BlockProperties.Names[i];
  o << std::endl;
  o.setf(std::ios::right, std::ios::adjustfield);
}

/// \todo clean up this method its a mess
void EstimatorManagerNew::startDriverRun()
{
  reset();
  RecordCount = 0;
  energyAccumulator.clear();
  varAccumulator.clear();
  BlockAverages.setValues(0.0);
  AverageCache.resize(BlockAverages.size());
  PropertyCache.resize(BlockProperties.size());
  // Now Estimatormanager New is actually valid i.e. in the state you would expect after the constructor.
  // Until the put is dropped this isn't feasible to fix.
#if defined(DEBUG_ESTIMATOR_ARCHIVE)
  if (!DebugArchive)
  {
    std::array<char, 128> fname;
    if (std::snprintf(fname.data(), fname.size(), "%s.p%03d.scalar.dat", my_comm_->getName().c_str(),
                      my_comm_->rank()) < 0)
      throw std::runtime_error("Error generating filename");
    DebugArchive = std::make_unique<std::ofstream>(fname.data());
    addHeader(*DebugArchive);
  }
#endif
  if (my_comm_->rank() == 0)
  {
    std::filesystem::path fname(my_comm_->getName());
    fname.concat(".scalar.dat");
    Archive = std::make_unique<std::ofstream>(fname);
    addHeader(*Archive);
    if (h5desc.size())
    {
      h5desc.clear();
    }
    fname  = my_comm_->getName() + ".stat.h5";
    h_file = std::make_unique<hdf_archive>();
    h_file->create(fname);
    main_estimator_->registerObservables(h5desc, *h_file);
    for (int i = 0; i < scalar_ests_.size(); i++)
      scalar_ests_[i]->registerObservables(h5desc, *h_file);
    for (auto& uope : operator_ests_)
      uope->registerOperatorEstimator(*h_file);
  }
}

void EstimatorManagerNew::stopDriverRun() { h_file.reset(); }

void EstimatorManagerNew::startBlock(int steps) { block_timer_.restart(); }

void EstimatorManagerNew::stopBlock(unsigned long accept, unsigned long reject, FullPrecRealType block_weight)
{
  /* Need a redesign of how accept, reject and block_weight are handled from driver to this manager.
   * DMC needs to add non-local move counters.
   * also need to add num_samples which differs from block_weight
   */
  //take block averages and update properties per block
  PropertyCache[weightInd] = block_weight;
  makeBlockAverages(accept, reject);
  reduceOperatorEstimators();
  writeOperatorEstimators();
  zeroOperatorEstimators();
  // intentionally put after all the estimator I/O
  PropertyCache[cpuInd] = block_timer_.elapsed();
  writeScalarH5();
  RecordCount++;
}

void EstimatorManagerNew::collectMainEstimators(const RefVector<ScalarEstimatorBase>& main_estimators)
{
  AverageCache = 0.0;
  for (ScalarEstimatorBase& est : main_estimators)
    est.addAccumulated(AverageCache.begin());
}

void EstimatorManagerNew::collectScalarEstimators(const std::vector<RefVector<ScalarEstimatorBase>>& crowd_scalar_ests)
{
  assert(crowd_scalar_ests[0].size() == scalar_ests_.size());
  for (int iop = 0; iop < scalar_ests_.size(); ++iop)
  {
    RefVector<ScalarEstimatorBase> this_scalar_est_for_all_crowds;
    for (int icrowd = 0; icrowd < crowd_scalar_ests.size(); ++icrowd)
      this_scalar_est_for_all_crowds.emplace_back(crowd_scalar_ests[icrowd][iop]);
    // There is actually state in each Scalar estimator that tells it what it's "first" index in the AverageCache
    // is, It's quite unclear to me if that would really work so I'm doing this in anticipation of having to rework that
    // if we don't drop multiple scalar estimators per crowd all together.
    for (ScalarEstimatorBase& est : this_scalar_est_for_all_crowds)
      est.addAccumulated(AverageCache.begin());
  }
}

void EstimatorManagerNew::collectOperatorEstimators(std::vector<RefVector<OperatorEstBase>>& crowd_op_ests)
{
  for (int iop = 0; iop < operator_ests_.size(); ++iop)
  {
    RefVector<OperatorEstBase> this_op_est_for_all_crowds;
    for (int icrowd = 0; icrowd < crowd_op_ests.size(); ++icrowd)
      this_op_est_for_all_crowds.emplace_back(crowd_op_ests[icrowd][iop]);
    operator_ests_[iop]->collect(this_op_est_for_all_crowds);
    operator_ests_[iop]->zero(this_op_est_for_all_crowds);
  }
}

void EstimatorManagerNew::reduceBlockData()
{
  //Transfer FullPrecisionRead data
  const size_t n1 = AverageCache.size();
  const size_t n2 = n1 + PropertyCache.size();

  // This is a hack but it needs to be the correct size

  std::vector<double> reduce_buffer(n2, 0.0);
  {
    auto cur = reduce_buffer.begin();
    copy(AverageCache.begin(), AverageCache.end(), cur);
    copy(PropertyCache.begin(), PropertyCache.end(), cur + n1);
  }

  // This is necessary to use mpi3's C++ style reduce
#ifdef HAVE_MPI
  my_comm_->comm.reduce_in_place_n(reduce_buffer.begin(), reduce_buffer.size(), std::plus<>{});
#endif
  if (my_comm_->rank() == 0)
  {
    auto cur = reduce_buffer.begin();
    copy(cur, cur + n1, AverageCache.begin());
    copy(cur + n1, cur + n2, PropertyCache.begin());
    const FullPrecRealType invTotWgt = 1.0 / PropertyCache[weightInd];
    AverageCache *= invTotWgt;
    //do not weight weightInd i.e. its index 0!
    for (int i = 1; i < PropertyCache.size(); i++)
      PropertyCache[i] *= invTotWgt;
  }
}

void EstimatorManagerNew::makeBlockAverages(unsigned long accepts, unsigned long rejects)
{
  // accumulate unsigned long counters over ranks.
  // Blocks ends are infrequent, two MPI transfers to preserve type
  // these could be replaced with a singple call MPI_struct_type some packing scheme or even
  // a pack into and out of an fp type that can be assured to hold the integral type exactly
  // IMHO they should not be primarily stored in a vector with magic indexes
  std::vector<unsigned long> accepts_and_rejects(my_comm_->size() * 2, 0);
  accepts_and_rejects[my_comm_->rank()]                    = accepts;
  accepts_and_rejects[my_comm_->size() + my_comm_->rank()] = rejects;
  my_comm_->allreduce(accepts_and_rejects);
  int64_t total_block_accept =
      std::accumulate(accepts_and_rejects.begin(), accepts_and_rejects.begin() + my_comm_->size(), int64_t(0));
  int64_t total_block_reject = std::accumulate(accepts_and_rejects.begin() + my_comm_->size(),
                                               accepts_and_rejects.begin() + my_comm_->size() * 2, int64_t(0));

  reduceBlockData();

  // now we put the correct accept ratio in
  PropertyCache[acceptRatioInd] = static_cast<FullPrecRealType>(total_block_accept) /
      static_cast<FullPrecRealType>(total_block_accept + total_block_reject);

  //add the block average to summarize
  energyAccumulator(AverageCache[0]);
  varAccumulator(AverageCache[1]);
}

void EstimatorManagerNew::writeScalarH5()
{
  //Do not assume h_file is valid
  if (h_file)
  {
    for (int o = 0; o < h5desc.size(); ++o)
      // cheating here, remove SquaredAverageCache from API
      h5desc[o].write(AverageCache.data(), *h_file);
    h_file->flush();
  }

  if (Archive)
  {
    *Archive << std::setw(10) << RecordCount;
    int maxobjs = std::min(BlockAverages.size(), max4ascii);
    for (int j = 0; j < maxobjs; j++)
      *Archive << std::setw(FieldWidth) << AverageCache[j];
    for (int j = 0; j < PropertyCache.size(); j++)
      *Archive << std::setw(FieldWidth) << PropertyCache[j];
    *Archive << std::endl;
  }
}

void EstimatorManagerNew::reduceOperatorEstimators()
{
  if (operator_ests_.size() > 0)
  {
    std::vector<size_t> operator_data_sizes(operator_ests_.size());
    RefVector<OperatorEstBase> ref_op_ests = convertUPtrToRefVector(operator_ests_);
    for (int iop = 0; iop < operator_data_sizes.size(); ++iop)
    {
      operator_data_sizes[iop] = operator_ests_[iop]->getFullDataSize();
    }

    // Allocate buffer to hold largest Operator est data size
    // +1 larger because we put the weight in to avoid dependence of the Scalar estimators being reduced first.
    size_t nops = *(std::max_element(operator_data_sizes.begin(), operator_data_sizes.end())) + 1;
    PooledData<RealType> operator_send_buffer;
    PooledData<RealType> operator_recv_buffer;
    operator_send_buffer.reserve(nops);
    operator_recv_buffer.reserve(nops);
    for (int iop = 0; iop < operator_ests_.size(); ++iop)
    {
      auto& estimator      = *operator_ests_[iop];
      size_t adjusted_size = estimator.getFullDataSize() + 1;
      operator_send_buffer.clear(); //(adjusted_size, 0.0);
      estimator.packData(operator_send_buffer);
      auto weight = static_cast<RealType>(estimator.get_walkers_weight());
      operator_send_buffer.add(weight);
      assert(operator_send_buffer.size() == adjusted_size);
      operator_recv_buffer.resize(adjusted_size);
      // This is necessary to use mpi3's C++ style reduce
#ifdef HAVE_MPI
      my_comm_->comm.reduce_n(operator_send_buffer.begin(), adjusted_size, operator_recv_buffer.begin(), std::plus<>{},
                              0);
#else
      operator_recv_buffer = operator_send_buffer;
#endif

      // This is a crucial step where summed over weighted observable is normalized by the total weight.  For correctness this should be done
      // only after the full weighted sum is done.
      // i.e.  (1 / Sum(w_1 + ... + w_n)) * (w_1 * S_1 + ... + w_n * S_n) != (1/w_1) * w_1 * S_1 + ... + (1/w_n) * w_n * S_n
      // for a general case where w's are not strictly ==
      // Assumptions lead rank is 0, true for Ensembles?
      if (my_comm_->rank() == 0)
      {
        estimator.unpackData(operator_recv_buffer);
        RealType reduced_walker_weights = 0.0;
        operator_recv_buffer.get(reduced_walker_weights);
        RealType invTotWgt = 1.0 / static_cast<RealType>(reduced_walker_weights);
        estimator.normalize(invTotWgt);
      }
    }
  }
}

void EstimatorManagerNew::writeOperatorEstimators()
{
  if (my_comm_->rank() == 0)
  {
    if (h_file)
    {
      for (auto& op_est : operator_ests_)
        op_est->write(*h_file);
      h_file->flush();
    }
  }
}

void EstimatorManagerNew::zeroOperatorEstimators()
{
  for (auto& op_est : operator_ests_)
    op_est->zero();
}

void EstimatorManagerNew::getApproximateEnergyVariance(FullPrecRealType& e, FullPrecRealType& var)
{
  RealType tmp[3];
  tmp[0] = energyAccumulator.count();
  tmp[1] = energyAccumulator.result();
  tmp[2] = varAccumulator.result();
  my_comm_->bcast(tmp, 3);
  e   = tmp[1] / tmp[0];
  var = tmp[2] / tmp[0] - e * e;
}

void EstimatorManagerNew::addMainEstimator(std::unique_ptr<ScalarEstimatorBase>&& estimator)
{
  if (main_estimator_ != nullptr)
    app_log() << "  EstimatorManagerNew replaced its main estimator with " << estimator->getSubTypeStr()
              << " estimator." << std::endl;
  main_estimator_ = std::move(estimator);
}

int EstimatorManagerNew::addScalarEstimator(std::unique_ptr<ScalarEstimatorBase>&& estimator)
{
  scalar_ests_.push_back(std::move(estimator));
  return scalar_ests_.size() - 1;
}

} // namespace qmcplusplus
