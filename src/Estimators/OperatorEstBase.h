//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2020 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File refactored from: OperatorBase.h
//////////////////////////////////////////////////////////////////////////////////////


/**@file
 */
#ifndef QMCPLUSPLUS_OPERATORESTBASE_H
#define QMCPLUSPLUS_OPERATORESTBASE_H

#include "Particle/ParticleSet.h"
#include "OhmmsData/RecordProperty.h"
#include "Utilities/RandomGenerator.h"
#include "QMCHamiltonians/ObservableHelper.h"
#include "QMCHamiltonians/QMCHamiltonian.h"
#include "QMCWaveFunctions/OrbitalSetTraits.h"
#include "type_traits/DataLocality.h"
#include "hdf/hdf_archive.h"
#include <bitset>

namespace qmcplusplus
{
class TrialWaveFunction;
namespace testing
{
class OEBAccessor;
}
/** @ingroup Estimators
 * @brief An abstract class for gridded estimators
 *
 */
class OperatorEstBase
{
public:
  using QMCT             = QMCTraits;
  using FullPrecRealType = QMCT::FullPrecRealType;
  using MCPWalker        = Walker<QMCTraits, PtclOnLatticeTraits>;
  using Real             = QMCT::RealType;
  using Data             = std::vector<Real>;

  ///constructor
  OperatorEstBase(DataLocality data_locality, const std::string& name, const std::string& type);
  ///virtual destructor
  virtual ~OperatorEstBase() = default;

  /** Accumulate whatever it is you are accumulating with respect to walkers
   *
   *  This method is assumed to be called from the crowd context
   *  It provides parallelism with respect to computational effort of the estimator
   *  without causing a global sync.
   *  Depending on data locality the accumlation of the result may be different from
   *  the single thread write directly into the OperatorEstimator data.
   *  \param[in]      walkers
   *  \param[inout]   pset_target   crowd scope target pset (should be returned to starting state after call)
   *  \param[in]      psets         per walker psets
   *  \param[in]      wnfs          per walker TrialWaveFunction
   *  \param[inout]   rng           crowd scope RandomGenerator
   */
  virtual void accumulate(const RefVector<MCPWalker>& walkers,
                          const RefVector<ParticleSet>& psets,
                          const RefVector<TrialWaveFunction>& wfns,
                          const RefVector<QMCHamiltonian>& hams,
                          RandomBase<FullPrecRealType>& rng) = 0;

  /** Reduce estimator result data from crowds to rank
   *
   *  This is assumed to be called from only from one thread per crowds->rank
   *  reduction. Implied is this is during a global sync or there is a guarantee
   *  that the crowd operator estimators accumulation data is not
   *  being written to.
   *
   *  It is assumed by derived classes and it is necessary to support
   *  derived type data locality schemes that the OperatorEstBase
   *  derived types in the refvector match.
   *
   *  The input operators are not zeroed after collect is called,
   *  the owner of the operators must handle the accumulated state explicitly.
   *
   *  A side effect is walker_weights_ are collected
   *  as well, so if this is not called from an override the
   *  walker_weights_ must be collected there.  If it is called they
   *  must not be collected there.
   *
   *  There could be concurrent operations inside the scope of the collect call.
   */
  virtual void collect(const RefVector<OperatorEstBase>& oebs);

  virtual void normalize(QMCT::RealType invToWgt);

  virtual void startBlock(int steps) = 0;

  const std::vector<QMCT::RealType>& get_data() const { return data_; }
  std::vector<QMCT::RealType>& get_data() { return data_; }

  virtual std::size_t getFullDataSize() const { return data_.size(); }

  /** @ingroup Functions to add or remove estimator data from PooledData<Real>
   *  @brief   used for MPI reduction.
   *           These are only used on the rank estimator owned by EstimatorManagerNew.
   *           The rank EstimatorManagerNew owns the buffer.
   *           It is not intended to store the state of the estimator.
   *           The packing and unpacking functions must follow the same sequence of adds or gets
   *           as PooledData is a stateful sequence of bytes with an internal position cursor.
   *  @{
   */

  /** Packs data from native container types in a subtype of Operator est base
   *  to buffer of type Real for reduction over MPI.
   *  I.e. writes to pooled data.
   */
  virtual void packData(PooledData<Real>& buffer) const;
  /** Unpacks data from mpi buffer of type Real into native container types
   *  after a reduction over MPI.
   *  i.e. reads from pooled data.
   */
  virtual void unpackData(PooledData<Real>& buffer);
  ///@}

  /*** create and tie OperatorEstimator's observable_helper hdf5 wrapper to stat.h5 file
   * @param gid hdf5 group to which the observables belong
   *
   * The default implementation does nothing. The derived classes which compute
   * big data, e.g. density, should overwrite this function.
   */
  virtual void registerOperatorEstimator(hdf_archive& file) {}

  virtual std::unique_ptr<OperatorEstBase> spawnCrowdClone() const = 0;

  /** Write to previously registered observable_helper hdf5 wrapper.
   *
   *  if you haven't registered Operator Estimator
   *  this will do nothing.
   */
  virtual void write(hdf_archive& file);

  /** Calls zero on every OperatorEstBase in refvector
   *
   *  like collect this is intended to be called with a refvector
   *  where the OperatorEstBase derived types are all the same.
   *  Derived types overriding this can assume this.
   */
  virtual void zero(RefVector<OperatorEstBase>& oebs) const;

  /** zero data appropriately for the DataLocality
   *
   *  Derived classes that don't solely rely on data_ for
   *  their accumulated data must override this function.
   */
  virtual void zero();

  /** Return the total walker weight for this block
   */
  QMCT::FullPrecRealType get_walkers_weight() const { return walkers_weight_; }

  const std::string& getMyName() const { return my_name_; }
  const std::string& getMyType() const { return my_type_; }

  /** Register 0-many listeners with a leading QMCHamiltonian instance i.e. a QMCHamiltonian
   *  that has acquired the crowd scope QMCHamiltonianMultiWalkerResource.
   *  This must be called for each crowd scope estimator that listens to register listeners into
   *  the crowd scope QMCHamiltonianMultiWalkerResource.
   *
   *  Many estimators don't need per particle values so the default implementation is no op.
   */
  virtual void registerListeners(QMCHamiltonian& ham_leader) {};

  bool isListenerRequired() { return requires_listener_; }

  DataLocality get_data_locality() const { return data_locality_; }

protected:
  /** Shallow copy constructor!
   *  This alows us to keep the default copy constructors for derived classes which
   *  is quite useful to the spawnCrowdClone design. But this is a
   *  code smell for sure.
   *
   *  Data is likely to be quite large and since the OperatorEstBase design is that the children
   *  reduce to the parent it is infact undesirable for them to copy the data the parent has.
   *  Initialization of Data (i.e. call to resize) if any is the responsibility of the derived class.
   */
  OperatorEstBase(const OperatorEstBase& oth);

  /** locality for accumulation of estimator data.
   *  This designates the memory scheme used for the estimator
   *  The default is:
   *  DataLocality::Crowd, each crowd and the rank level estimator have a full representation of the data
   *  Memory Savings Schemes:
   *  One:
   *  DataLocality::Rank,  This estimator has the full representation of the data but its crowd spawn will have
   *  One per crowd:
   *  DataLocality::Queue  This estimator accumulates queue of values to collect to the Rank estimator data
   *  DataLocality::?      Another way to reduce memory use on thread/crowd local estimators.
   */
  DataLocality data_locality_;

  ///name of this object -- only used for debugging and h5 output
  std::string my_name_;
  std::string my_type_;

  QMCT::FullPrecRealType walkers_weight_;

  // convenient Descriptors hdf5 for Operator Estimators only populated for rank scope OperatorEstimator
  std::vector<ObservableHelper> h5desc_;

  Data data_;

  bool requires_listener_ = false;

  friend testing::OEBAccessor;
};
} // namespace qmcplusplus
#endif
