//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2025 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File refactored from: OperatorBase.cpp
//////////////////////////////////////////////////////////////////////////////////////

/**@file
 */
#include "Message/Communicate.h"
#include "OperatorEstBase.h"
#include "QMCHamiltonians/QMCHamiltonian.h"

namespace qmcplusplus
{
using Real = OperatorEstBase::Real;

OperatorEstBase::OperatorEstBase(DataLocality data_locality, const std::string& name, const std::string& type)
    : data_locality_(data_locality), my_name_(name), my_type_(type), walkers_weight_(0)
{}

OperatorEstBase::OperatorEstBase(const OperatorEstBase& oth)
    : data_locality_(oth.data_locality_), my_name_(oth.my_name_), my_type_(oth.my_type_), walkers_weight_(0)
{}

void OperatorEstBase::collect(const RefVector<OperatorEstBase>& type_erased_operator_estimators)
{
  for (OperatorEstBase& crowd_oeb : type_erased_operator_estimators)
  {
    std::transform(data_.begin(), data_.end(), crowd_oeb.get_data().begin(), data_.begin(), std::plus<>{});
    walkers_weight_ += crowd_oeb.walkers_weight_;
  }
}

void OperatorEstBase::normalize(QMCT::RealType invTotWgt)
{
  for (QMCT::RealType& elem : data_)
    elem *= invTotWgt;
}

void OperatorEstBase::write(hdf_archive& file)
{
  if (h5desc_.empty())
    return;
  // We have to do this to deal with the legacy design that Observables using
  // collectables in mixed precision were accumulated in float but always written
  // to hdf5 in double.
#ifdef MIXED_PRECISION
  std::vector<QMCT::FullPrecRealType> expanded_data(data_.size(), 0.0);
  std::copy_n(data_.begin(), data_.size(), expanded_data.begin());
  assert(!data_.empty());
  // auto total = std::accumulate(data_->begin(), data_->end(), 0.0);
  // std::cout << "data size: " << data_->size() << " : " << total << '\n';
  for (auto& h5d : h5desc_)
    h5d.write(expanded_data.data(), file);
#else
  for (auto& h5d : h5desc_)
    h5d.write(data_.data(), file);
#endif
  file.pop();
}

void OperatorEstBase::packData(PooledData<Real>& buffer) const { buffer.add(data_.begin(), data_.end()); }

void OperatorEstBase::unpackData(PooledData<Real>& buffer) { buffer.get(data_.begin(), data_.end()); }

void OperatorEstBase::zero(RefVector<OperatorEstBase>& type_erased_operator_estimators) const
{
  for (OperatorEstBase& crowd_oeb : type_erased_operator_estimators)
    crowd_oeb.zero();
}

void OperatorEstBase::zero()
{
  if (data_locality_ == DataLocality::rank || data_locality_ == DataLocality::crowd)
    std::fill(data_.begin(), data_.end(), 0.0);
  else
    data_.clear();
  walkers_weight_ = 0;
}

} // namespace qmcplusplus
