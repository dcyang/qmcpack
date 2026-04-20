//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2026 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//////////////////////////////////////////////////////////////////////////////////////


#include "MultiBsplineOffloadMapper.hpp"

namespace qmcplusplus
{
template<typename T>
MultiBsplineOffloadMapper<T>::MultiBsplineOffloadMapper(const HostBspline& host_bsplines)
    : host_bsplines_(host_bsplines)
{
  block_coefs_.reserve(host_bsplines_.getNumBlocks());
  for (int ib = 0; ib < host_bsplines_.getNumBlocks(); ib++)
  {
    auto* coefs = host_bsplines_.getBlock(ib).coefs;
    block_coefs_.push_back(coefs);
  }
}

template<typename T>
void MultiBsplineOffloadMapper<T>::mapToDevice()
{
  for (int ib = 0; ib < host_bsplines_.getNumBlocks(); ib++)
  {
    auto* spline_m = &host_bsplines_.getBlock(ib);
    auto* coefs    = block_coefs_[ib];
    PRAGMA_OFFLOAD("omp target enter data map(to: spline_m[:1]) map(alloc: coefs[:spline_m->coefs_size])")
  }
}

template<typename T>
MultiBsplineOffloadMapper<T>::~MultiBsplineOffloadMapper()
{
  for (int ib = 0; ib < host_bsplines_.getNumBlocks(); ib++)
  {
    auto* spline_m = &host_bsplines_.getBlock(ib);
    auto* coefs    = block_coefs_[ib];
    PRAGMA_OFFLOAD("omp target exit data map(delete: spline_m[:1]) map(delete: coefs[:spline_m->coefs_size])")
  }
}

template<typename T>
void MultiBsplineOffloadMapper<T>::updateToDevice()
{
  for (int ib = 0; ib < host_bsplines_.getNumBlocks(); ib++)
  {
    auto* spline_m = &host_bsplines_.getBlock(ib);
    auto* coefs    = block_coefs_[ib];
    PRAGMA_OFFLOAD("omp target update to(coefs[:spline_m->coefs_size])")
  }
}

template class MultiBsplineOffloadMapper<float>;
template class MultiBsplineOffloadMapper<double>;
} // namespace qmcplusplus
