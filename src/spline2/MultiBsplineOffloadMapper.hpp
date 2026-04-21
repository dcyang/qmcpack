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
// -*- C++ -*-
/**@file MultiBsplineOffload.hpp
 *
 * define classes MultiBsplineOffload
 * The evaluation functions are defined in MultiBsplineOffloadEval.hpp
 */
#ifndef QMCPLUSPLUS_MULTIEINSPLINEOFFLOAD_HPP
#define QMCPLUSPLUS_MULTIEINSPLINEOFFLOAD_HPP

#include "MultiBsplineBase.hpp"
#include <vector>
#include "SoAFields3D.hpp"

namespace qmcplusplus
{
/** container class to hold a 3D multi spline pointer and BsplineAllocator
 * @tparam T the precision of splines
 */
template<typename T>
class MultiBsplineOffloadMapper
{
  using HostBspline = MultiBsplineBase<T>;

  const HostBspline& host_bsplines_;
  std::vector<const T*> block_coefs_;

public:
  MultiBsplineOffloadMapper(const HostBspline& host_bsplines);

  ~MultiBsplineOffloadMapper();

  virtual void mapToDevice();

  void updateToDevice();

  void mw_evaluate_v(int num_pos, T* pos_arr, T* spline_v);
  void mw_evaluate_vgh(int num_pos, T* pos_arr, T* spline_vgh);
};

extern template class MultiBsplineOffloadMapper<float>;
extern template class MultiBsplineOffloadMapper<double>;
} // namespace qmcplusplus

#endif
