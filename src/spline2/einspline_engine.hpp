//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_EINSPLINE_ENGINE_HPP
#define QMCPLUSPLUS_EINSPLINE_ENGINE_HPP

#include "bspline_traits.hpp"

namespace qmcplusplus
{
/** einspline_engine
   *
   * The copy constructor is disabled.
   */
template<typename T, unsigned D>
class einspline_engine
{
public:
  using real_type  = typename bspline_traits<T, D>::real_type;
  using value_type  = typename bspline_traits<T, D>::value_type;
  using SplineType  = typename bspline_traits<T, D>::SplineType;
  ///spline engine
  SplineType& spliner;

  einspline_engine(SplineType& s) : spliner(s)
  {
  }
};
} // namespace qmcplusplus

#endif
