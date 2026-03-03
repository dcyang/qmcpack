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

#include "SplineUtils.h"
#include <sstream>
#include "spline/einspline_engine.hpp"
#include "spline/einspline_util.hpp"

namespace qmcplusplus
{
template<typename ST>
bool SplineUtils<ST>::read(MultiBsplineBase<ST>& spline, hdf_archive& h5f)
{
  std::ostringstream o;
  o << "spline_" << my_index;
  einspline_engine bigtable(spline.getSplinePtr());
  return h5f.readEntry(bigtable, o.str());
}

template<typename ST>
bool SplineUtils<ST>::write(MultiBsplineBase<ST>& spline, hdf_archive& h5f)
{
  std::ostringstream o;
  o << "spline_" << my_index;
  einspline_engine bigtable(spline.getSplinePtr());
  return h5f.writeEntry(bigtable, o.str());
}

template<typename ST>
bool SplineUtils<ST>::read(MultiBspline1D<ST>& spline, hdf_archive& h5f)
{
  einspline_engine bigtable(spline.getSplinePtr());
  return h5f.readEntry(bigtable, "radial_spline");
}

template<typename ST>
bool SplineUtils<ST>::write(MultiBspline1D<ST>& spline, hdf_archive& h5f)
{
  einspline_engine bigtable(spline.getSplinePtr());
  return h5f.writeEntry(bigtable, "radial_spline");
}

template class SplineUtils<float>;
template class SplineUtils<double>;
} // namespace qmcplusplus
