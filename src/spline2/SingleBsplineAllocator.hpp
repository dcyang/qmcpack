//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2025 QMCPACK developers.
//
// File developed by: Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@intel.com, Intel Corp.
//////////////////////////////////////////////////////////////////////////////////////
// -*- C++ -*-
/** @file SingleBsplineAllocator.hpp
 * @brief BsplineAllocator and management classes
 */
#ifndef QMCPLUSPLUS_SINGLE_BSPLINE_ALLOCATOR_H
#define QMCPLUSPLUS_SINGLE_BSPLINE_ALLOCATOR_H

#include "spline2/bspline_traits.hpp"

extern "C"
{
  //forward declaration of the functions
  void find_coefs_1d_d(Ugrid grid, BCtype_d bc, double* data, intptr_t dstride, double* coefs, intptr_t cstride);

  void find_coefs_1d_s(Ugrid grid, BCtype_s bc, float* data, intptr_t dstride, float* coefs, intptr_t cstride);
}

namespace qmcplusplus
{
namespace einspline
{
inline void find_coefs_1d(Ugrid grid, BCtype_d bc, double* data, intptr_t dstride, double* coefs, intptr_t cstride)
{
  find_coefs_1d_d(grid, bc, data, dstride, coefs, cstride);
}

inline void find_coefs_1d(Ugrid grid, BCtype_s bc, float* data, intptr_t dstride, float* coefs, intptr_t cstride)
{
  find_coefs_1d_s(grid, bc, data, dstride, coefs, cstride);
}
} // namespace einspline

template<typename T>
class SingleBsplineAllocator
{
  using SingleSplineType = typename bspline_traits<T, 3>::SingleSplineType;
  using BCType           = typename bspline_traits<T, 3>::BCType;

public:
  static void destroy(SingleSplineType* spline)
  {
    delete[] spline->coefs;
    delete spline;
  }

  ///allocate a UBspline_3d_d, it can be made template to support UBspline_3d_s
  static SingleSplineType* allocateUBspline(Ugrid x_grid,
                                            Ugrid y_grid,
                                            Ugrid z_grid,
                                            BCType xBC,
                                            BCType yBC,
                                            BCType zBC,
                                            T* data = nullptr);
};

template<typename T>
typename SingleBsplineAllocator<T>::SingleSplineType* SingleBsplineAllocator<T>::allocateUBspline(Ugrid x_grid,
                                                                                                  Ugrid y_grid,
                                                                                                  Ugrid z_grid,
                                                                                                  BCType xBC,
                                                                                                  BCType yBC,
                                                                                                  BCType zBC,
                                                                                                  T* data)
{
  // Create new spline
  auto* spline = new SingleSplineType;

  spline->spcode = bspline_traits<T, 3>::single_spcode;
  spline->tcode  = bspline_traits<T, 3>::tcode;
  spline->xBC    = xBC;
  spline->yBC    = yBC;
  spline->zBC    = zBC;

  // Setup internal variables
  int Mx = x_grid.num;
  int My = y_grid.num;
  int Mz = z_grid.num;
  int Nx, Ny, Nz;

  if (xBC.lCode == PERIODIC || xBC.lCode == ANTIPERIODIC)
    Nx = Mx + 3;
  else
    Nx = Mx + 2;
  x_grid.delta     = (x_grid.end - x_grid.start) / (double)(Nx - 3);
  x_grid.delta_inv = 1.0 / x_grid.delta;
  spline->x_grid   = x_grid;

  if (yBC.lCode == PERIODIC || yBC.lCode == ANTIPERIODIC)
    Ny = My + 3;
  else
    Ny = My + 2;
  y_grid.delta     = (y_grid.end - y_grid.start) / (double)(Ny - 3);
  y_grid.delta_inv = 1.0 / y_grid.delta;
  spline->y_grid   = y_grid;

  if (zBC.lCode == PERIODIC || zBC.lCode == ANTIPERIODIC)
    Nz = Mz + 3;
  else
    Nz = Mz + 2;
  z_grid.delta     = (z_grid.end - z_grid.start) / (double)(Nz - 3);
  z_grid.delta_inv = 1.0 / z_grid.delta;
  spline->z_grid   = z_grid;

  spline->x_stride = Ny * Nz;
  spline->y_stride = Nz;

  spline->coefs_size = (size_t)Nx * (size_t)Ny * (size_t)Nz;

  spline->coefs = new T[spline->coefs_size];

  if (data != NULL) // only data is provided
  {
// First, solve in the X-direction
#pragma omp parallel for
    for (int iy = 0; iy < My; iy++)
      for (int iz = 0; iz < Mz; iz++)
      {
        intptr_t doffset = iy * Mz + iz;
        intptr_t coffset = iy * Nz + iz;
        einspline::find_coefs_1d(spline->x_grid, xBC, data + doffset, My * Mz, spline->coefs + coffset, Ny * Nz);
      }

// Now, solve in the Y-direction
#pragma omp parallel for
    for (intptr_t ix = 0; ix < Nx; ix++)
      for (intptr_t iz = 0; iz < Nz; iz++)
      {
        intptr_t doffset = ix * Ny * Nz + iz;
        intptr_t coffset = ix * Ny * Nz + iz;
        einspline::find_coefs_1d(spline->y_grid, yBC, spline->coefs + doffset, Nz, spline->coefs + coffset, Nz);
      }

// Now, solve in the Z-direction
#pragma omp parallel for
    for (intptr_t ix = 0; ix < Nx; ix++)
      for (intptr_t iy = 0; iy < Ny; iy++)
      {
        intptr_t doffset = (ix * Ny + iy) * Nz;
        intptr_t coffset = (ix * Ny + iy) * Nz;
        einspline::find_coefs_1d(spline->z_grid, zBC, spline->coefs + doffset, 1, spline->coefs + coffset, 1);
      }
  }
  return spline;
}

} // namespace qmcplusplus
#endif
