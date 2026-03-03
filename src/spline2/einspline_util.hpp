//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/** @file einspline_util.hpp
 * @brief utility functions for i/o and bcast of einspline objects
 *
 */
#ifndef QMCPLUSPLUS_EINSPLINE_UTILITIES_H
#define QMCPLUSPLUS_EINSPLINE_UTILITIES_H

#include "Message/CommOperators.h"
#include "mpi/mpi_datatype.h"
#include "Host/OutputManager.h"
#include "OhmmsData/FileUtility.h"
#include "hdf/hdf_archive.h"
#include "einspline/multi_bspline_copy.h"
#include "bspline_traits.hpp"
#include <limits>

namespace qmcplusplus
{
///handles i/o and bcast, testing for now
template<typename T>
inline void chunked_bcast(Communicate* comm, T* buffer, size_t ntot)
{
  if (comm->size() == 1)
    return;

  size_t chunk_size = (1 << 30) / sizeof(T); //256 MB
  int n             = static_cast<int>(ntot / chunk_size);

  size_t offset = 0;
  for (int i = 0; i < n; ++i, offset += chunk_size)
  {
    comm->bcast(buffer + offset, static_cast<int>(chunk_size));
  }

  if (offset < ntot)
  {
    comm->bcast(buffer + offset, static_cast<int>(ntot - offset));
  }
}

template<typename ENGT>
inline void chunked_bcast(Communicate* comm, ENGT* buffer)
{ chunked_bcast(comm, buffer->coefs, buffer->coefs_size); }

template<typename T, unsigned D>
inline void gatherv(Communicate* comm,
                    typename bspline_traits<T, D>::SplineType* buffer,
                    const int ncol,
                    const std::vector<int>& offset)
{
  std::vector<int> counts(offset.size() - 1, 0);
  for (size_t ib = 0; ib < counts.size(); ib++)
    counts[ib] = offset[ib + 1] - offset[ib];
  const auto& counts_const(counts);
  const size_t coef_type_bytes = sizeof(T);
  if (buffer->coefs_size * coef_type_bytes > std::numeric_limits<int>::max())
  {
    // Some MPI libraries have problems when message sizes exceed range of integer (2^31-1)
    // Perform the gatherv in columns to reduce risk
    const size_t xs = buffer->x_stride;
    if (xs * coef_type_bytes >= std::numeric_limits<int>::max())
      app_warning() << "Large single message even after splitting by the number of grid points in x direction! "
                    << "Some MPI libraries may not work!" << std::endl;
    const size_t nx         = buffer->coefs_size / xs;
    const int nrow          = buffer->coefs_size / (ncol * nx);
    MPI_Datatype columntype = mpi::construct_column_type(buffer->coefs, nrow, ncol);
    for (size_t iz = 0; iz < nx; iz++)
      comm->gatherv_in_place(buffer->coefs + xs * iz, columntype, counts_const, offset);
    mpi::free_column_type(columntype);
  }
  else
  {
    const int nrow          = buffer->coefs_size / ncol;
    MPI_Datatype columntype = mpi::construct_column_type(buffer->coefs, nrow, ncol);
    comm->gatherv_in_place(buffer->coefs, columntype, counts_const, offset);
    mpi::free_column_type(columntype);
  }
}

template<unsigned DIM>
struct dim_traits
{};

// for 3D multi
template<>
struct dim_traits<4>
{
  template<typename data_type>
  static void setdim(data_type& a, hsize_t* dims)
  {
    dims[0] = a.spliner.x_grid.num + 3;
    dims[1] = a.spliner.y_grid.num + 3;
    dims[2] = a.spliner.z_grid.num + 3;
    dims[3] = a.spliner.z_stride;
  }
};

// for 1D multi
template<>
struct dim_traits<2>
{
  template<typename data_type>
  static void setdim(data_type& a, hsize_t* dims)
  {
    dims[0] = a.spliner.x_grid.num + 2;
    dims[1] = a.spliner.x_stride;
  }
};

/** specialization of h5data_proxy for einspline_engine
   */
template<typename T, unsigned D>
struct h5data_proxy<einspline_engine<T, D>>
    : public h5_space_type<T, D + 1>
{
  using value_type = T;
  using Base       = h5_space_type<T, D + 1>;
  using Base::dims;
  using Base::get_address;
  using data_type = einspline_engine<T, D>;

  inline h5data_proxy(const data_type& a) { dim_traits<D + 1>::setdim(a, dims); }

  inline bool read(data_type& ref, hid_t grp, const std::string& aname, hid_t xfer_plist = H5P_DEFAULT)
  {
    return h5d_read(grp, aname, get_address(ref.spliner.coefs), xfer_plist);
  }

  inline bool write(const data_type& ref, hid_t grp, const std::string& aname, hid_t xfer_plist = H5P_DEFAULT) const
  { return h5d_write(grp, aname.c_str(), Base::rank, dims, get_address(ref.spliner.coefs), xfer_plist); }
};

namespace einspline
{
template<typename IV>
bool outOfBound(const IV& a)
{
  for (int i = 0; i < a.size(); ++i)
    if (a[i] < 0.0 || a[i] >= 1.0)
      return true;
  return false;
}

template<typename IV>
bool validRange(const IV& low, const IV& up)
{
  bool ok = low[0] < up[0];
  for (int i = 1; i < low.size(); ++i)
    ok &= (low[i] < up[i]);
  return ok;
}
} // namespace einspline

} // namespace qmcplusplus
#endif
