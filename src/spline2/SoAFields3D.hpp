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


#ifndef QMCPLUSPLUS_SOAFIELDS3D_HPP
#define QMCPLUSPLUS_SOAFIELDS3D_HPP

namespace qmcplusplus
{
enum SoAFields3D
{
  VAL = 0,
  GRAD0,
  GRAD1,
  GRAD2,
  HESS00,
  HESS01,
  HESS02,
  HESS11,
  HESS12,
  HESS22,
  LAPL,
  NUM_FIELDS
};
} // namespace qmcplusplus
#endif
