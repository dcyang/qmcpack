//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2025 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////

#include "ApproximateEquality.hpp"

namespace qmcplusplus
{
template bool approxEquality<float>(float val_a, float val_b, std::optional<double> eps);
template bool approxEquality<std::complex<float>>(std::complex<float> val_a,
                                                  std::complex<float> val_b,
                                                  std::optional<double> eps);
template bool approxEquality<double>(double val_a, double val_b, std::optional<double> eps);
template bool approxEquality<std::complex<double>>(std::complex<double> val_a,
                                                   std::complex<double> val_b,
                                                   std::optional<double> eps);
} // namespace qmcplusplus
