//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_ESTIMATOR_MANAGER_INPUT_TEST_H
#define QMCPLUSPLUS_ESTIMATOR_MANAGER_INPUT_TEST_H

#include "OhmmsData/Libxml2Doc.h"

namespace qmcplusplus
{
namespace testing
{

Libxml2Document createEstimatorManagerNewGlobalInputXML();
Libxml2Document createEstimatorNewGlobalOperatorEstInputXML();
constexpr int n_opest_new_global_input_xml{1};
Libxml2Document createEstimatorManagerNewInputXML();
constexpr int n_opest_new_input_xml{4};
Libxml2Document createEstimatorManagerNewVMCInputXML();
constexpr int n_opest_new_vmc_xml{3};
Libxml2Document createEstimatorManagerEnergyDenistyInputXML();

} // namespace testing
} // namespace qmcplusplus
#endif
