//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2021 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////


#include "catch.hpp"

#include "OneBodyDensityMatricesInput.h"
#include "ValidOneBodyDensityMatricesInput.h"
#include "InvalidOneBodyDensityMatricesInput.h"
#include "EstimatorTesting.h"
#include "ParticleSet.h"
#include "OhmmsData/Libxml2Doc.h"
#include "Message/UniformCommunicateError.h"

#include <iostream>

namespace qmcplusplus
{
TEST_CASE("OneBodyDensityMatricesInput::from_xml", "[estimators]")
{
  using Input = testing::ValidOneBodyDensityMatricesInput;
  Input valid_input;
  for (auto input_xml : valid_input)
  {
    Libxml2Document doc;
    bool okay = doc.parseFromString(input_xml);
    REQUIRE(okay);
    xmlNodePtr node = doc.getRoot();
    OneBodyDensityMatricesInput obdmi(node);
  }

  using invalid_input = testing::InvalidOneBodyDensityMatricesInput;
  for (auto input_xml : invalid_input::xml)
  {
    Libxml2Document doc;
    bool okay = doc.parseFromString(input_xml);
    REQUIRE(okay);
    xmlNodePtr node = doc.getRoot();

    CHECK_THROWS_AS(OneBodyDensityMatricesInput(node), UniformCommunicateError);
  }
}

TEST_CASE("OneBodyDensityMatricesInput::copy_construction", "[estimators]")
{
  using Input = testing::ValidOneBodyDensityMatricesInput;
  Libxml2Document doc;
  bool okay       = doc.parseFromString(Input::getXml(Input::valid::SCALE));
  xmlNodePtr node = doc.getRoot();
  OneBodyDensityMatricesInput obdmi(node);
  static_assert(std::is_copy_constructible_v<OneBodyDensityMatricesInput>);
}

} // namespace qmcplusplus
