//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2024 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//
// File created by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Lab
//////////////////////////////////////////////////////////////////////////////////////

#ifndef QMCPLUSPLUS_VALID_STRUCTURE_FACTOR_INPUT_H
#define QMCPLUSPLUS_VALID_STRUCTURE_FACTOR_INPUT_H

namespace qmcplusplus
{
namespace testing
{

class ValidStructureFactorInput
{
  /** The first will produce the legacy sk estimator behavior
   *  The second will produce the skall behavior
   */
  static constexpr std::array<std::string_view, 3> xml{
      R"XML(
<estimator type="StructureFactor" name="sk1" hdf5="yes"/>
)XML",
      R"XML(
<estimator type="StructureFactor" name="sk1" source="ion" target="e" writerho="yes" hdf5="yes"/>
)XML",
      R"XML(
<estimator type="StructureFactor" name="sk1" source="ion0" target="e" hdf5="yes" writerho="yes" writeionion="yes"/>
)XML"};

public:
  enum class valid
  {
    SK = 0,
    SKALL,
    SKALL_IONION
  };
  static std::string_view getXml(valid val) { return xml[static_cast<std::size_t>(val)]; }
  auto begin() { return xml.begin(); }
  auto end() { return xml.end(); }
};

class InvalidStructureFactorInput
{
  /** Invalid input sections
   *  writerho without explicit source and target
   *  writeionion without explicit source and target
   */
  static constexpr std::array<std::string_view, 2> xml{
      R"XML(
<estimator type="StructureFactor" name="sk1" writerho="yes" hdf5="yes"/>
)XML",
      R"XML(
<estimator type="StructureFactor" name="sk1" hdf5="yes" writeionion="yes"/>
)XML"};

public:
  auto begin() { return xml.begin(); }
  auto end() { return xml.end(); }
};

} // namespace testing
} // namespace qmcplusplus
#endif
