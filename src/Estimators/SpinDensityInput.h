//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2022 QMCPACK developers.
//
// File developed by: Peter Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// Some code refactored from: SpinDensity.h
//////////////////////////////////////////////////////////////////////////////////////
#ifndef QMCPLUSPLUS_SPINDENSITYINPUT_H
#define QMCPLUSPLUS_SPINDENSITYINPUT_H

#include "Configuration.h"
#include "OhmmsData/ParameterSet.h"
#include "Containers/OhmmsPETE/TinyVector.h"

namespace qmcplusplus
{

class SpinDensityNew;

/** Native representation for Spin Density Estimators inputs
 *
 *  This class servers three purposes all related to properly handling
 *  and verifying the spin density input.
 *  1. An immutable representation of actual user input
 *  2. Parse the xml node of SpinDensityNew input.
 *  3. Hold the logic of calculating derived parameters.
 *
 */
class SpinDensityInput
{
public:
  static constexpr std::string_view type_tag{"SpinDensity"};
  using Real               = QMCTraits::RealType;
  using PosType            = QMCTraits::PosType;
  using Consumer           = SpinDensityNew;
  static constexpr int DIM = QMCTraits::DIM;

  SpinDensityInput(xmlNodePtr node);
  /** default copy constructor
   *  This is required due to SDI being part of a variant used as a vector element.
   */
  SpinDensityInput(const SpinDensityInput&) = default;
  Lattice get_cell() const { return cell_; }
  PosType get_corner() const { return corner_; }
  TinyVector<int, DIM> get_grid() const { return grid_; }
  int get_npoints() const { return npoints_; }
  bool get_write_report() const { return write_report_; }
  bool get_save_memory() const { return save_memory_; }
  const std::string& get_name() const { return name_; }
  const std::string& get_type() const { return type_; }

  struct DerivedParameters
  {
    PosType corner;
    TinyVector<int, DIM> grid;
    TinyVector<int, DIM> gdims;
    size_t npoints;
  };

  /** Derived parameters of SpinDensity
   *
   *  These require the cell the SpinDensity is evaluated over,
   *  the caller (SpinDensityNew) either gets this from the input and
   *  passes it back or passes in the cell from the relevant ParticleSet.
   *
   */
  DerivedParameters calculateDerivedParameters(const Lattice& lattice) const;

private:
  void readXML(xmlNodePtr cur);

  ///name of this Estimator
  std::string name_{type_tag};
  std::string type_{type_tag};

  Lattice cell_;
  PosType corner_;
  PosType dr_;
  PosType center_;
  TinyVector<int, DIM> grid_;
  int npoints_;
  bool write_report_;
  bool save_memory_;
  /** these are necessary for calculateDerivedParameters
   *
   *  If we are going to later write out a canonical input for
   *  this input then they are needed as well.
   */
  bool have_dr_     = false;
  bool have_grid_   = false;
  bool have_center_ = false;
  bool have_corner_ = false;
  bool have_cell_   = false;
};

} // namespace qmcplusplus
#endif /* SPINDENSITYINPUT_H */
