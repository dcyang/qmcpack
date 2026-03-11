//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2023 QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Jeongnim Kim, jeongnim.kim@inte.com, Intel Corp.
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#include "SplineSetReader.h"
#include "OneSplineOrbData.hpp"
#include "Utilities/FairDivide.h"
#include "einspline_helper.hpp"
#include <Timer.h>
#if defined(QMC_COMPLEX)
#include "SplineC2C.h"
#include "SplineC2COMPTarget.h"
#else
#include "SplineR2R.h"
#include "SplineC2R.h"
#include "SplineC2ROMPTarget.h"
#endif
#include "Message/CommOperators.h"
#include "spline2/SplineUtils.h"
#include "spline2/MultiBsplineBase.hpp"


namespace qmcplusplus
{
template<typename ST>
SplineSetReader<ST>::SplineSetReader(EinsplineSetBuilder* e, bool use_duplex_splines)
    : BsplineReader(e, use_duplex_splines)
{}

template<typename ST>
std::unique_ptr<SPOSet> SplineSetReader<ST>::create_spline_set(const std::string& my_name,
                                                               int spin,
                                                               const BandInfoGroup& bandgroup)
{
  const int N = bandgroup.getNumDistinctOrbitals();

  if (use_duplex_splines_)
    app_log() << "  Using complex einspline table" << std::endl;
  else
    app_log() << "  Using real einspline table" << std::endl;
  const TinyVector<int, 3> half_g = use_duplex_splines_
      ? TinyVector<int, 3>(0, 0, 0)
      : computeHalfG(mybuilder->TargetPtcl.getLattice().BoxBConds, mybuilder->primcell_kpoints,
                     bandgroup.myBands[0].TwistIndex);

  Ugrid xyz_grid[3];
  typename bspline_traits<ST, 3>::BCType xyz_bc[3];
  set_grid(mybuilder->MeshSize, half_g, xyz_grid, xyz_bc);

  const size_t num_splines = getAlignedSize<ST>(use_duplex_splines_ ? N * 2 : N);
  std::unique_ptr<MultiBsplineBase<ST>> multi_splines_ptr;
  if (use_offload)
    multi_splines_ptr = std::make_unique<MultiBsplineOffload<ST>>(xyz_grid, xyz_bc, num_splines);
  else
    multi_splines_ptr = std::make_unique<MultiBspline<ST>>(xyz_grid, xyz_bc, num_splines);

  auto& multi_splines(*multi_splines_ptr);
  app_log() << "MEMORY " << multi_splines.sizeInByte() / (1 << 20) << " MB allocated "
            << "for the coefficients in 3D spline orbital representation" << std::endl;

  std::unique_ptr<BsplineSet> bspline;
#if defined(QMC_COMPLEX)
  if (use_offload)
    bspline = std::make_unique<SplineC2COMPTarget<ST>>(my_name, bandgroup.getNumSPOs(), mybuilder->PrimCell,
                                                       std::move(multi_splines_ptr), use_offload);
  else
    bspline = std::make_unique<SplineC2C<ST>>(my_name, bandgroup.getNumSPOs(), mybuilder->PrimCell,
                                              std::move(multi_splines_ptr), use_offload);
#else
  if (use_duplex_splines_)
  {
    if (use_offload)
      bspline = std::make_unique<SplineC2ROMPTarget<ST>>(my_name, bandgroup.getNumSPOs(), mybuilder->PrimCell,
                                                         std::move(multi_splines_ptr), use_offload);
    else
      bspline = std::make_unique<SplineC2R<ST>>(my_name, bandgroup.getNumSPOs(), mybuilder->PrimCell,
                                                std::move(multi_splines_ptr), use_offload);
  }
  else
    bspline = std::make_unique<SplineR2R<ST>>(my_name, bandgroup.getNumSPOs(), mybuilder->PrimCell,
                                              std::move(multi_splines_ptr), use_offload);
#endif

  app_log() << "  ClassName = " << bspline->getClassName() << std::endl;

  bspline->resizeStorage(N);

  bspline->HalfG = half_g;

  if (use_duplex_splines_)
  {
    //baseclass handles twists
    check_twists(*bspline, bandgroup);
  }

  bool foundspline = lookforSplineDataDumpFile(bandgroup, bspline->getKeyword(), sizeof(ST));
  if (foundspline && myComm->rank() == 0)
  {
    Timer now;
    hdf_archive h5f(myComm);
    const auto splinefile = getSplineDumpFileName(bandgroup);
    h5f.open(splinefile, H5F_ACC_RDONLY);
    foundspline = SplineUtils<ST>::read(multi_splines, h5f);
    if (foundspline)
      app_log() << "  Successfully restored 3D B-spline coefficients from " << splinefile << ". The reading time is "
                << now.elapsed() << " sec." << std::endl;
  }

  if (!foundspline)
  {
    multi_splines.flush_zero();

    Timer now;
    initialize_spline_pio_gather(spin, bandgroup, half_g, bspline->BandIndexMap, multi_splines);
    app_log() << "  SplineSetReader initialize_spline_pio " << now.elapsed() << " sec" << std::endl;

    if (saveSplineCoefs && myComm->rank() == 0)
    {
      Timer now;
      const std::string splinefile(getSplineDumpFileName(bandgroup));
      hdf_archive h5f;
      h5f.create(splinefile);
      std::string classname = bspline->getClassName();
      h5f.write(classname, "class_name");
      int sizeD = sizeof(ST);
      h5f.write(sizeD, "sizeof");
      SplineUtils<ST>::write(multi_splines, h5f);
      h5f.close();
      app_log() << "  Stored spline coefficients in " << splinefile << " for potential reuse. The writing time is "
                << now.elapsed() << " sec." << std::endl;
    }
  }

  {
    myComm->barrier();
    Timer now;
    SplineUtils<ST>::bcast(multi_splines, *myComm);
    app_log() << "  Time to bcast the table = " << now.elapsed() << std::endl;
  }

  return bspline;
}

template<typename ST>
bool SplineSetReader<ST>::lookforSplineDataDumpFile(const BandInfoGroup& bandgroup,
                                                    const std::string& keyword,
                                                    size_t datatype_size) const
{
  int foundspline = 0;
  Timer now;
  if (myComm->rank() == 0)
  {
    now.restart();
    hdf_archive h5f(myComm);
    foundspline = h5f.open(getSplineDumpFileName(bandgroup), H5F_ACC_RDONLY);
    if (foundspline)
    {
      std::string aname("none");
      foundspline = h5f.readEntry(aname, "class_name");
      foundspline = (aname.find(keyword) != std::string::npos);
    }
    if (foundspline)
    {
      int sizeD   = 0;
      foundspline = h5f.readEntry(sizeD, "sizeof");
      foundspline = (sizeD == datatype_size);
    }
    h5f.close();
  }
  myComm->bcast(foundspline);
  return foundspline;
}

template<typename ST>
void SplineSetReader<ST>::readOneOrbitalCoefs(const std::string& s,
                                              hdf_archive& h5f,
                                              Vector<std::complex<double>>& cG) const
{
  if (!h5f.readEntry(cG, s))
  {
    std::ostringstream msg;
    msg << "SplineSetReader Failed to read band(s) from h5 file. " << "Attempted dataset " << s << " with " << cG.size()
        << " complex numbers." << std::endl;
    throw std::runtime_error(msg.str());
  }
  double total_norm = compute_norm(cG);
  if ((checkNorm) && (std::abs(total_norm - 1.0) > PW_COEFF_NORM_TOLERANCE))
  {
    std::ostringstream msg;
    msg << "SplineSetReader The orbital dataset " << s << " has a wrong norm " << total_norm
        << ", computed from plane wave coefficients!" << std::endl
        << "This may indicate a problem with the HDF5 library versions used "
        << "during wavefunction conversion or read." << std::endl;
    throw std::runtime_error(msg.str());
  }
}

template<typename ST>
void SplineSetReader<ST>::initialize_spline_pio_gather(const int spin,
                                                       const BandInfoGroup& bandgroup,
                                                       const TinyVector<int, 3>& half_g,
                                                       const aligned_vector<int>& band_index_map,
                                                       MultiBsplineBase<ST>& multi_splines) const
{
  //distribute bands over processor groups
  int Nbands            = bandgroup.getNumDistinctOrbitals();
  const int Nprocs      = myComm->size();
  const int Nbandgroups = std::min(Nbands, Nprocs);
  Communicate band_group_comm(*myComm, Nbandgroups);
  std::vector<int> band_groups(Nbandgroups + 1, 0);
  FairDivideLow(Nbands, Nbandgroups, band_groups);
  int iorb_first = band_groups[band_group_comm.getGroupID()];
  int iorb_last  = band_groups[band_group_comm.getGroupID() + 1];

  app_log() << "Start transforming plane waves to 3D B-Splines." << std::endl;
  OneSplineOrbData oneband(mybuilder->MeshSize, half_g, use_duplex_splines_);
  hdf_archive h5f(&band_group_comm, false);
  Vector<std::complex<double>> cG(mybuilder->Gvecs[0].size());
  const std::vector<BandInfo>& cur_bands = bandgroup.myBands;
  if (band_group_comm.isGroupLeader())
  {
    h5f.open(mybuilder->H5FileName, H5F_ACC_RDONLY);
    for (int iorb = iorb_first; iorb < iorb_last; iorb++)
    {
      const auto& cur_band = cur_bands[band_index_map[iorb]];
      const int ti         = cur_band.TwistIndex;
      readOneOrbitalCoefs(psi_g_path(ti, spin, cur_band.BandIndex), h5f, cG);
      oneband.fft_spline(cG, mybuilder->Gvecs[0], mybuilder->primcell_kpoints[ti], rotate);
      if (use_duplex_splines_)
      {
        copy_spline<double, ST>(oneband.get_spline_r(), *multi_splines.getSplinePtr(), iorb * 2);
        copy_spline<double, ST>(oneband.get_spline_i(), *multi_splines.getSplinePtr(), iorb * 2 + 1);
      }
      else
        copy_spline<double, ST>(oneband.get_spline_r(), *multi_splines.getSplinePtr(), iorb);
    }

    {
      band_group_comm.getGroupLeaderComm()->barrier();
      Timer now;
      if (use_duplex_splines_)
      {
        std::vector<int> offset(band_groups.size());
        for (int i = 0; i < offset.size(); i++)
          offset[i] = band_groups[i] * 2;
        SplineUtils<ST>::gatherv(multi_splines, offset, *band_group_comm.getGroupLeaderComm());
      }
      else
        SplineUtils<ST>::gatherv(multi_splines, band_groups, *band_group_comm.getGroupLeaderComm());
      app_log() << "  Time to gather the table = " << now.elapsed() << std::endl;
    }
  }
}

template class SplineSetReader<float>;
#if !defined(QMC_MIXED_PRECISION)
template class SplineSetReader<double>;
#endif
} // namespace qmcplusplus
