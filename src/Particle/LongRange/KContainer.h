//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2025 QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Peter W. Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_KCONTAINER_H
#define QMCPLUSPLUS_KCONTAINER_H

#include "Configuration.h"
#include "OhmmsSoA/VectorSoaContainer.h"
#include "OMPTarget/OffloadAlignedAllocators.hpp"

namespace qmcplusplus
{
/** Container for k-points
 *
 * It generates a set of k-points that are unit-translations of the reciprocal-space
 * cell. K-points are generated within a spherical cutoff set by the supercell
 */
template<typename REAL>
class KContainerT
{
public:
  using Real               = REAL;
  using FullPrecReal       = QMCTraits::FullPrecRealType;
  using PositionFull       = QMCTraits::QTFull::PosType;
  static constexpr int DIM = OHMMS_DIM;
  using AppPosition        = typename QMCTraits::PosType;
  using Position           = typename QMCTypes<Real, DIM>::PosType;
  using Kpts               = std::vector<TinyVector<int, DIM>>;

private:
  /// The cutoff up to which k-vectors are generated.
  FullPrecReal kcutoff;

public:
  const auto& get_kpts_cart_soa() const { return kpts_cart_soa_; }
  const std::vector<TinyVector<int, DIM>>& getKpts() const { return kpts_; }

  /** @ingroup Working Precision accessors
   *  @brief   return reference kpt values in the working precision of build
   *           these encapsulate the cached reduced precision or pass through of the
   *           full precision values.
   *  @{
   */
  ///get cartesian kpt representation in the working precision of the application.
  const std::vector<AppPosition>& getKptsCartWorking() const;
  ///get ksqr in the working precision of the application.
  const std::vector<Real>& getKSQWorking() const;
  /// @}

  const std::vector<int>& getKShell() const { return kshell; };


  int getMinusK(int k) const;

  int getNumK() const { return numk; }

  /** update k-vectors
   * @param sc supercell
   * @param kc cutoff radius in the K
   * @param twist shifts the center of the grid of k-vectors
   * @param useSphere if true, use the |K|
   */
  void updateKLists(const Lattice& lattice,
                    FullPrecReal kc,
                    unsigned ndim,
                    const Position& twist = Position(),
                    bool useSphere        = true);

private:
  ///number of k-points
  int numk;

  /** maximum integer translations of reciprocal cell within kc.
   *
   * Last index is max. of first dimension+1
   */
  TinyVector<int, DIM + 1> mmax;

  /** k-vector in reduced coordinates
   */
  Kpts kpts_;
  /// k-vectors in Cartesian coordinates
  std::vector<Position> kpts_cart_;
  /// k squared at full precision
  std::vector<FullPrecReal> ksq_;
  /** @ingroup Cached Working Precision values
   *  @brief   only used or initialized for mixed precision
   *  @{
   */
  std::vector<AppPosition> kpts_cart_working_;
  std::vector<Real> ksq_working_;
  /// @}

  /** Given a k index, return index to -k
   */
  std::vector<int> minusk;
  /** kpts which belong to the ith-shell [kshell[i], kshell[i+1]) */
  std::vector<int> kshell;

  /** compute approximate parallelpiped that surrounds kc
   * @param lattice supercell
   */
  void findApproxMMax(const Lattice& lattice, unsigned ndim);
  /** construct the container for k-vectors */
  void BuildKLists(const Lattice& lattice, const Position& twist, bool useSphere);

  /** K-vector in Cartesian coordinates in SoA layout
   */
  VectorSoaContainer<FullPrecReal, DIM, OffloadAllocator<FullPrecReal>> kpts_cart_soa_;
};

using KContainer = KContainerT<QMCTraits::RealType>;

#ifdef MIXED_PRECISION
extern template class KContainerT<float>;
extern template class KContainerT<double>;
#else
extern template class KContainerT<double>;
#endif

} // namespace qmcplusplus

#endif
