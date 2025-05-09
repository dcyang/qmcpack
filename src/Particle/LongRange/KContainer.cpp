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


#include "KContainer.h"
#include <map>
#include <cstdint>
#include "Message/Communicate.h"
#include "LRCoulombSingleton.h"
#include "Utilities/qmc_common.h"

namespace qmcplusplus
{

template<typename REAL>
const std::vector<typename KContainerT<REAL>::AppPosition>& KContainerT<
    REAL>::getKptsCartWorking() const
{
  // This is an `if constexpr` so it should not cost a branch at runtime.
  if constexpr (std::is_same_v<decltype(kpts_cart_), decltype(kpts_cart_working_)>)
    return kpts_cart_;
  else
    return kpts_cart_working_;
}

template<typename REAL>
const std::vector<REAL>& KContainerT<REAL>::getKSQWorking() const
{
  // This is an `if constexpr` so it should not cost a branch at runtime.
  if constexpr (std::is_same<decltype(ksq_), decltype(ksq_working_)>::value)
    return ksq_;
  else
    return ksq_working_;
}

template<typename REAL>
int KContainerT<REAL>::getMinusK(int k) const
{
  assert(k < minusk.size());
  return minusk[k];
}

template<typename REAL>
void KContainerT<REAL>::updateKLists(const Lattice& lattice,
                                     FullPrecReal kc,
                                     unsigned ndim,
                                     const Position& twist,
                                     bool useSphere)
{
  kcutoff = kc;
  if (kcutoff <= 0.0)
  {
    APP_ABORT("  Illegal cutoff for KContainer");
  }
  findApproxMMax(lattice, ndim);
  BuildKLists(lattice, twist, useSphere);

  app_log() << "  KContainer initialised with cutoff " << kcutoff << std::endl;
  app_log() << "   # of K-shell  = " << kshell.size() << std::endl;
  app_log() << "   # of K points = " << kpts_.size() << std::endl;
  app_log() << std::endl;
}

template<typename REAL>
void KContainerT<REAL>::findApproxMMax(const Lattice& lattice, unsigned ndim)
{
  //Estimate the size of the parallelpiped that encompasses a sphere of kcutoff.
  //mmax is stored as integer translations of the reciprocal cell vectors.
  //Does not require an orthorhombic cell.
  /* Old method.
  //2pi is not included in lattice.b
  Matrix<RealType> mmat;
  mmat.resize(3,3);
  for(int j=0;j<3;j++)
    for(int i=0;i<3;i++){
      mmat[i][j] = 0.0;
      for(int k=0;k<3;k++)
  mmat[i][j] = mmat[i][j] + 4.0*M_PI*M_PI*lattice.b(k)[i]*lattice.b(j)[k];
    }

  TinyVector<RealType,3> x,temp;
  RealType tempr;
  for(int idim=0;idim<3;idim++){
    int i = ((idim)%3);
    int j = ((idim+1)%3);
    int k = ((idim+2)%3);

    x[i] = 1.0;
    x[j] = (mmat[j][k]*mmat[k][i] - mmat[k][k]*mmat[i][j]);
    x[j]/= (mmat[j][j]*mmat[k][k] - mmat[j][k]*mmat[j][k]);
    x[k] = -(mmat[k][i] + mmat[j][k]*x[j])/mmat[k][k];

    for(i=0;i<3;i++){
      temp[i] = 0.0;
  for(j=0;j<3;j++)
    temp[i] += mmat[i][j]*x[j];
    }

    tempr = dot(x,temp);
    mmax[idim] = static_cast<int>(sqrt(4.0*kcut2/tempr)) + 1;
  }
  */
  // see rmm, Electronic Structure, p. 85 for details
  for (int i = 0; i < DIM; i++)
    mmax[i] = static_cast<int>(std::floor(std::sqrt(dot(lattice.a(i), lattice.a(i))) * kcutoff / (2 * M_PI))) + 1;

  mmax[DIM] = mmax[0];
  for (int i = 1; i < DIM; ++i)
    mmax[DIM] = std::max(mmax[i], mmax[DIM]);

  //overwrite the non-periodic directon to be zero
  if (LRCoulombSingleton::isQuasi2D())
  {
    app_log() << "  No kspace sum perpendicular to slab " << std::endl;
    mmax[2] = 0;
  }
  if (ndim < 3)
  {
    app_log() << "  No kspace sum along z " << std::endl;
    mmax[2] = 0;
  }
  if (ndim < 2)
    mmax[1] = 0;
}

template<typename REAL>
void KContainerT<REAL>::BuildKLists(const Lattice& lattice,
                                    const Position& twist,
                                    bool useSphere)
{
  TinyVector<int, DIM + 1> TempActualMax;
  TinyVector<int, DIM> kvec;
  TinyVector<FullPrecReal, DIM> kvec_cart;
  FullPrecReal modk2;
  std::vector<TinyVector<int, DIM>> kpts_tmp;
  std::vector<PositionFull> kpts_cart_tmp;
  std::vector<FullPrecReal> ksq_tmp;
  // reserve the space for memory efficiency
  if (useSphere)
  {
    const FullPrecReal kcut2 = kcutoff * kcutoff;
    //Loop over guesses for valid k-points.
    for (int i = -mmax[0]; i <= mmax[0]; i++)
    {
      kvec[0] = i;
      for (int j = -mmax[1]; j <= mmax[1]; j++)
      {
        kvec[1] = j;
        for (int k = -mmax[2]; k <= mmax[2]; k++)
        {
          kvec[2] = k;
          //Do not include k=0 in evaluations.
          if (i == 0 && j == 0 && k == 0)
            continue;
          //Convert kvec to Cartesian
          kvec_cart = lattice.k_cart(kvec + twist);
          //Find modk
          modk2 = dot(kvec_cart, kvec_cart);
          if (modk2 > kcut2)
            continue; //Inside cutoff?
          //This k-point should be added to the list
          kpts_tmp.push_back(kvec);
          kpts_cart_tmp.push_back(kvec_cart);
          ksq_tmp.push_back(modk2);
          //Update record of the allowed maximum translation.
          for (int idim = 0; idim < 3; idim++)
            if (std::abs(kvec[idim]) > TempActualMax[idim])
              TempActualMax[idim] = std::abs(kvec[idim]);
        }
      }
    }
  }
  else
  {
    // Loop over all k-points in the parallelpiped and add them to kcontainer
    // note layout is for interfacing with fft, so for each dimension, the
    // positive indexes come first then the negative indexes backwards
    // e.g.    0, 1, .... mmax, -mmax+1, -mmax+2, ... -1
    const int idimsize = mmax[0] * 2;
    const int jdimsize = mmax[1] * 2;
    const int kdimsize = mmax[2] * 2;
    for (int i = 0; i < idimsize; i++)
    {
      kvec[0] = i;
      if (kvec[0] > mmax[0])
        kvec[0] -= idimsize;
      for (int j = 0; j < jdimsize; j++)
      {
        kvec[1] = j;
        if (kvec[1] > mmax[1])
          kvec[1] -= jdimsize;
        for (int k = 0; k < kdimsize; k++)
        {
          kvec[2] = k;
          if (kvec[2] > mmax[2])
            kvec[2] -= kdimsize;
          // get cartesian location and modk2
          kvec_cart = lattice.k_cart(kvec);
          modk2     = dot(kvec_cart, kvec_cart);
          // add k-point to lists
          kpts_tmp.push_back(kvec);
          kpts_cart_tmp.push_back(kvec_cart);
          ksq_tmp.push_back(modk2);
        }
      }
    }
    // set allowed maximum translation
    TempActualMax[0] = mmax[0];
    TempActualMax[1] = mmax[1];
    TempActualMax[2] = mmax[2];
  }

  //Update a record of the number of k vectors
  numk = kpts_tmp.size();
  std::map<int64_t, std::vector<int>*> kpts_sorted;
  //create the map: use simple integer with resolution of 0.00000001 in ksq
  for (int ik = 0; ik < numk; ik++)
  {
    //This is a workaround for ewald bug (Issue #2105).  Basically, 1e-7 is the resolution of |k|^2 for doubles,
    //so we jack up the tolerance to match that.
    const int64_t k_ind = static_cast<int64_t>(ksq_tmp[ik] * 10000000);
    auto it(kpts_sorted.find(k_ind));
    if (it == kpts_sorted.end())
    {
      std::vector<int>* newSet = new std::vector<int>;
      kpts_sorted[k_ind]       = newSet;
      newSet->push_back(ik);
    }
    else
    {
      (*it).second->push_back(ik);
    }
  }
  std::map<int64_t, std::vector<int>*>::iterator it(kpts_sorted.begin());
  kpts_.resize(numk);
  kpts_cart_.resize(numk);
  kpts_cart_soa_.resize(numk);
  ksq_.resize(numk);
  kshell.resize(kpts_sorted.size() + 1, 0);
  int ok = 0, ish = 0;
  while (it != kpts_sorted.end())
  {
    std::vector<int>::iterator vit((*it).second->begin());
    while (vit != (*it).second->end())
    {
      int ik             = (*vit);
      kpts_[ok]          = kpts_tmp[ik];
      kpts_cart_[ok]     = kpts_cart_tmp[ik];
      kpts_cart_soa_(ok) = kpts_cart_tmp[ik];
      ksq_[ok]           = ksq_tmp[ik];
      ++vit;
      ++ok;
    }
    kshell[ish + 1] = kshell[ish] + (*it).second->size();
    ++it;
    ++ish;
  }
  kpts_cart_soa_.updateTo();
  if constexpr (!std::is_same<Real, FullPrecReal>::value)
  {
    // This copy implicity does the precision reduction.
    // the working vectors are not used or initialized for full precision builds.
    std::copy(kpts_cart_.begin(), kpts_cart_.end(), std::back_inserter(kpts_cart_working_));
    std::copy(ksq_.begin(), ksq_.end(), std::back_inserter(ksq_working_));
  }
  it = kpts_sorted.begin();
  std::map<int64_t, std::vector<int>*>::iterator e_it(kpts_sorted.end());
  while (it != e_it)
  {
    delete it->second;
    it++;
  }
  //Finished searching k-points. Copy list of maximum translations.
  mmax[DIM] = 0;
  for (int idim = 0; idim < DIM; idim++)
  {
    mmax[idim] = TempActualMax[idim];
    mmax[DIM]  = std::max(mmax[idim], mmax[DIM]);
    //if(mmax[idim] > mmax[DIM]) mmax[DIM] = mmax[idim];
  }
  //Now fill the array that returns the index of -k when given the index of k.
  minusk.resize(numk);

  //Assigns a unique hash value to each kpoint.
  auto getHashOfVec = [](const auto& inpv, int hashparam) -> int64_t {
    int64_t hash = 0; // this will cause integral promotion below
    for (int i = 0; i < inpv.Size; ++i)
      hash += inpv[i] + hash * hashparam;
    return hash;
  };

  // Create a map from the hash value for each k vector to the index
  std::map<int64_t, int> hashToIndex;
  for (int ki = 0; ki < numk; ki++)
  {
    hashToIndex[getHashOfVec(kpts_[ki], numk)] = ki;
  }
  // Use the map to find the index of -k from the index of k
  for (int ki = 0; ki < numk; ki++)
  {
    minusk[ki] = hashToIndex[getHashOfVec(-1 * kpts_[ki], numk)];
  }
}

#ifdef MIXED_PRECISION
template class KContainerT<float>;
template class KContainerT<double>;
#else
template class KContainerT<double>;
#endif
} // namespace qmcplusplus
