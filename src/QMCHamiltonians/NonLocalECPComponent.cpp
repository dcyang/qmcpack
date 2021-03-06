//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jaron T. Krogel, krogeljt@ornl.gov, Oak Ridge National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////
    
    
#include "Particle/DistanceTableData.h"
#include "QMCHamiltonians/NonLocalECPComponent.h"

namespace qmcplusplus
{

NonLocalECPComponent::NonLocalECPComponent():
  lmax(0), nchannel(0), nknot(0), Rmax(-1), VP(0)
{
#if !defined(REMOVE_TRACEMANAGER)
  streaming_particles = false;
#endif
}

NonLocalECPComponent::~NonLocalECPComponent()
{
  for(int ip=0; ip<nlpp_m.size(); ip++)
    delete nlpp_m[ip];
  if(VP) delete VP;
}

NonLocalECPComponent* NonLocalECPComponent::makeClone(const ParticleSet &qp)
{
  NonLocalECPComponent* myclone=new NonLocalECPComponent(*this);
  for(int i=0; i<nlpp_m.size(); ++i)
    myclone->nlpp_m[i]=nlpp_m[i]->makeClone();
  if(VP) myclone->VP=new VirtualParticleSet(qp,nknot);
  return myclone;
}

void NonLocalECPComponent::initVirtualParticle(const ParticleSet &qp)
{
  assert(VP==0);
  VP=new VirtualParticleSet(qp,nknot);
}

void NonLocalECPComponent::add(int l, RadialPotentialType* pp)
{
  angpp_m.push_back(l);
  wgt_angpp_m.push_back(static_cast<RealType>(2*l+1));
  nlpp_m.push_back(pp);
}

void NonLocalECPComponent::resize_warrays(int n,int m,int l)
{
  psiratio.resize(n);
  psigrad.resize(n);
  psigrad_source.resize(n);
  vrad.resize(m);
  dvrad.resize(m);
  wvec.resize(m);
  Amat.resize(n*m);
  dAmat.resize(n*m);
  lpol.resize(l+1,1.0);
  dlpol.resize(l+1,0.0);
  rrotsgrid_m.resize(n);
  nchannel=nlpp_m.size();
  nknot=sgridxyz_m.size();
  //This is just to check
  //for(int nl=1; nl<nlpp_m.size(); nl++) nlpp_m[nl]->setGridManager(false);
  if(lmax)
  {
    if(lmax>7) 
    {
      APP_ABORT("Increase the maximum angular momentum implemented.");
    }
    //Lfactor1.resize(lmax);
    //Lfactor2.resize(lmax);
    for(int nl=0; nl<lmax; nl++)
    {
      Lfactor1[nl]=static_cast<RealType>(2*nl+1);
      Lfactor2[nl]=1.0e0/static_cast<RealType>(nl+1);
    }
  }
}

void NonLocalECPComponent::print(std::ostream& os)
{
  os << "    Maximum angular mementum = "<<  lmax << std::endl;
  os << "    Number of non-local channels = " << nchannel << std::endl;
  for(int l=0; l <nchannel; l++)
    os << "       l(" << l << ")=" << angpp_m[l] << std::endl;
  os << "    Cutoff radius = " << Rmax << std::endl;
  os << "    Spherical grids and weights: " << std::endl;
  for(int ik=0; ik<nknot; ik++)
    os << "       " << sgridxyz_m[ik] << std::setw(20) << sgridweight_m[ik] << std::endl;
}

NonLocalECPComponent::RealType
NonLocalECPComponent::evaluateOne(ParticleSet& W, int iat, TrialWaveFunction& psi, 
    int iel, RealType r, const PosType& dr, 
    bool Tmove, std::vector<NonLocalData>& Txy) const
{
  constexpr RealType czero(0);
  constexpr RealType cone(1);

  RealType lpol_[lmax+1];
  RealType vrad_[nchannel];
  std::vector<RealType> psiratio_(nknot);
  PosType deltaV[nknot];

  if(VP)
  {
    // Compute ratios with VP
    ParticleSet::ParticlePos_t VPos(nknot);
    for (int j=0; j<nknot; j++)
    {
      deltaV[j]=r*rrotsgrid_m[j]-dr;
      VPos[j]=deltaV[j]+W.R[iel];
    }
    VP->makeMoves(iel,VPos,true,iat);
    psi.evaluateRatios(*VP,psiratio_);
    for (int j=0; j<nknot; j++)
      psiratio_[j]*=sgridweight_m[j];
  }
  else
  {
    // Compute ratio of wave functions
    for (int j=0; j<nknot; j++)
    {
      deltaV[j]=r*rrotsgrid_m[j]-dr;
      W.makeMoveOnSphere(iel,deltaV[j]);
#if defined(QMC_COMPLEX)
      psiratio_[j]=psi.ratio(W,iel)*sgridweight_m[j]*std::cos(psi.getPhaseDiff());
#else
      psiratio_[j]=psi.ratio(W,iel)*sgridweight_m[j];
#endif
      W.rejectMove(iel);
      psi.resetPhaseDiff();
      //psi.rejectMove(iel);
    }
  }

  // Compute radial potential, multiplied by (2l+1) factor.
  for(int ip=0; ip< nchannel; ip++)
    vrad_[ip]=nlpp_m[ip]->splint(r)*wgt_angpp_m[ip];

  const RealType rinv=cone/r;
  RealType pairpot=0; 
  // Compute spherical harmonics on grid
  for (int j=0; j<nknot ; j++)
  {
    RealType zz=dot(dr,rrotsgrid_m[j])*rinv;
    // Forming the Legendre polynomials
    lpol_[0]=cone;
    RealType lpolprev=czero;
    for (int l=0 ; l< lmax ; l++)
    {
      //Not a big difference
      //lpol[l+1]=(2*l+1)*zz*lpol[l]-l*lpolprev;
      //lpol[l+1]/=(l+1);
      lpol_[l+1]=Lfactor1[l]*zz*lpol_[l]-l*lpolprev;
      lpol_[l+1]*=Lfactor2[l];
      lpolprev=lpol_[l];
    }

    RealType lsum=czero;
    for(int l=0; l <nchannel; l++)
      lsum += vrad_[l]*lpol_[ angpp_m[l] ];
    lsum *= psiratio_[j];
    if(Tmove) Txy.push_back(NonLocalData(iel,lsum,deltaV[j]));
    pairpot+=lsum;
  }

#if !defined(REMOVE_TRACEMANAGER)
  if( streaming_particles)
  {
    (*Vi_sample)(iat) += .5*pairpot;
    (*Ve_sample)(iel) += .5*pairpot;
  }
#endif
  return pairpot;
}

NonLocalECPComponent::RealType
NonLocalECPComponent::evaluateOneWithForces(ParticleSet& W, int iat, TrialWaveFunction& psi, 
    int iel, RealType r, const PosType& dr, 
    PosType & force_iat, bool Tmove, std::vector<NonLocalData>& Txy) const
{
  constexpr RealType czero(0);
  constexpr RealType cone(1);

  //Array for P_l[cos(theta)].
  RealType lpol_[lmax+1];
  //Array for P'_l[cos(theta)]
  RealType dlpol_[lmax+1];

  //Array for v_l(r).
  RealType vrad_[nchannel];
  //Array for (2l+1)*v'_l(r)/r.
  RealType dvrad_[nchannel];

  //$\Psi(...q...)/\Psi(...r...)$ for all quadrature points q.
  std::vector<RealType> psiratio_(nknot);
  //$\nabla \Psi(...q...)/\Psi(...r...)$ for all quadrature points q.
  //  $\nabla$ is w.r.t. the electron coordinates involved in the quadrature.
  std::vector<PosType> gradpsiratio_(nknot);

  //This stores gradient of v(r):
  std::vector<PosType> vgrad_(nchannel);  
  //This stores the gradient of the cos(theta) term in force expression.
  std::vector<PosType> cosgrad_(nknot);
  //This stores grad psi/psi - dot(u,grad psi)
  std::vector<PosType> wfngrad_(nknot);

  //^^^ "Why not GradType?" you ask.  Because GradType can be complex.
  PosType deltaV[nknot];

  GradType gradtmp_(0);
  PosType realgradtmp_(0);

  //Pseudopotential derivative w.r.t. ions can be split up into 3 contributions:
  // term coming from the gradient of the radial potential
  PosType gradpotterm_(0);
  // term coming from gradient of legendre polynomial
  PosType gradlpolyterm_(0);
  // term coming from dependence of quadrature grid on ion position.
  PosType gradwfnterm_(0);

  if(VP)
  {
    APP_ABORT("NonLocalECPComponent::evaluateOneWithForces(...): Forces not implemented with virtual particle moves\n");
    // Compute ratios with VP
    ParticleSet::ParticlePos_t VPos(nknot);
    for (int j=0; j<nknot; j++)
    {
      deltaV[j]=r*rrotsgrid_m[j]-dr;
      VPos[j]=deltaV[j]+W.R[iel];
    }
    VP->makeMoves(iel,VPos,true,iat);
    psi.evaluateRatios(*VP,psiratio_);
    for (int j=0; j<nknot; j++)
      psiratio_[j]*=sgridweight_m[j];
  }
  else
  { 
    ValueType ratio(0);
    // Compute ratio of wave functions
    for (int j=0; j<nknot; j++)
    {
      deltaV[j]=r*rrotsgrid_m[j]-dr;
      W.makeMoveOnSphere(iel,deltaV[j]);
#if defined(QMC_COMPLEX)
      gradtmp_=0;

      RealType ratio_mag=psi.ratioGrad(W,iel,gradtmp_);
      RealType ratio_r=ratio_mag*std::cos(psi.getPhaseDiff());
      RealType ratio_i=ratio_mag*std::sin(psi.getPhaseDiff());

      ratio=ValueType(ratio_r,ratio_i);
     
      //Only need the real part.
      psiratio_[j]=ratio_r*sgridweight_m[j]; 
     
      //QMCPACK spits out $\nabla\Psi(q)/\Psi(q)$.
      //Multiply times $\Psi(q)/\Psi(r)$ to get 
      // $\nabla\Psi(q)/\Psi(r)
      gradtmp_*=ratio;

      //And now we take the real part and save it.  
      convert(gradtmp_,gradpsiratio_[j]);
     
    
#else
      //Real nonlocalpp forces seem to differ from those in the complex build.  Since
      //complex build has been validated against QE, that indicates there's a bug for the real build.
      
      APP_ABORT("NonLocalECPComponent::evaluateOneWithForces(...): Forces not implemented for real wavefunctions yet.");

      gradtmp_=0;
      ratio=psi.ratioGrad(W,iel,gradtmp_);
      
      gradtmp_*=ratio;
      psiratio_[j]=ratio*sgridweight_m[j]; 
     
      gradpsiratio_[j]=gradtmp_;
#endif
      W.rejectMove(iel);
      psi.resetPhaseDiff();
      //psi.rejectMove(iel);
    }
  }

  // This is just a temporary variable to dump d2/dr2 into for spline evaluation.
  RealType secondderiv(0);

  const RealType rinv=cone/r;

  // Compute radial potential and its derivative times (2l+1)
  for(int ip=0; ip< nchannel; ip++)
  {
    //fun fact.  NLPComponent stores v(r) as v(r), and not as r*v(r) like in other places.  
    vrad_[ip]=nlpp_m[ip]->splint(r,dvrad_[ip],secondderiv)*wgt_angpp_m[ip];
    vgrad_[ip]=dvrad_[ip]*dr*wgt_angpp_m[ip]*rinv;
  }

  RealType pairpot=0; 
  // Compute spherical harmonics on grid
  for (int j=0;  j<nknot ; j++)
  {
    RealType zz=dot(dr,rrotsgrid_m[j])*rinv;
    PosType uminusrvec=rrotsgrid_m[j]-zz*dr*rinv;
    
    cosgrad_[j]=rinv*uminusrvec;

    RealType udotgradpsi=dot(gradpsiratio_[j],rrotsgrid_m[j]);
    wfngrad_[j]=gradpsiratio_[j] - dr*(udotgradpsi*rinv);
    wfngrad_[j]*=sgridweight_m[j];

    // Forming the Legendre polynomials
    //P_0(x)=1; P'_0(x)=0.
    lpol_[0]=cone;
    dlpol_[0]=czero;

    RealType lpolprev=czero;
    RealType dlpolprev=czero;

    for (int l=0 ; l< lmax ; l++)
    {
      //Legendre polynomial recursion formula.
      lpol_[l+1]=Lfactor1[l]*zz*lpol_[l]-l*lpolprev;
      lpol_[l+1]*=Lfactor2[l];
       
      //and for the derivative...
      dlpol_[l+1]=Lfactor1[l]*(zz*dlpol_[l]+lpol_[l])-l*dlpolprev;
      dlpol_[l+1]*=Lfactor2[l];

      lpolprev=lpol_[l];
      dlpolprev=dlpol_[l];
    }

    RealType lsum=czero;
   // Now to compute the forces:
    gradpotterm_  =0;
    gradlpolyterm_=0;
    gradwfnterm_  =0;
    
    for(int l=0; l<nchannel; l++)
    {
      lsum          += vrad_[l] * lpol_[ angpp_m[l] ] * psiratio_[j];
      gradpotterm_  += vgrad_[l] * lpol_[ angpp_m[l] ] * psiratio_[j];
      gradlpolyterm_ += vrad_[l] * dlpol_[ angpp_m[l] ] * cosgrad_[j] * psiratio_[j];
      gradwfnterm_  += vrad_[l] * lpol_[ angpp_m[l] ] * wfngrad_[j]; 
    }
    
    if(Tmove) Txy.push_back(NonLocalData(iel,lsum,deltaV[j]));
    pairpot+=lsum;
    force_iat+=  gradpotterm_ + gradlpolyterm_ - gradwfnterm_;
    
  }

#if !defined(REMOVE_TRACEMANAGER)
  if( streaming_particles)
  {
    (*Vi_sample)(iat) += .5*pairpot;
    (*Ve_sample)(iel) += .5*pairpot;
  }
#endif
  return pairpot;
}

///Randomly rotate sgrid_m
void NonLocalECPComponent::randomize_grid(RandomGenerator_t& myRNG)
{
  RealType phi(TWOPI*myRNG()), psi(TWOPI*myRNG()), cth(myRNG()-0.5);
  RealType sph(std::sin(phi)),cph(std::cos(phi)),
           sth(std::sqrt(1.0-cth*cth)),sps(std::sin(psi)),
           cps(std::cos(psi));
  TensorType rmat( cph*cth*cps-sph*sps, sph*cth*cps+cph*sps,-sth*cps,
                   -cph*cth*sps-sph*cps,-sph*cth*sps+cph*cps, sth*sps,
                   cph*sth,             sph*sth,             cth     );
  for(int i=0; i<sgridxyz_m.size(); i++)
    rrotsgrid_m[i] = dot(rmat,sgridxyz_m[i]);
}

template<typename T>
void NonLocalECPComponent::randomize_grid(std::vector<T> &sphere, RandomGenerator_t& myRNG)
{
  RealType phi(TWOPI*myRNG()), psi(TWOPI*myRNG()), cth(myRNG()-0.5);
  RealType sph(std::sin(phi)),cph(std::cos(phi)),
           sth(std::sqrt(1.0-cth*cth)),sps(std::sin(psi)),
           cps(std::cos(psi));
  TensorType rmat( cph*cth*cps-sph*sps, sph*cth*cps+cph*sps,-sth*cps,
                   -cph*cth*sps-sph*cps,-sph*cth*sps+cph*cps, sth*sps,
                   cph*sth,             sph*sth,             cth     );
  SpherGridType::iterator it(sgridxyz_m.begin());
  SpherGridType::iterator it_end(sgridxyz_m.end());
  SpherGridType::iterator jt(rrotsgrid_m.begin());
  while(it != it_end)
  {
    *jt = dot(rmat,*it);
    ++it;
    ++jt;
  }
  //copy the randomized grid to sphere
  for (int i=0; i<rrotsgrid_m.size(); i++)
    for (int j=0; j<OHMMS_DIM; j++)
      sphere[OHMMS_DIM*i+j] = rrotsgrid_m[i][j];
}

template void NonLocalECPComponent::randomize_grid(std::vector<float> &sphere, RandomGenerator_t& myRNG);
template void NonLocalECPComponent::randomize_grid(std::vector<double> &sphere, RandomGenerator_t& myRNG);


}
