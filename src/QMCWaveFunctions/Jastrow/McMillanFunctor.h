//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: ChangMo Yang, dcyang@unist.ac.kr, UNIST, Korea
//
// File created by: ChangMo Yang, dcyang@unist.ac.kr, UNIST, Korea
//////////////////////////////////////////////////////////////////////////////////////
    
    
/** @file McMillanFunctor.h
 * @brief Functor which implements McMillan two-body correlations
 */
#ifndef QMCPLUSPLUS_MCMILLAN1965_H
#define QMCPLUSPLUS_MCMILLAN1965_H
#include "Numerics/OptimizableFunctorBase.h"
#include "OhmmsData/AttributeSet.h"
#include <cmath>
// #include <vector>
#include "OhmmsPETE/TinyVector.h"


namespace qmcplusplus
{

  /** Implements a function \f$ u(r) = \frac{1}{2*A}\left[1-\exp(-A*r)\right], \f$
   *
   * McMillan, PR 138 A442 (1965)
   *
   */
  // XXX: old code below
  template<class T>
    struct McMillanFunctor: public OptimizableFunctorBase
  {
    ///whether the parameters are optimizable
    bool Opt_A, Opt_B;
    ///input A, B
    real_type A, B;

    /// internal variables
    real_type rcusp, c1, c2;

    ///id of A, B
    std::string ID_A, ID_B;

    ///default constructor
    explicit McMillanFunctor(real_type a = 5.0, real_type b = 5.7448, real_type r0 = 2.5)
      : A(a), B(b), rcusp(r0), Opt_A(true), Opt_B(true), ID_A("0"), ID_B("0")
    {
      reset();
    }

    OptimizableFunctorBase* makeClone() const
    {
      return new McMillanFunctor<T>(*this);
    }

    void reset() {
      reset(A, B);
      // rcusp = 2.5;
      real_type Y;
      c1 = evaluate(rcusp, Y, c2);
      c2 = -0.5*Y/(rcusp*c1);
      c1 *= std::exp(-c2*rcusp*rcusp);
    }
    inline void reset(real_type a, real_type b) { A = a; B = b; }

    inline real_type evaluate(real_type r) const
    {
      return ( r < rcusp ? c1*std::exp(-c2*r*r) :
          ( r < cutoff_radius ? std::pow(B,A)*(std::pow(r,-A) + std::pow(2.0*cutoff_radius-r,-A) - 2.0*std::         pow(cutoff_radius,-A)) :
            0.0 ) );
    }

    inline real_type
      evaluate(real_type r, real_type& dudr, real_type& d2udr2) const
      {
	real_type val;
	if (r < rcusp) {
	  val = c1*std::exp(-c2*r*r);
	  dudr = -2.0*c1*c2*r*std::exp(-c2*r*r);
	  d2udr2 = -2.0*c1*c2*std::exp(-c2*r*r) + 4.0*c1*c2*c2*r*r*std::exp(-c2*r*r);
	  /*
	     val = 17318.0912945689*std::exp(-1.00041061270143*r*r);
	     dudr = -34650.4046456380*r*std::exp(-1.00041061270143*r*r);
	     d2udr2 = -34650.4046456380*std::exp(-1.00041061270143*r*r) + 69329.2650837902*r*r*std::exp(-1.00041061270143*r*r);
	     */
	}
	else if (r < cutoff_radius) {
	  real_type rd = 2.0*cutoff_radius - r;
	  val = std::pow(B,A)*(std::pow(r,-A) + std::pow(rd,-A) - 2.0*std::pow(cutoff_radius,-A));
	  dudr = A*std::pow(B,A)*(std::pow(rd,-(A+1.0))-std::pow(r,-(A+1.0)));
	  d2udr2 = (A+1.0)*A*std::pow(B,A)*(std::pow(rd,-(A+2.0))+std::pow(r,-(A+2.0)));
	}
	else {
	  val = 0.0;
	  dudr = 0.0;
	  d2udr2 = 0.0;
	}
	return val;
      }

    /*
    inline real_type
      evaluate(real_type r, real_type& dudr, real_type& d2udr2, real_type& d3udr3) const
      {
        real_type u = 1.0/(1.0+B*r);
        dudr = A*u*u;
        d2udr2 = -B*B*dudr*u;
        d3udr3 = -3.0*B*d2udr2*u;
        return A*u*r;
      }
    */

    inline real_type evaluateV(const int iat, const int iStart, const int iEnd,
        const T* restrict _distArray, T* restrict distArrayCompressed ) const
    {
      real_type sum(0);
      for(int idx=iStart; idx<iEnd; idx++)
        if (idx!=iat) sum += evaluate(_distArray[idx]);
      return sum;
    }

    inline void evaluateVGL(const int iat, const int iStart, const int iEnd,
        const T* distArray,  T* restrict valArray,
        T* restrict gradArray, T* restrict laplArray,
        T* restrict distArrayCompressed, int* restrict distIndices ) const
    {
      for(int idx=iStart; idx<iEnd; idx++)
      {
        valArray[idx] = evaluate(distArray[idx], gradArray[idx], laplArray[idx]);
        gradArray[idx] /= distArray[idx];
      }
      if ( iat>=iStart && iat<iEnd )
        valArray[iat] = gradArray[iat] = laplArray[iat] = T(0);
    }

    inline real_type f(real_type r)
    {
      return evaluate(r);
    }

    inline real_type df(real_type r)
    {
      real_type dudr,d2udr2;
      real_type res=evaluate(r,dudr,d2udr2);
      return dudr;
    }

    // TODO: evaluateDerivatives() are still unchanged from the Pade version
    //       change this as necessary for opt
    /// compute derivatives with respect to A and B
    inline bool evaluateDerivatives (real_type r, std::vector<TinyVector<real_type,3> >& derivs)
    {
      int i=0;
      real_type u = 1.0/(1.0+B*r);
      if(Opt_A)
      {
        derivs[i][0]= r*u-1/B; //du/da
        derivs[i][1]= u*u; //d(du/da)/dr
        derivs[i][2]= -B*B*u*u*u; //d^2 (du/da)/dr
        ++i;
      }
      if(Opt_B)
      {
        derivs[i][0]= -A*r*r*u*u+A/(B*B); //du/db
        derivs[i][1]= -2.0*A*r*u*u*u; //d(du/db)/dr
        derivs[i][2]=  2.0*A*(B*B*r-1)*u*u*u*u; //d^2(du/db)/dr^2
      }
      return true;
    }

    /// compute derivatives with respect to A and B
    inline bool evaluateDerivatives (real_type r, std::vector<real_type>& derivs)
    {
      int i=0;
      real_type u = 1.0/(1.0+B*r);
      if(Opt_A)
      {
        derivs[i]= r*u-1/B; //du/da
        ++i;
      }
      if(Opt_B)
      {
        derivs[i]= -A*r*r*u*u+A/(B*B); //du/db
      }
      return true;
    }

    bool put(xmlNodePtr cur) {
      xmlNodePtr tcur = cur->xmlChildrenNode;
      while(tcur != NULL) {
	//@todo Var -> <param(eter) role="opt"/>
        std::string cname((const char*)(tcur->name));
        tolower(cname);
	if(cname == "parameter" || cname == "var") {
          std::string pName((const char*)(xmlGetProp(tcur,(const xmlChar *)"name")));
	  //            string idname((const char*)(xmlGetProp(tcur,(const xmlChar *)"id")));
	  if(pName == "a") {
	    putContent(A,tcur);
	  } else if(pName == "b") {
	    ID_B = (const char*)(xmlGetProp(tcur,(const xmlChar *)"id"));
	    putContent(B,tcur);
	  }
	}
	tcur = tcur->next;
      }
      reset(A,B);
      // TODO: figure out what optimize::* do
      if (Opt_A) myVars.insert(ID_A, A, Opt_A, optimize::OTHER_P);
      if (Opt_B) myVars.insert(ID_B, B, Opt_B, optimize::OTHER_P);

      real_type Y;
      c1 = evaluate(rcusp, Y, c2);
      c2 = -Y/(2.0*rcusp*c1);
      c1 *= std::exp(c2*rcusp*rcusp);
      // std::cout << "ChangMo's test: " << rcusp << ", " << cutoff_radius << ", " << c1 << ", " << c2 << std::endl << evaluate(rcusp-0.0001) << ", " << evaluate(rcusp+0.0001) << std::endl;

      return true;
    }

    void checkInVariables(opt_variables_type& active)
    {
      active.insertFrom(myVars);
      //myVars.print(std::cout);
    }

    void checkOutVariables(const opt_variables_type& active)
    {
      myVars.getIndex(active);
      //myVars.print(std::cout);
    }

    void resetParameters(const opt_variables_type& active)
    {
      int ia=myVars.where(0); if(ia>-1) A=active[ia];
      int ib=myVars.where(1); if(ib>-1) B=active[ib];
      reset(A,B);
      real_type Y;
      c1 = evaluate(rcusp, Y, c2);
      c2 = -Y/(2.0*rcusp*c1);
      c1 *= std::exp(-c2*rcusp*rcusp);
    }
  };
}
#endif

