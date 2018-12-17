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
    
    
/** @file GaussianFunctors.h
 * @brief Functors which implement 1-body correlation 
 */
#ifndef QMCPLUSPLUS_GAUSSIANFUNCTORS_H
#define QMCPLUSPLUS_GAUSSIANFUNCTORS_H
#include "Numerics/OptimizableFunctorBase.h"
#include "OhmmsData/AttributeSet.h"
#include <cmath>
// #include <vector>
#include "OhmmsPETE/TinyVector.h"


namespace qmcplusplus
{

  /** Implements a Gaussian Function \f$u[r]=A*r*r/2\f$
   *  Useful for harmonic oscillators, or near-harmonic oscillators such as solid helium
   *  Ref: PRL 2, 290 (1959)
   */
  template<class T>
    struct GaussianFunctor: public OptimizableFunctorBase
  {
    ///true, if A is optimizable
    bool Opt_A;
    ///input A
    real_type A;
    ///id of A
    std::string ID_A;

    ///default constructor
    GaussianFunctor(): ID_A("0"), A(1.0), Opt_A(true) {
      reset();
    }

    ///constructor
    explicit GaussianFunctor(real_type a):
      A(a), Opt_A(true)
    {
      reset();
    }

    /** constructor with A
     * @param a value of A
     * @param ida id of A
     */
    explicit GaussianFunctor(real_type a, const std::string& ida)
      :A(a),ID_A(ida), Opt_A(true)
    {
      reset();
    }

    OptimizableFunctorBase* makeClone() const
    {
      return new GaussianFunctor(*this);
    }

    void reset()
    {
      cutoff_radius=1.0e4; //some big range
    }

    void reset(double A_rhs)
    {
      cutoff_radius=1.0e4; //some big range
      A=A_rhs;
    }

    inline real_type evaluate(real_type r) const
    {
      return A*A*r*r;
    }

    inline real_type
      evaluate(real_type r, real_type& dudr, real_type& d2udr2) const
      {
        dudr = 2.0*A*A*r;
        d2udr2 = 2.0*A*A;
        return A*A*r*r;
      }

    inline real_type
      evaluate(real_type r, real_type& dudr, real_type& d2udr2, real_type& d3udr3) const
      {
        dudr = 2.0*A*A*r;
        d2udr2 = 2.0*A*A;
        d3udr3 = 0.0;
        return A*A*r*r;
      }

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

    /// compute derivatives with respect to A
    inline bool evaluateDerivatives (real_type r, std::vector<TinyVector<real_type,3> >& derivs)
    {
      int i=0;
      if (Opt_A)
      {
        derivs[i][0]= (A+A)*r*r;        // ∂u/∂a
        derivs[i][1]= 4.0*A*r;          // d(∂u/∂a)/dr
        derivs[i][2]= 4.0*A;            // d²(∂u/∂a)/dr²
        ++i;
      }
      return true;
    }

    /// compute derivatives with respect to A
    inline bool evaluateDerivatives (real_type r, std::vector<real_type>& derivs)
    {
      int i=0;
      if (Opt_A)
      {
        derivs[i] = (A+A)*r*r;           // ∂u/∂a
        ++i;
      }
      return true;
    }

    bool put(xmlNodePtr cur)
    {
      real_type Atemp(A);
      xmlNodePtr tcur = cur->xmlChildrenNode;
      while(tcur != NULL)
      {
        std::string cname((const char*)(tcur->name));
        tolower(cname);
        if(cname.find("param") != std::string::npos || cname.find("var") != std::string::npos)
        {
          std::string id_in("0");
          std::string p_name("A");
          OhmmsAttributeSet rAttrib;
          rAttrib.add(id_in, "id");
          rAttrib.add(p_name, "name");
          rAttrib.put(tcur);
          tolower(p_name);
          if (p_name=="a")
          {
            ID_A = id_in;
            putContent(Atemp,tcur);

            std::string optimize("yes");
            OhmmsAttributeSet pAttrib;
            pAttrib.add(optimize, "optimize");
            pAttrib.put(tcur);

            if (optimize.find("no") != std::string::npos
                || optimize.find("false") != std::string::npos
                || optimize.find("0") != std::string::npos) Opt_A = false;
          }
        }
        tcur = tcur->next;
      }
      A=Atemp;
      reset();
      // myVars.clear();        // XXX: WHY is this needed?
      app_log() << "  Gaussian Jastrow parameter A = " << A << std::endl;
      if (Opt_A) {
        myVars.insert(ID_A,A, Opt_A,optimize::LOGLINEAR_P);
        app_log() << "  A is set to be optimizable." << std::endl;
      }

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
      int i = 0;
      double A_rhs = A;
      if (Opt_A) {
        int j = myVars.where(i); if (j > -1) A_rhs = active[j];
        ++i;
      }
      reset(A_rhs);
    }
  };
}
#endif
