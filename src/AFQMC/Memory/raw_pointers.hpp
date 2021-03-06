//////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source
// License.  See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by:
//    Lawrence Livermore National Laboratory 
//
// File created by:
// Miguel A. Morales, moralessilva2@llnl.gov 
//    Lawrence Livermore National Laboratory 
////////////////////////////////////////////////////////////////////////////////

#ifndef AFQMC_RAW_POINTERS_DETAIL_HPP 
#define AFQMC_RAW_POINTERS_DETAIL_HPP

#include <type_traits>
#include <complex>
#include "AFQMC/Utilities/type_conversion.hpp"

namespace qmcplusplus {
namespace afqmc {

  template<class T,
           typename = typename std::enable_if_t<std::is_fundamental<T>::value>>
  inline static T* to_address(T* p) { return p; }

  template<class T>
  inline static std::complex<T>* to_address(std::complex<T>* p) { return p; }

  template<class Q, class T,
           typename = typename std::enable_if_t<std::is_fundamental<Q>::value>,
           typename = typename std::enable_if_t<std::is_fundamental<T>::value>>
  inline static Q* pointer_cast(T* p) { return reinterpret_cast<Q*>(p); }

  template<class Q, class T,
           typename = typename std::enable_if_t<std::is_fundamental<Q>::value>>
  inline static Q* pointer_cast(std::complex<T>* p) { return reinterpret_cast<Q*>(p); }


  /************* copy_n_cast ****************/
  template<class T, class Q, class Size>
  Q* copy_n_cast(T const* A, Size n, Q* B) {
    for(Size i=0; i<n; i++, ++A, ++B)
      *B = static_cast<Q>(*A);
    return B;
  }

  /************* print ****************/
  template<typename T>
  void print(std::string str, T const* p, int n) {
    std::cout<<str <<" ";
    for(int i=0; i<n; i++)
      std::cout<<*(p+i) <<" ";
    std::cout<<std::endl;
  }

}
}

#endif
