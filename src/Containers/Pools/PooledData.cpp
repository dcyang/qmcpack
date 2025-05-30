//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2025 QMCPACK developers
//
// File developed by: Peter W. Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Peter W. Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//////////////////////////////////////////////////////////////////////////////////////

#include "PooledData.h"
#include <cassert>

/** @file Pools/PooledData.cpp
 * @brief Explicitly instantiate some PooledData's and some templated
 * member functions.
 *
 * JK: Remove iterator version on 2016-01-04
 */

template<typename T>
void PooledData<T>::add(const std::complex<T>* first, const std::complex<T>* last)
{
  static_assert(!(qmcplusplus::IsComplex_t<T>::value));
  long long d_n = 2 * (last - first);
  assert(d_n >= 0);
  myData.insert(myData.end(), d_n, T{});
  while (first != last)
  {
    myData[Current++] = (*first).real();
    myData[Current++] = (*first).imag();
    ++first;
  }
}

template<typename T>
void PooledData<T>::add(std::complex<T>* first, std::complex<T>* last)
{
  static_assert(!(qmcplusplus::IsComplex_t<T>::value));
  long long d_n = 2 * (last - first);
  assert(d_n >= 0);
  myData.insert(myData.end(), d_n, T{});
  while (first != last)
  {
    myData[Current++] = (*first).real();
    myData[Current++] = (*first).imag();
    ++first;
  }
}

template<typename T>
template<typename T1>
void PooledData<T>::add(std::complex<T1>* first, std::complex<T1>* last)
{
  static_assert(!(qmcplusplus::IsComplex_t<T>::value));
  // only widening coversions are ok.
  static_assert(!std::is_same_v<T, T1>);
  static_assert(sizeof(T) > sizeof(T1));

  long long t_n = 2 * (last - first);
  assert(t_n > 0);
  myData.insert(myData.end(), t_n, T{});
  while (first != last)
  {
    myData[Current++] = (*first).real();
    myData[Current++] = (*first).imag();
    ++first;
  }
}

template<typename T>
void PooledData<T>::get(std::complex<T>* first, std::complex<T>* last)
{
  auto getptr = first;
  static_assert(!(qmcplusplus::IsComplex_t<T>::value));
  while (getptr != last)
  {
    getptr->real(myData[Current++]);
    getptr->imag(myData[Current++]);
    ++getptr;
  }
}

template<typename T>
template<typename T1>
void PooledData<T>::get(std::complex<T1>* first, std::complex<T1>* last)
{
  auto getptr = first;
  static_assert(!qmcplusplus::IsComplex_t<T>::value);
  static_assert(!std::is_same_v<T, T1>);
  while (getptr != last)
  {
    // Generally this is going to be a narrowing conversion.
    // However there should be symmetry between the add's and get's
    // applied to a PooledData so this should not result in a loss of
    // precision when sizeof(T) > sizeof(T1).
    getptr->real(myData[Current++]);
    getptr->imag(myData[Current++]);
    ++getptr;
  }
}

template struct PooledData<float>;
template struct PooledData<double>;
template void PooledData<double>::add<float>(std::complex<float>* first, std::complex<float>* last);
template void PooledData<double>::get<float>(std::complex<float>* first, std::complex<float>* last);
