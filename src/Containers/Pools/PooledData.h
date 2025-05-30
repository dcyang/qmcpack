//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2025 QMCPACK developers.
//
// File developed by: Ken Esler, kpesler@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Ye Luo, yeluo@anl.gov, Argonne National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//                    Peter W. Doak, doakpw@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


/** @file Pools/PooledData.h
 * @brief Define a serialized buffer to store anonymous data
 *
 * JK: Remove iterator version on 2016-01-04
 */
#ifndef QMCPLUSPLUS_POOLEDDATA_H
#define QMCPLUSPLUS_POOLEDDATA_H
#include <vector>
#include <complex>
#include <limits>
#include <iterator>
#include <cassert>
#include "type_traits/complex_help.hpp"

/** Dynamically pack data into a buffer of type T.

    Somehow callers must ensure that the order of add and get commands
    is the same or
    This class is defined to allow dynamically */
template<class T>
struct PooledData
{
  using value_type = T;
  using size_type  = typename std::vector<T>::size_type;

  size_type Current;
  std::vector<T> myData;

  ///default constructor
  inline PooledData() : Current(0) {}

  ///constructor with a size
  inline PooledData(size_type n) : Current(0) { myData.resize(n, T()); }

  ///return the size of the data
  inline size_type byteSize() const { return sizeof(T) * myData.size(); }

  ///return the size of the data
  inline size_type size() const { return myData.size(); }

  //inline void clear() { std::vector<T>::clear(); Current=0;}
  inline size_type current() const { return Current; }

  /** set the Current to a cursor
   * @param cur locator to which Current is assigned
   */
  inline void rewind(size_type cur = 0) { Current = cur; }

  ///return the starting iterator
  inline typename std::vector<T>::iterator begin() { return myData.begin(); }
  ///return the ending iterator
  inline typename std::vector<T>::iterator end() { return myData.end(); }

  /*@{ matching functions to std::vector functions */
  ///clear the data and set Current=0
  inline void clear()
  {
    myData.clear();
    Current = 0;
  }
  ///reserve the memory using std::vector<T>::reserve
  inline void reserve(size_type n)
  {
    myData.reserve(n);
    Current = 0;
  }
  ///resize
  inline void resize(size_type n, T val = T())
  {
    myData.resize(n, val);
    Current = 0;
  }
  ///return i-th value
  inline T operator[](size_type i) const { return myData[i]; }
  ///return i-th value to assign
  inline T& operator[](size_type i) { return myData[i]; }
  /*@}*/

  inline void add(T& x)
  {
    static_assert(!qmcplusplus::IsComplex_t<T>::value);
    Current++;
    myData.push_back(x);
  }

  inline void add(std::complex<T>& x)
  {
    Current += 2;
    myData.push_back(x.real());
    myData.push_back(x.imag());
  }

  template<class T1>
  inline void add(T1& x)
  {
    static_assert(!qmcplusplus::IsComplex_t<T1>::value);
    Current++;
    myData.push_back(static_cast<T>(x));
  }

  template<class _InputIterator>
  inline void add(_InputIterator first, _InputIterator last)
  {
    Current += last - first;
    myData.insert(myData.end(), first, last);
  }

  void add(std::complex<T>* first, std::complex<T>* last);

  void add(const std::complex<T>* first, const std::complex<T>* last);

  template<typename T1>
  void add(std::complex<T1>* first, std::complex<T1>* last);

  inline void get(T& x) { x = myData[Current++]; }

  inline void get(std::complex<T>& x)
  {
    x = std::complex<T>(myData[Current], myData[Current + 1]);
    Current += 2;
  }

  template<class T1>
  inline void get(T1& x)
  {
    x = static_cast<T1>(myData[Current++]);
  }

  template<class _OutputIterator>
  inline void get(_OutputIterator first, _OutputIterator last)
  {
    size_type now = Current;
    Current += last - first;
    copy(myData.begin() + now, myData.begin() + Current, first);
  }

  void get(std::complex<T>* first, std::complex<T>* last);

  template<typename T1>
  void get(std::complex<T1>* first, std::complex<T1>* last);

  template<typename T1>
  inline void get(T1* first, T1* last)
  {
    size_type now = Current;
    Current += last - first;
    std::copy(myData.begin() + now, myData.begin() + Current, first);
  }

  inline void put(T& x) { myData[Current++] = x; }

  inline void put(std::complex<T>& x)
  {
    myData[Current++] = x.real();
    myData[Current++] = x.imag();
  }

  template<typename T1>
  inline void put(T1& x)
  {
    myData[Current++] = static_cast<T>(x);
  }

  template<class _InputIterator>
  inline void put(_InputIterator first, _InputIterator last)
  {
    copy(first, last, myData.begin() + Current);
    Current += last - first;
  }

  template<typename T1>
  inline void put(T1* first, T1* last)
  {
    std::copy(first, last, myData.begin() + Current);
    Current += last - first;
  }

  template<typename T1>
  inline void put(std::complex<T1>* first, std::complex<T1>* last)
  {
    while (first != last)
    {
      myData[Current++] = (*first).real();
      myData[Current++] = (*first).imag();
      ++first;
    }
  }

  /** return the address of the first element **/
  inline T* data() { return myData.data(); }

  inline void print(std::ostream& os) { copy(myData.begin(), myData.end(), std::ostream_iterator<T>(os, " ")); }

  template<class Msg>
  inline Msg& putMessage(Msg& m)
  {
    m.Pack(myData.data(), myData.size());
    return m;
  }

  template<class Msg>
  inline Msg& getMessage(Msg& m)
  {
    m.Unpack(myData.data(), myData.size());
    return m;
  }

  /** Accumulate element wise values into this from other
   */
  inline PooledData<T>& operator+=(const PooledData<T>& other)
  {
    assert(myData.size() >= other.size());
    for (int i = 0; i < myData.size(); ++i)
      myData[i] += other[i];
    return *this;
  }

  inline PooledData<T>& operator*=(T scale)
  {
    for (int i = 0; i < myData.size(); ++i)
      myData[i] *= scale;
    return *this;
  }
};

/** operator to check if two buffers are identical */
template<class T>
bool operator==(const PooledData<T>& a, const PooledData<T>& b)
{
  if (a.size() != b.size())
    return false;
  //if(a.Current != b.Current) return false;
  for (typename PooledData<T>::size_type i = 0; i < a.size(); ++i)
  {
    if (std::abs(a[i] - b[i]) > std::numeric_limits<T>::epsilon())
      return false;
  }
  return true;
}

/** operator to check if two buffers are different */
template<class T>
bool operator!=(const PooledData<T>& a, const PooledData<T>& b)
{
  if (a.size() != b.size())
    return true;
  //if(a.Current != b.Current) return true;
  for (typename PooledData<T>::size_type i = 0; i < a.size(); ++i)
  {
    if (std::abs(a[i] - b[i]) > std::numeric_limits<T>::epsilon())
      return true;
  }
  return false;
}

extern template struct PooledData<float>;
extern template struct PooledData<double>;
extern template void PooledData<double>::add<float>(std::complex<float>* first, std::complex<float>* last);
extern template void PooledData<double>::get<float>(std::complex<float>* first, std::complex<float>* last);

#endif
