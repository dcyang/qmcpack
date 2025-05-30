//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//                    Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//                    Raymond Clay III, j.k.rofling@gmail.com, Lawrence Livermore National Laboratory
//                    Mark A. Berrill, berrillma@ornl.gov, Oak Ridge National Laboratory
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////


#ifndef QMCPLUSPLUS_OPTIMIZE_VARIABLESET_H
#define QMCPLUSPLUS_OPTIMIZE_VARIABLESET_H
#include "config.h"
#include <map>
#include <vector>
#include <iostream>
#include <complex>
#include "Configuration.h"

namespace qmcplusplus
{
class hdf_archive;
}

namespace optimize
{
/** An enum useful for determining the type of parameter is being optimized.
*   knowing this in the opt routine can reduce the computational load.
*/
enum
{
  OTHER_P = 0,
  LOGLINEAR_P, //B-spline Jastrows
  LOGLINEAR_K, //K space Jastrows
  LINEAR_P,    //Multi-determinant coefficients
  SPO_P,       //SPO set Parameters
  BACKFLOW_P   //Backflow parameters
};

/** class to handle a set of variables that can be modified during optimizations
 *
 * A serialized container of named variables.
 */
struct VariableSet
{
  using real_type       = qmcplusplus::QMCTraits::RealType;
  using pair_type       = std::pair<std::string, real_type>;
  using index_pair_type = std::pair<std::string, int>;
  using iterator        = std::vector<pair_type>::iterator;
  using const_iterator  = std::vector<pair_type>::const_iterator;
  using size_type       = std::vector<pair_type>::size_type;

private:
  ///number of active variables
  int num_active_vars;
  std::vector<pair_type> NameAndValue;
  std::vector<index_pair_type> ParameterType;

public:
  /** store locator of the named variable
   *
   * if(Index[i]  == -1), the named variable is not active
   */
  std::vector<int> Index;
  ///default constructor
  inline VariableSet() : num_active_vars(0) {}
  ///viturval destructor for safety
  virtual ~VariableSet() = default;
  /** if any of Index value is not zero, return true
   */
  inline bool is_optimizable() const { return num_active_vars > 0; }
  ///return the number of active variables
  inline int size_of_active() const { return num_active_vars; }
  ///return the first const_iterator
  inline const_iterator begin() const { return NameAndValue.begin(); }
  ///return the last const_iterator
  inline const_iterator end() const { return NameAndValue.end(); }
  ///return the first iterator
  inline iterator begin() { return NameAndValue.begin(); }
  ///return the last iterator
  inline iterator end() { return NameAndValue.end(); }
  ///return the size
  inline size_type size() const { return NameAndValue.size(); }
  ///return the locator of the i-th Index
  inline int where(int i) const { return Index[i]; }
  /** return the iterator of a named parameter
   * @param vname name of a parameter
   * @return the locator of vname
   *
   * If vname is not found among the Names, return NameAndValue.end()
   * so that ::end() member function can be used to validate the iterator.
   */
  inline iterator find(const std::string& vname)
  {
    return std::find_if(NameAndValue.begin(), NameAndValue.end(),
                        [&vname](const auto& value) { return value.first == vname; });
  }

  /** return the Index vaule for the named parameter
   * @param vname name of the variable
   *
   * If vname is not found in this variables, return -1;
   */
  int getIndex(const std::string& vname) const;

  /* return the NameAndValue index for the named parameter
   * @ param vname name of the variable
   *
   * Differs from getIndex by not relying on the indices cached in Index
   * myVars[i] will always return the value of the parameter if it is stored
   * regardless of whether or not the Index array has been correctly reset
   *
   * if vname is not found, return -1
   *
   */
  inline int getLoc(const std::string& vname) const
  {
    int loc = 0;
    while (loc != NameAndValue.size())
    {
      if (NameAndValue[loc].first == vname)
        return loc;
      ++loc;
    }
    return -1;
  }

  inline void insert(const std::string& vname, real_type v, bool enable = true, int type = OTHER_P)
  {
    iterator loc = find(vname);
    int ind_loc  = loc - NameAndValue.begin();
    if (loc == NameAndValue.end()) //  && enable==true)
    {
      Index.push_back(ind_loc);
      NameAndValue.push_back(pair_type(vname, v));
      ParameterType.push_back(index_pair_type(vname, type));
    }
    //disable it if enable == false
    if (!enable)
      Index[ind_loc] = -1;
  }

  inline void setParameterType(int type)
  {
    std::vector<index_pair_type>::iterator PTit(ParameterType.begin()), PTend(ParameterType.end());
    while (PTit != PTend)
    {
      (*PTit).second = type;
      PTit++;
    }
  }

  inline void getParameterTypeList(std::vector<int>& types) const
  {
    auto ptit(ParameterType.begin()), ptend(ParameterType.end());
    types.resize(ptend - ptit);
    auto tit(types.begin());
    while (ptit != ptend)
      (*tit++) = (*ptit++).second;
  }


  /** equivalent to std::map<std::string,T>[string] operator
   */
  inline real_type& operator[](const std::string& vname)
  {
    iterator loc = find(vname);
    if (loc == NameAndValue.end())
    {
      Index.push_back(-1);
      NameAndValue.push_back(pair_type(vname, 0));
      ParameterType.push_back(index_pair_type(vname, 0));
      return NameAndValue.back().second;
    }
    return (*loc).second;
  }


  /** return the name of i-th variable
   * @param i index
   */
  const std::string& name(int i) const { return NameAndValue[i].first; }

  /** return the i-th value
   * @param i index
   */
  inline real_type operator[](int i) const { return NameAndValue[i].second; }

  /** assign the i-th value
   * @param i index
   */
  inline real_type& operator[](int i) { return NameAndValue[i].second; }

  /** get the i-th parameter's type
  * @param i index
  */
  inline int getType(int i) const { return ParameterType[i].second; }

  /** clear the variable set
   *
   * Remove all the data.
   */
  void clear();

  /** insert a VariableSet to the list
   * @param input variables
   */
  void insertFrom(const VariableSet& input);

  /** reset Index of active parameters
   */
  void resetIndex();

  /** set the index table of this VariableSet
   * @param selected input variables
   *
   * This VariableSet is a subset of selected.
   */
  void getIndex(const VariableSet& selected);

  /** find the index of the first parameter of *this set in the selection
   * return -1 if not found.
   */
  int findIndexOfFirstParam(const VariableSet& selected) const;

  /** set default Indices, namely all the variables are active
   */
  void setIndexDefault();

  void print(std::ostream& os, int leftPadSpaces = 0, bool printHeader = false) const;

  // Save variational parameters to an HDF file
  void writeToHDF(const std::string& filename, qmcplusplus::hdf_archive& hout) const;

  /// Read variational parameters from an HDF file.
  /// This assumes VariableSet is already set up.
  void readFromHDF(const std::string& filename, qmcplusplus::hdf_archive& hin);
};
} // namespace optimize

#endif
