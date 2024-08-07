//////////////////////////////////////////////////////////////////////////////////////
// This file is distributed under the University of Illinois/NCSA Open Source License.
// See LICENSE file in top directory for details.
//
// Copyright (c) 2016 Jeongnim Kim and QMCPACK developers.
//
// File developed by: Jeremy McMinnis, jmcminis@gmail.com, University of Illinois at Urbana-Champaign
//
// File created by: Jeongnim Kim, jeongnim.kim@gmail.com, University of Illinois at Urbana-Champaign
//////////////////////////////////////////////////////////////////////////////////////

#ifndef OHMMS_MATRIXOPERATOR_H
#define OHMMS_MATRIXOPERATOR_H

///////////////////////////////////////////////////////////////////////////////
//
// WARNING: THIS FILE WAS GENERATED AUTOMATICALLY!
// YOU SHOULD MODIFY THE INPUT FILES INSTEAD OF CHANGING THIS FILE DIRECTLY!
//
// THE FOLLOWING INPUT FILES WERE USED TO MAKE THIS FILE:
//
// MakeOperators
// matrixOps.in
// MatrixDefs.in
///////////////////////////////////////////////////////////////////////////////

namespace qmcplusplus
{
template<class T1, class C1>
inline typename MakeReturn<UnaryNode<OpUnaryMinus, typename CreateLeaf<Matrix<T1, C1>>::Leaf_t>>::Expression_t
    operator-(const Matrix<T1, C1>& l)
{
  using Tree_t = UnaryNode<OpUnaryMinus, typename CreateLeaf<Matrix<T1, C1>>::Leaf_t>;
  return MakeReturn<Tree_t>::make(Tree_t(CreateLeaf<Matrix<T1, C1>>::make(l)));
}

template<class T1, class C1>
inline typename MakeReturn<UnaryNode<OpUnaryPlus, typename CreateLeaf<Matrix<T1, C1>>::Leaf_t>>::Expression_t operator+(
    const Matrix<T1, C1>& l)
{
  using Tree_t = UnaryNode<OpUnaryPlus, typename CreateLeaf<Matrix<T1, C1>>::Leaf_t>;
  return MakeReturn<Tree_t>::make(Tree_t(CreateLeaf<Matrix<T1, C1>>::make(l)));
}

template<class T1, class C1>
inline typename MakeReturn<UnaryNode<OpBitwiseNot, typename CreateLeaf<Matrix<T1, C1>>::Leaf_t>>::Expression_t
    operator~(const Matrix<T1, C1>& l)
{
  using Tree_t = UnaryNode<OpBitwiseNot, typename CreateLeaf<Matrix<T1, C1>>::Leaf_t>;
  return MakeReturn<Tree_t>::make(Tree_t(CreateLeaf<Matrix<T1, C1>>::make(l)));
}

template<class T1, class C1>
inline typename MakeReturn<UnaryNode<OpIdentity, typename CreateLeaf<Matrix<T1, C1>>::Leaf_t>>::Expression_t
    PETE_identity(const Matrix<T1, C1>& l)
{
  using Tree_t = UnaryNode<OpIdentity, typename CreateLeaf<Matrix<T1, C1>>::Leaf_t>;
  return MakeReturn<Tree_t>::make(Tree_t(CreateLeaf<Matrix<T1, C1>>::make(l)));
}

template<class T1, class T2, class C2>
inline typename MakeReturn<UnaryNode<OpCast<T1>, typename CreateLeaf<Matrix<T2, C2>>::Leaf_t>>::Expression_t peteCast(
    const T1&,
    const Matrix<T2, C2>& l)
{
  using Tree_t = UnaryNode<OpCast<T1>, typename CreateLeaf<Matrix<T2, C2>>::Leaf_t>;
  return MakeReturn<Tree_t>::make(Tree_t(CreateLeaf<Matrix<T2, C2>>::make(l)));
}

template<class T1, class C1, class T2, class C2>
inline typename MakeReturn<
    BinaryNode<OpAdd, typename CreateLeaf<Matrix<T1, C1>>::Leaf_t, typename CreateLeaf<Matrix<T2, C2>>::Leaf_t>>::
    Expression_t
    operator+(const Matrix<T1, C1>& l, const Matrix<T2, C2>& r)
{
  typedef BinaryNode<OpAdd, typename CreateLeaf<Matrix<T1, C1>>::Leaf_t, typename CreateLeaf<Matrix<T2, C2>>::Leaf_t>
      Tree_t;
  return MakeReturn<Tree_t>::make(Tree_t(CreateLeaf<Matrix<T1, C1>>::make(l), CreateLeaf<Matrix<T2, C2>>::make(r)));
}

template<class T1, class C1, class T2, class C2>
inline typename MakeReturn<
    BinaryNode<OpSubtract, typename CreateLeaf<Matrix<T1, C1>>::Leaf_t, typename CreateLeaf<Matrix<T2, C2>>::Leaf_t>>::
    Expression_t
    operator-(const Matrix<T1, C1>& l, const Matrix<T2, C2>& r)
{
  typedef BinaryNode<OpSubtract, typename CreateLeaf<Matrix<T1, C1>>::Leaf_t,
                     typename CreateLeaf<Matrix<T2, C2>>::Leaf_t>
      Tree_t;
  return MakeReturn<Tree_t>::make(Tree_t(CreateLeaf<Matrix<T1, C1>>::make(l), CreateLeaf<Matrix<T2, C2>>::make(r)));
}

/////////////////////////////////////////////////////////
template<class T1, class C1, class T2, class C2>
inline typename MakeReturn<
    BinaryNode<OpMultiply, typename CreateLeaf<Matrix<T1, C1>>::Leaf_t, typename CreateLeaf<Matrix<T2, C2>>::Leaf_t>>::
    Expression_t
    operator*(const Matrix<T1, C1>& l, const Matrix<T2, C2>& r)
{
  typedef BinaryNode<OpMultiply, typename CreateLeaf<Matrix<T1, C1>>::Leaf_t,
                     typename CreateLeaf<Matrix<T2, C2>>::Leaf_t>
      Tree_t;
  return MakeReturn<Tree_t>::make(Tree_t(CreateLeaf<Matrix<T1, C1>>::make(l), CreateLeaf<Matrix<T2, C2>>::make(r)));
}
////////////////////////////////////////////////////////
template<class T1, class C1, class T2, class C2>
inline typename MakeReturn<
    BinaryNode<OpMod, typename CreateLeaf<Matrix<T1, C1>>::Leaf_t, typename CreateLeaf<Matrix<T2, C2>>::Leaf_t>>::
    Expression_t
    operator%(const Matrix<T1, C1>& l, const Matrix<T2, C2>& r)
{
  typedef BinaryNode<OpMod, typename CreateLeaf<Matrix<T1, C1>>::Leaf_t, typename CreateLeaf<Matrix<T2, C2>>::Leaf_t>
      Tree_t;
  return MakeReturn<Tree_t>::make(Tree_t(CreateLeaf<Matrix<T1, C1>>::make(l), CreateLeaf<Matrix<T2, C2>>::make(r)));
}

template<class T1, class C1, class T2, class C2>
inline typename MakeReturn<BinaryNode<OpBitwiseAnd,
                                      typename CreateLeaf<Matrix<T1, C1>>::Leaf_t,
                                      typename CreateLeaf<Matrix<T2, C2>>::Leaf_t>>::Expression_t
    operator&(const Matrix<T1, C1>& l, const Matrix<T2, C2>& r)
{
  typedef BinaryNode<OpBitwiseAnd, typename CreateLeaf<Matrix<T1, C1>>::Leaf_t,
                     typename CreateLeaf<Matrix<T2, C2>>::Leaf_t>
      Tree_t;
  return MakeReturn<Tree_t>::make(Tree_t(CreateLeaf<Matrix<T1, C1>>::make(l), CreateLeaf<Matrix<T2, C2>>::make(r)));
}

template<class T1, class C1, class T2, class C2>
inline typename MakeReturn<
    BinaryNode<OpBitwiseOr, typename CreateLeaf<Matrix<T1, C1>>::Leaf_t, typename CreateLeaf<Matrix<T2, C2>>::Leaf_t>>::
    Expression_t
    operator|(const Matrix<T1, C1>& l, const Matrix<T2, C2>& r)
{
  typedef BinaryNode<OpBitwiseOr, typename CreateLeaf<Matrix<T1, C1>>::Leaf_t,
                     typename CreateLeaf<Matrix<T2, C2>>::Leaf_t>
      Tree_t;
  return MakeReturn<Tree_t>::make(Tree_t(CreateLeaf<Matrix<T1, C1>>::make(l), CreateLeaf<Matrix<T2, C2>>::make(r)));
}

template<class T1, class C1, class T2, class C2>
inline typename MakeReturn<BinaryNode<OpBitwiseXor,
                                      typename CreateLeaf<Matrix<T1, C1>>::Leaf_t,
                                      typename CreateLeaf<Matrix<T2, C2>>::Leaf_t>>::Expression_t
    operator^(const Matrix<T1, C1>& l, const Matrix<T2, C2>& r)
{
  typedef BinaryNode<OpBitwiseXor, typename CreateLeaf<Matrix<T1, C1>>::Leaf_t,
                     typename CreateLeaf<Matrix<T2, C2>>::Leaf_t>
      Tree_t;
  return MakeReturn<Tree_t>::make(Tree_t(CreateLeaf<Matrix<T1, C1>>::make(l), CreateLeaf<Matrix<T2, C2>>::make(r)));
}

template<class T1, class C1, class T2>
inline typename MakeReturn<
    BinaryNode<OpAdd, typename CreateLeaf<Matrix<T1, C1>>::Leaf_t, typename CreateLeaf<Expression<T2>>::Leaf_t>>::
    Expression_t
    operator+(const Matrix<T1, C1>& l, const Expression<T2>& r)
{
  typedef BinaryNode<OpAdd, typename CreateLeaf<Matrix<T1, C1>>::Leaf_t, typename CreateLeaf<Expression<T2>>::Leaf_t>
      Tree_t;
  return MakeReturn<Tree_t>::make(Tree_t(CreateLeaf<Matrix<T1, C1>>::make(l), CreateLeaf<Expression<T2>>::make(r)));
}

template<class T1, class C1, class T2>
inline typename MakeReturn<
    BinaryNode<OpSubtract, typename CreateLeaf<Matrix<T1, C1>>::Leaf_t, typename CreateLeaf<Expression<T2>>::Leaf_t>>::
    Expression_t
    operator-(const Matrix<T1, C1>& l, const Expression<T2>& r)
{
  typedef BinaryNode<OpSubtract, typename CreateLeaf<Matrix<T1, C1>>::Leaf_t,
                     typename CreateLeaf<Expression<T2>>::Leaf_t>
      Tree_t;
  return MakeReturn<Tree_t>::make(Tree_t(CreateLeaf<Matrix<T1, C1>>::make(l), CreateLeaf<Expression<T2>>::make(r)));
}

template<class T1, class C1, class T2>
inline typename MakeReturn<
    BinaryNode<OpMultiply, typename CreateLeaf<Matrix<T1, C1>>::Leaf_t, typename CreateLeaf<Expression<T2>>::Leaf_t>>::
    Expression_t
    operator*(const Matrix<T1, C1>& l, const Expression<T2>& r)
{
  typedef BinaryNode<OpMultiply, typename CreateLeaf<Matrix<T1, C1>>::Leaf_t,
                     typename CreateLeaf<Expression<T2>>::Leaf_t>
      Tree_t;
  return MakeReturn<Tree_t>::make(Tree_t(CreateLeaf<Matrix<T1, C1>>::make(l), CreateLeaf<Expression<T2>>::make(r)));
}

template<class T1, class C1, class T2>
inline typename MakeReturn<
    BinaryNode<OpMod, typename CreateLeaf<Matrix<T1, C1>>::Leaf_t, typename CreateLeaf<Expression<T2>>::Leaf_t>>::
    Expression_t
    operator%(const Matrix<T1, C1>& l, const Expression<T2>& r)
{
  typedef BinaryNode<OpMod, typename CreateLeaf<Matrix<T1, C1>>::Leaf_t, typename CreateLeaf<Expression<T2>>::Leaf_t>
      Tree_t;
  return MakeReturn<Tree_t>::make(Tree_t(CreateLeaf<Matrix<T1, C1>>::make(l), CreateLeaf<Expression<T2>>::make(r)));
}

template<class T1, class C1, class T2>
inline typename MakeReturn<BinaryNode<OpBitwiseAnd,
                                      typename CreateLeaf<Matrix<T1, C1>>::Leaf_t,
                                      typename CreateLeaf<Expression<T2>>::Leaf_t>>::Expression_t
    operator&(const Matrix<T1, C1>& l, const Expression<T2>& r)
{
  typedef BinaryNode<OpBitwiseAnd, typename CreateLeaf<Matrix<T1, C1>>::Leaf_t,
                     typename CreateLeaf<Expression<T2>>::Leaf_t>
      Tree_t;
  return MakeReturn<Tree_t>::make(Tree_t(CreateLeaf<Matrix<T1, C1>>::make(l), CreateLeaf<Expression<T2>>::make(r)));
}

template<class T1, class C1, class T2>
inline typename MakeReturn<
    BinaryNode<OpBitwiseOr, typename CreateLeaf<Matrix<T1, C1>>::Leaf_t, typename CreateLeaf<Expression<T2>>::Leaf_t>>::
    Expression_t
    operator|(const Matrix<T1, C1>& l, const Expression<T2>& r)
{
  typedef BinaryNode<OpBitwiseOr, typename CreateLeaf<Matrix<T1, C1>>::Leaf_t,
                     typename CreateLeaf<Expression<T2>>::Leaf_t>
      Tree_t;
  return MakeReturn<Tree_t>::make(Tree_t(CreateLeaf<Matrix<T1, C1>>::make(l), CreateLeaf<Expression<T2>>::make(r)));
}

template<class T1, class C1, class T2>
inline typename MakeReturn<BinaryNode<OpBitwiseXor,
                                      typename CreateLeaf<Matrix<T1, C1>>::Leaf_t,
                                      typename CreateLeaf<Expression<T2>>::Leaf_t>>::Expression_t
    operator^(const Matrix<T1, C1>& l, const Expression<T2>& r)
{
  typedef BinaryNode<OpBitwiseXor, typename CreateLeaf<Matrix<T1, C1>>::Leaf_t,
                     typename CreateLeaf<Expression<T2>>::Leaf_t>
      Tree_t;
  return MakeReturn<Tree_t>::make(Tree_t(CreateLeaf<Matrix<T1, C1>>::make(l), CreateLeaf<Expression<T2>>::make(r)));
}

template<class T1, class T2, class C2>
inline typename MakeReturn<
    BinaryNode<OpAdd, typename CreateLeaf<Expression<T1>>::Leaf_t, typename CreateLeaf<Matrix<T2, C2>>::Leaf_t>>::
    Expression_t
    operator+(const Expression<T1>& l, const Matrix<T2, C2>& r)
{
  typedef BinaryNode<OpAdd, typename CreateLeaf<Expression<T1>>::Leaf_t, typename CreateLeaf<Matrix<T2, C2>>::Leaf_t>
      Tree_t;
  return MakeReturn<Tree_t>::make(Tree_t(CreateLeaf<Expression<T1>>::make(l), CreateLeaf<Matrix<T2, C2>>::make(r)));
}

template<class T1, class T2, class C2>
inline typename MakeReturn<
    BinaryNode<OpSubtract, typename CreateLeaf<Expression<T1>>::Leaf_t, typename CreateLeaf<Matrix<T2, C2>>::Leaf_t>>::
    Expression_t
    operator-(const Expression<T1>& l, const Matrix<T2, C2>& r)
{
  typedef BinaryNode<OpSubtract, typename CreateLeaf<Expression<T1>>::Leaf_t,
                     typename CreateLeaf<Matrix<T2, C2>>::Leaf_t>
      Tree_t;
  return MakeReturn<Tree_t>::make(Tree_t(CreateLeaf<Expression<T1>>::make(l), CreateLeaf<Matrix<T2, C2>>::make(r)));
}

template<class T1, class T2, class C2>
inline typename MakeReturn<
    BinaryNode<OpMultiply, typename CreateLeaf<Expression<T1>>::Leaf_t, typename CreateLeaf<Matrix<T2, C2>>::Leaf_t>>::
    Expression_t
    operator*(const Expression<T1>& l, const Matrix<T2, C2>& r)
{
  typedef BinaryNode<OpMultiply, typename CreateLeaf<Expression<T1>>::Leaf_t,
                     typename CreateLeaf<Matrix<T2, C2>>::Leaf_t>
      Tree_t;
  return MakeReturn<Tree_t>::make(Tree_t(CreateLeaf<Expression<T1>>::make(l), CreateLeaf<Matrix<T2, C2>>::make(r)));
}

template<class T1, class T2, class C2>
inline typename MakeReturn<
    BinaryNode<OpMod, typename CreateLeaf<Expression<T1>>::Leaf_t, typename CreateLeaf<Matrix<T2, C2>>::Leaf_t>>::
    Expression_t
    operator%(const Expression<T1>& l, const Matrix<T2, C2>& r)
{
  typedef BinaryNode<OpMod, typename CreateLeaf<Expression<T1>>::Leaf_t, typename CreateLeaf<Matrix<T2, C2>>::Leaf_t>
      Tree_t;
  return MakeReturn<Tree_t>::make(Tree_t(CreateLeaf<Expression<T1>>::make(l), CreateLeaf<Matrix<T2, C2>>::make(r)));
}

template<class T1, class T2, class C2>
inline typename MakeReturn<BinaryNode<OpBitwiseAnd,
                                      typename CreateLeaf<Expression<T1>>::Leaf_t,
                                      typename CreateLeaf<Matrix<T2, C2>>::Leaf_t>>::Expression_t
    operator&(const Expression<T1>& l, const Matrix<T2, C2>& r)
{
  typedef BinaryNode<OpBitwiseAnd, typename CreateLeaf<Expression<T1>>::Leaf_t,
                     typename CreateLeaf<Matrix<T2, C2>>::Leaf_t>
      Tree_t;
  return MakeReturn<Tree_t>::make(Tree_t(CreateLeaf<Expression<T1>>::make(l), CreateLeaf<Matrix<T2, C2>>::make(r)));
}

template<class T1, class T2, class C2>
inline typename MakeReturn<
    BinaryNode<OpBitwiseOr, typename CreateLeaf<Expression<T1>>::Leaf_t, typename CreateLeaf<Matrix<T2, C2>>::Leaf_t>>::
    Expression_t
    operator|(const Expression<T1>& l, const Matrix<T2, C2>& r)
{
  typedef BinaryNode<OpBitwiseOr, typename CreateLeaf<Expression<T1>>::Leaf_t,
                     typename CreateLeaf<Matrix<T2, C2>>::Leaf_t>
      Tree_t;
  return MakeReturn<Tree_t>::make(Tree_t(CreateLeaf<Expression<T1>>::make(l), CreateLeaf<Matrix<T2, C2>>::make(r)));
}

template<class T1, class T2, class C2>
inline typename MakeReturn<BinaryNode<OpBitwiseXor,
                                      typename CreateLeaf<Expression<T1>>::Leaf_t,
                                      typename CreateLeaf<Matrix<T2, C2>>::Leaf_t>>::Expression_t
    operator^(const Expression<T1>& l, const Matrix<T2, C2>& r)
{
  typedef BinaryNode<OpBitwiseXor, typename CreateLeaf<Expression<T1>>::Leaf_t,
                     typename CreateLeaf<Matrix<T2, C2>>::Leaf_t>
      Tree_t;
  return MakeReturn<Tree_t>::make(Tree_t(CreateLeaf<Expression<T1>>::make(l), CreateLeaf<Matrix<T2, C2>>::make(r)));
}

template<class T1, class C1, class T2>
inline typename MakeReturn<
    BinaryNode<OpAdd, typename CreateLeaf<Matrix<T1, C1>>::Leaf_t, typename CreateLeaf<T2>::Leaf_t>>::Expression_t
    operator+(const Matrix<T1, C1>& l, const T2& r)
{
  using Tree_t = BinaryNode<OpAdd, typename CreateLeaf<Matrix<T1, C1>>::Leaf_t, typename CreateLeaf<T2>::Leaf_t>;
  return MakeReturn<Tree_t>::make(Tree_t(CreateLeaf<Matrix<T1, C1>>::make(l), CreateLeaf<T2>::make(r)));
}

template<class T1, class C1, class T2>
inline typename MakeReturn<
    BinaryNode<OpSubtract, typename CreateLeaf<Matrix<T1, C1>>::Leaf_t, typename CreateLeaf<T2>::Leaf_t>>::Expression_t
    operator-(const Matrix<T1, C1>& l, const T2& r)
{
  using Tree_t = BinaryNode<OpSubtract, typename CreateLeaf<Matrix<T1, C1>>::Leaf_t, typename CreateLeaf<T2>::Leaf_t>;
  return MakeReturn<Tree_t>::make(Tree_t(CreateLeaf<Matrix<T1, C1>>::make(l), CreateLeaf<T2>::make(r)));
}

template<class T1, class C1, class T2>
inline typename MakeReturn<
    BinaryNode<OpMultiply, typename CreateLeaf<Matrix<T1, C1>>::Leaf_t, typename CreateLeaf<T2>::Leaf_t>>::Expression_t
    operator*(const Matrix<T1, C1>& l, const T2& r)
{
  using Tree_t = BinaryNode<OpMultiply, typename CreateLeaf<Matrix<T1, C1>>::Leaf_t, typename CreateLeaf<T2>::Leaf_t>;
  return MakeReturn<Tree_t>::make(Tree_t(CreateLeaf<Matrix<T1, C1>>::make(l), CreateLeaf<T2>::make(r)));
}

template<class T1, class C1, class T2>
inline typename MakeReturn<
    BinaryNode<OpMod, typename CreateLeaf<Matrix<T1, C1>>::Leaf_t, typename CreateLeaf<T2>::Leaf_t>>::Expression_t
    operator%(const Matrix<T1, C1>& l, const T2& r)
{
  using Tree_t = BinaryNode<OpMod, typename CreateLeaf<Matrix<T1, C1>>::Leaf_t, typename CreateLeaf<T2>::Leaf_t>;
  return MakeReturn<Tree_t>::make(Tree_t(CreateLeaf<Matrix<T1, C1>>::make(l), CreateLeaf<T2>::make(r)));
}

template<class T1, class C1, class T2>
inline typename MakeReturn<
    BinaryNode<OpBitwiseAnd, typename CreateLeaf<Matrix<T1, C1>>::Leaf_t, typename CreateLeaf<T2>::Leaf_t>>::
    Expression_t
    operator&(const Matrix<T1, C1>& l, const T2& r)
{
  using Tree_t = BinaryNode<OpBitwiseAnd, typename CreateLeaf<Matrix<T1, C1>>::Leaf_t, typename CreateLeaf<T2>::Leaf_t>;
  return MakeReturn<Tree_t>::make(Tree_t(CreateLeaf<Matrix<T1, C1>>::make(l), CreateLeaf<T2>::make(r)));
}

template<class T1, class C1, class T2>
inline typename MakeReturn<
    BinaryNode<OpBitwiseOr, typename CreateLeaf<Matrix<T1, C1>>::Leaf_t, typename CreateLeaf<T2>::Leaf_t>>::Expression_t
    operator|(const Matrix<T1, C1>& l, const T2& r)
{
  using Tree_t = BinaryNode<OpBitwiseOr, typename CreateLeaf<Matrix<T1, C1>>::Leaf_t, typename CreateLeaf<T2>::Leaf_t>;
  return MakeReturn<Tree_t>::make(Tree_t(CreateLeaf<Matrix<T1, C1>>::make(l), CreateLeaf<T2>::make(r)));
}

template<class T1, class C1, class T2>
inline typename MakeReturn<
    BinaryNode<OpBitwiseXor, typename CreateLeaf<Matrix<T1, C1>>::Leaf_t, typename CreateLeaf<T2>::Leaf_t>>::
    Expression_t
    operator^(const Matrix<T1, C1>& l, const T2& r)
{
  using Tree_t = BinaryNode<OpBitwiseXor, typename CreateLeaf<Matrix<T1, C1>>::Leaf_t, typename CreateLeaf<T2>::Leaf_t>;
  return MakeReturn<Tree_t>::make(Tree_t(CreateLeaf<Matrix<T1, C1>>::make(l), CreateLeaf<T2>::make(r)));
}

template<class T1, class T2, class C2>
inline typename MakeReturn<
    BinaryNode<OpAdd, typename CreateLeaf<T1>::Leaf_t, typename CreateLeaf<Matrix<T2, C2>>::Leaf_t>>::Expression_t
    operator+(const T1& l, const Matrix<T2, C2>& r)
{
  using Tree_t = BinaryNode<OpAdd, typename CreateLeaf<T1>::Leaf_t, typename CreateLeaf<Matrix<T2, C2>>::Leaf_t>;
  return MakeReturn<Tree_t>::make(Tree_t(CreateLeaf<T1>::make(l), CreateLeaf<Matrix<T2, C2>>::make(r)));
}

template<class T1, class T2, class C2>
inline typename MakeReturn<
    BinaryNode<OpSubtract, typename CreateLeaf<T1>::Leaf_t, typename CreateLeaf<Matrix<T2, C2>>::Leaf_t>>::Expression_t
    operator-(const T1& l, const Matrix<T2, C2>& r)
{
  using Tree_t = BinaryNode<OpSubtract, typename CreateLeaf<T1>::Leaf_t, typename CreateLeaf<Matrix<T2, C2>>::Leaf_t>;
  return MakeReturn<Tree_t>::make(Tree_t(CreateLeaf<T1>::make(l), CreateLeaf<Matrix<T2, C2>>::make(r)));
}

template<class T1, class T2, class C2>
inline typename MakeReturn<
    BinaryNode<OpMultiply, typename CreateLeaf<T1>::Leaf_t, typename CreateLeaf<Matrix<T2, C2>>::Leaf_t>>::Expression_t
    operator*(const T1& l, const Matrix<T2, C2>& r)
{
  using Tree_t = BinaryNode<OpMultiply, typename CreateLeaf<T1>::Leaf_t, typename CreateLeaf<Matrix<T2, C2>>::Leaf_t>;
  return MakeReturn<Tree_t>::make(Tree_t(CreateLeaf<T1>::make(l), CreateLeaf<Matrix<T2, C2>>::make(r)));
}

template<class T1, class T2, class C2>
inline typename MakeReturn<
    BinaryNode<OpMod, typename CreateLeaf<T1>::Leaf_t, typename CreateLeaf<Matrix<T2, C2>>::Leaf_t>>::Expression_t
    operator%(const T1& l, const Matrix<T2, C2>& r)
{
  using Tree_t = BinaryNode<OpMod, typename CreateLeaf<T1>::Leaf_t, typename CreateLeaf<Matrix<T2, C2>>::Leaf_t>;
  return MakeReturn<Tree_t>::make(Tree_t(CreateLeaf<T1>::make(l), CreateLeaf<Matrix<T2, C2>>::make(r)));
}

template<class T1, class T2, class C2>
inline typename MakeReturn<
    BinaryNode<OpBitwiseAnd, typename CreateLeaf<T1>::Leaf_t, typename CreateLeaf<Matrix<T2, C2>>::Leaf_t>>::
    Expression_t
    operator&(const T1& l, const Matrix<T2, C2>& r)
{
  using Tree_t = BinaryNode<OpBitwiseAnd, typename CreateLeaf<T1>::Leaf_t, typename CreateLeaf<Matrix<T2, C2>>::Leaf_t>;
  return MakeReturn<Tree_t>::make(Tree_t(CreateLeaf<T1>::make(l), CreateLeaf<Matrix<T2, C2>>::make(r)));
}

template<class T1, class T2, class C2>
inline typename MakeReturn<
    BinaryNode<OpBitwiseOr, typename CreateLeaf<T1>::Leaf_t, typename CreateLeaf<Matrix<T2, C2>>::Leaf_t>>::Expression_t
    operator|(const T1& l, const Matrix<T2, C2>& r)
{
  using Tree_t = BinaryNode<OpBitwiseOr, typename CreateLeaf<T1>::Leaf_t, typename CreateLeaf<Matrix<T2, C2>>::Leaf_t>;
  return MakeReturn<Tree_t>::make(Tree_t(CreateLeaf<T1>::make(l), CreateLeaf<Matrix<T2, C2>>::make(r)));
}

template<class T1, class T2, class C2>
inline typename MakeReturn<
    BinaryNode<OpBitwiseXor, typename CreateLeaf<T1>::Leaf_t, typename CreateLeaf<Matrix<T2, C2>>::Leaf_t>>::
    Expression_t
    operator^(const T1& l, const Matrix<T2, C2>& r)
{
  using Tree_t = BinaryNode<OpBitwiseXor, typename CreateLeaf<T1>::Leaf_t, typename CreateLeaf<Matrix<T2, C2>>::Leaf_t>;
  return MakeReturn<Tree_t>::make(Tree_t(CreateLeaf<T1>::make(l), CreateLeaf<Matrix<T2, C2>>::make(r)));
}
#ifdef PETE_ALLOW_SCALAR_SHIFT
#endif // PETE_ALLOW_SCALAR_SHIFT

template<class T1, class C1, class T2, class T3>
inline typename MakeReturn<TrinaryNode<FnWhere,
                                       typename CreateLeaf<Matrix<T1, C1>>::Leaf_t,
                                       typename CreateLeaf<T2>::Leaf_t,
                                       typename CreateLeaf<T3>::Leaf_t>>::Expression_t
    where(const Matrix<T1, C1>& c, const T2& t, const T3& f)
{
  typedef TrinaryNode<FnWhere, typename CreateLeaf<Matrix<T1, C1>>::Leaf_t, typename CreateLeaf<T2>::Leaf_t,
                      typename CreateLeaf<T3>::Leaf_t>
      Tree_t;
  return MakeReturn<Tree_t>::make(
      Tree_t(CreateLeaf<Matrix<T1, C1>>::make(c), CreateLeaf<T2>::make(t), CreateLeaf<T3>::make(f)));
}
#ifndef PETE_EXPRESSION_OPERATORS
#define PETE_EXPRESSION_OPERATORS

template<class T1>
inline typename MakeReturn<UnaryNode<OpUnaryMinus, typename CreateLeaf<Expression<T1>>::Leaf_t>>::Expression_t
    operator-(const Expression<T1>& l)
{
  using Tree_t = UnaryNode<OpUnaryMinus, typename CreateLeaf<Expression<T1>>::Leaf_t>;
  return MakeReturn<Tree_t>::make(Tree_t(CreateLeaf<Expression<T1>>::make(l)));
}

template<class T1>
inline typename MakeReturn<UnaryNode<OpUnaryPlus, typename CreateLeaf<Expression<T1>>::Leaf_t>>::Expression_t operator+(
    const Expression<T1>& l)
{
  using Tree_t = UnaryNode<OpUnaryPlus, typename CreateLeaf<Expression<T1>>::Leaf_t>;
  return MakeReturn<Tree_t>::make(Tree_t(CreateLeaf<Expression<T1>>::make(l)));
}

template<class T1>
inline typename MakeReturn<UnaryNode<OpBitwiseNot, typename CreateLeaf<Expression<T1>>::Leaf_t>>::Expression_t
    operator~(const Expression<T1>& l)
{
  using Tree_t = UnaryNode<OpBitwiseNot, typename CreateLeaf<Expression<T1>>::Leaf_t>;
  return MakeReturn<Tree_t>::make(Tree_t(CreateLeaf<Expression<T1>>::make(l)));
}

template<class T1>
inline typename MakeReturn<UnaryNode<OpIdentity, typename CreateLeaf<Expression<T1>>::Leaf_t>>::Expression_t
    PETE_identity(const Expression<T1>& l)
{
  using Tree_t = UnaryNode<OpIdentity, typename CreateLeaf<Expression<T1>>::Leaf_t>;
  return MakeReturn<Tree_t>::make(Tree_t(CreateLeaf<Expression<T1>>::make(l)));
}

template<class T1, class T2>
inline typename MakeReturn<UnaryNode<OpCast<T1>, typename CreateLeaf<Expression<T2>>::Leaf_t>>::Expression_t peteCast(
    const T1&,
    const Expression<T2>& l)
{
  using Tree_t = UnaryNode<OpCast<T1>, typename CreateLeaf<Expression<T2>>::Leaf_t>;
  return MakeReturn<Tree_t>::make(Tree_t(CreateLeaf<Expression<T2>>::make(l)));
}

template<class T1, class T2>
inline typename MakeReturn<
    BinaryNode<OpAdd, typename CreateLeaf<Expression<T1>>::Leaf_t, typename CreateLeaf<Expression<T2>>::Leaf_t>>::
    Expression_t
    operator+(const Expression<T1>& l, const Expression<T2>& r)
{
  typedef BinaryNode<OpAdd, typename CreateLeaf<Expression<T1>>::Leaf_t, typename CreateLeaf<Expression<T2>>::Leaf_t>
      Tree_t;
  return MakeReturn<Tree_t>::make(Tree_t(CreateLeaf<Expression<T1>>::make(l), CreateLeaf<Expression<T2>>::make(r)));
}

template<class T1, class T2>
inline typename MakeReturn<
    BinaryNode<OpSubtract, typename CreateLeaf<Expression<T1>>::Leaf_t, typename CreateLeaf<Expression<T2>>::Leaf_t>>::
    Expression_t
    operator-(const Expression<T1>& l, const Expression<T2>& r)
{
  typedef BinaryNode<OpSubtract, typename CreateLeaf<Expression<T1>>::Leaf_t,
                     typename CreateLeaf<Expression<T2>>::Leaf_t>
      Tree_t;
  return MakeReturn<Tree_t>::make(Tree_t(CreateLeaf<Expression<T1>>::make(l), CreateLeaf<Expression<T2>>::make(r)));
}

template<class T1, class T2>
inline typename MakeReturn<
    BinaryNode<OpMultiply, typename CreateLeaf<Expression<T1>>::Leaf_t, typename CreateLeaf<Expression<T2>>::Leaf_t>>::
    Expression_t
    operator*(const Expression<T1>& l, const Expression<T2>& r)
{
  typedef BinaryNode<OpMultiply, typename CreateLeaf<Expression<T1>>::Leaf_t,
                     typename CreateLeaf<Expression<T2>>::Leaf_t>
      Tree_t;
  return MakeReturn<Tree_t>::make(Tree_t(CreateLeaf<Expression<T1>>::make(l), CreateLeaf<Expression<T2>>::make(r)));
}

template<class T1, class T2>
inline typename MakeReturn<
    BinaryNode<OpMod, typename CreateLeaf<Expression<T1>>::Leaf_t, typename CreateLeaf<Expression<T2>>::Leaf_t>>::
    Expression_t
    operator%(const Expression<T1>& l, const Expression<T2>& r)
{
  typedef BinaryNode<OpMod, typename CreateLeaf<Expression<T1>>::Leaf_t, typename CreateLeaf<Expression<T2>>::Leaf_t>
      Tree_t;
  return MakeReturn<Tree_t>::make(Tree_t(CreateLeaf<Expression<T1>>::make(l), CreateLeaf<Expression<T2>>::make(r)));
}

template<class T1, class T2>
inline typename MakeReturn<BinaryNode<OpBitwiseAnd,
                                      typename CreateLeaf<Expression<T1>>::Leaf_t,
                                      typename CreateLeaf<Expression<T2>>::Leaf_t>>::Expression_t
    operator&(const Expression<T1>& l, const Expression<T2>& r)
{
  typedef BinaryNode<OpBitwiseAnd, typename CreateLeaf<Expression<T1>>::Leaf_t,
                     typename CreateLeaf<Expression<T2>>::Leaf_t>
      Tree_t;
  return MakeReturn<Tree_t>::make(Tree_t(CreateLeaf<Expression<T1>>::make(l), CreateLeaf<Expression<T2>>::make(r)));
}

template<class T1, class T2>
inline typename MakeReturn<
    BinaryNode<OpBitwiseOr, typename CreateLeaf<Expression<T1>>::Leaf_t, typename CreateLeaf<Expression<T2>>::Leaf_t>>::
    Expression_t
    operator|(const Expression<T1>& l, const Expression<T2>& r)
{
  typedef BinaryNode<OpBitwiseOr, typename CreateLeaf<Expression<T1>>::Leaf_t,
                     typename CreateLeaf<Expression<T2>>::Leaf_t>
      Tree_t;
  return MakeReturn<Tree_t>::make(Tree_t(CreateLeaf<Expression<T1>>::make(l), CreateLeaf<Expression<T2>>::make(r)));
}

template<class T1, class T2>
inline typename MakeReturn<BinaryNode<OpBitwiseXor,
                                      typename CreateLeaf<Expression<T1>>::Leaf_t,
                                      typename CreateLeaf<Expression<T2>>::Leaf_t>>::Expression_t
    operator^(const Expression<T1>& l, const Expression<T2>& r)
{
  typedef BinaryNode<OpBitwiseXor, typename CreateLeaf<Expression<T1>>::Leaf_t,
                     typename CreateLeaf<Expression<T2>>::Leaf_t>
      Tree_t;
  return MakeReturn<Tree_t>::make(Tree_t(CreateLeaf<Expression<T1>>::make(l), CreateLeaf<Expression<T2>>::make(r)));
}

template<class T1, class T2>
inline typename MakeReturn<
    BinaryNode<OpAdd, typename CreateLeaf<Expression<T1>>::Leaf_t, typename CreateLeaf<T2>::Leaf_t>>::Expression_t
    operator+(const Expression<T1>& l, const T2& r)
{
  using Tree_t = BinaryNode<OpAdd, typename CreateLeaf<Expression<T1>>::Leaf_t, typename CreateLeaf<T2>::Leaf_t>;
  return MakeReturn<Tree_t>::make(Tree_t(CreateLeaf<Expression<T1>>::make(l), CreateLeaf<T2>::make(r)));
}

template<class T1, class T2>
inline typename MakeReturn<
    BinaryNode<OpSubtract, typename CreateLeaf<Expression<T1>>::Leaf_t, typename CreateLeaf<T2>::Leaf_t>>::Expression_t
    operator-(const Expression<T1>& l, const T2& r)
{
  using Tree_t = BinaryNode<OpSubtract, typename CreateLeaf<Expression<T1>>::Leaf_t, typename CreateLeaf<T2>::Leaf_t>;
  return MakeReturn<Tree_t>::make(Tree_t(CreateLeaf<Expression<T1>>::make(l), CreateLeaf<T2>::make(r)));
}

template<class T1, class T2>
inline typename MakeReturn<
    BinaryNode<OpMultiply, typename CreateLeaf<Expression<T1>>::Leaf_t, typename CreateLeaf<T2>::Leaf_t>>::Expression_t
    operator*(const Expression<T1>& l, const T2& r)
{
  using Tree_t = BinaryNode<OpMultiply, typename CreateLeaf<Expression<T1>>::Leaf_t, typename CreateLeaf<T2>::Leaf_t>;
  return MakeReturn<Tree_t>::make(Tree_t(CreateLeaf<Expression<T1>>::make(l), CreateLeaf<T2>::make(r)));
}

template<class T1, class T2>
inline typename MakeReturn<
    BinaryNode<OpMod, typename CreateLeaf<Expression<T1>>::Leaf_t, typename CreateLeaf<T2>::Leaf_t>>::Expression_t
    operator%(const Expression<T1>& l, const T2& r)
{
  using Tree_t = BinaryNode<OpMod, typename CreateLeaf<Expression<T1>>::Leaf_t, typename CreateLeaf<T2>::Leaf_t>;
  return MakeReturn<Tree_t>::make(Tree_t(CreateLeaf<Expression<T1>>::make(l), CreateLeaf<T2>::make(r)));
}

template<class T1, class T2>
inline typename MakeReturn<
    BinaryNode<OpBitwiseAnd, typename CreateLeaf<Expression<T1>>::Leaf_t, typename CreateLeaf<T2>::Leaf_t>>::
    Expression_t
    operator&(const Expression<T1>& l, const T2& r)
{
  using Tree_t = BinaryNode<OpBitwiseAnd, typename CreateLeaf<Expression<T1>>::Leaf_t, typename CreateLeaf<T2>::Leaf_t>;
  return MakeReturn<Tree_t>::make(Tree_t(CreateLeaf<Expression<T1>>::make(l), CreateLeaf<T2>::make(r)));
}

template<class T1, class T2>
inline typename MakeReturn<
    BinaryNode<OpBitwiseOr, typename CreateLeaf<Expression<T1>>::Leaf_t, typename CreateLeaf<T2>::Leaf_t>>::Expression_t
    operator|(const Expression<T1>& l, const T2& r)
{
  using Tree_t = BinaryNode<OpBitwiseOr, typename CreateLeaf<Expression<T1>>::Leaf_t, typename CreateLeaf<T2>::Leaf_t>;
  return MakeReturn<Tree_t>::make(Tree_t(CreateLeaf<Expression<T1>>::make(l), CreateLeaf<T2>::make(r)));
}

template<class T1, class T2>
inline typename MakeReturn<
    BinaryNode<OpBitwiseXor, typename CreateLeaf<Expression<T1>>::Leaf_t, typename CreateLeaf<T2>::Leaf_t>>::
    Expression_t
    operator^(const Expression<T1>& l, const T2& r)
{
  using Tree_t = BinaryNode<OpBitwiseXor, typename CreateLeaf<Expression<T1>>::Leaf_t, typename CreateLeaf<T2>::Leaf_t>;
  return MakeReturn<Tree_t>::make(Tree_t(CreateLeaf<Expression<T1>>::make(l), CreateLeaf<T2>::make(r)));
}

template<class T1, class T2>
inline typename MakeReturn<
    BinaryNode<OpAdd, typename CreateLeaf<T1>::Leaf_t, typename CreateLeaf<Expression<T2>>::Leaf_t>>::Expression_t
    operator+(const T1& l, const Expression<T2>& r)
{
  using Tree_t = BinaryNode<OpAdd, typename CreateLeaf<T1>::Leaf_t, typename CreateLeaf<Expression<T2>>::Leaf_t>;
  return MakeReturn<Tree_t>::make(Tree_t(CreateLeaf<T1>::make(l), CreateLeaf<Expression<T2>>::make(r)));
}

template<class T1, class T2>
inline typename MakeReturn<
    BinaryNode<OpSubtract, typename CreateLeaf<T1>::Leaf_t, typename CreateLeaf<Expression<T2>>::Leaf_t>>::Expression_t
    operator-(const T1& l, const Expression<T2>& r)
{
  using Tree_t = BinaryNode<OpSubtract, typename CreateLeaf<T1>::Leaf_t, typename CreateLeaf<Expression<T2>>::Leaf_t>;
  return MakeReturn<Tree_t>::make(Tree_t(CreateLeaf<T1>::make(l), CreateLeaf<Expression<T2>>::make(r)));
}

template<class T1, class T2>
inline typename MakeReturn<
    BinaryNode<OpMultiply, typename CreateLeaf<T1>::Leaf_t, typename CreateLeaf<Expression<T2>>::Leaf_t>>::Expression_t
    operator*(const T1& l, const Expression<T2>& r)
{
  using Tree_t = BinaryNode<OpMultiply, typename CreateLeaf<T1>::Leaf_t, typename CreateLeaf<Expression<T2>>::Leaf_t>;
  return MakeReturn<Tree_t>::make(Tree_t(CreateLeaf<T1>::make(l), CreateLeaf<Expression<T2>>::make(r)));
}

template<class T1, class T2>
inline typename MakeReturn<
    BinaryNode<OpMod, typename CreateLeaf<T1>::Leaf_t, typename CreateLeaf<Expression<T2>>::Leaf_t>>::Expression_t
    operator%(const T1& l, const Expression<T2>& r)
{
  using Tree_t = BinaryNode<OpMod, typename CreateLeaf<T1>::Leaf_t, typename CreateLeaf<Expression<T2>>::Leaf_t>;
  return MakeReturn<Tree_t>::make(Tree_t(CreateLeaf<T1>::make(l), CreateLeaf<Expression<T2>>::make(r)));
}

template<class T1, class T2>
inline typename MakeReturn<
    BinaryNode<OpBitwiseAnd, typename CreateLeaf<T1>::Leaf_t, typename CreateLeaf<Expression<T2>>::Leaf_t>>::
    Expression_t
    operator&(const T1& l, const Expression<T2>& r)
{
  using Tree_t = BinaryNode<OpBitwiseAnd, typename CreateLeaf<T1>::Leaf_t, typename CreateLeaf<Expression<T2>>::Leaf_t>;
  return MakeReturn<Tree_t>::make(Tree_t(CreateLeaf<T1>::make(l), CreateLeaf<Expression<T2>>::make(r)));
}

template<class T1, class T2>
inline typename MakeReturn<
    BinaryNode<OpBitwiseOr, typename CreateLeaf<T1>::Leaf_t, typename CreateLeaf<Expression<T2>>::Leaf_t>>::Expression_t
    operator|(const T1& l, const Expression<T2>& r)
{
  using Tree_t = BinaryNode<OpBitwiseOr, typename CreateLeaf<T1>::Leaf_t, typename CreateLeaf<Expression<T2>>::Leaf_t>;
  return MakeReturn<Tree_t>::make(Tree_t(CreateLeaf<T1>::make(l), CreateLeaf<Expression<T2>>::make(r)));
}

template<class T1, class T2>
inline typename MakeReturn<
    BinaryNode<OpBitwiseXor, typename CreateLeaf<T1>::Leaf_t, typename CreateLeaf<Expression<T2>>::Leaf_t>>::
    Expression_t
    operator^(const T1& l, const Expression<T2>& r)
{
  using Tree_t = BinaryNode<OpBitwiseXor, typename CreateLeaf<T1>::Leaf_t, typename CreateLeaf<Expression<T2>>::Leaf_t>;
  return MakeReturn<Tree_t>::make(Tree_t(CreateLeaf<T1>::make(l), CreateLeaf<Expression<T2>>::make(r)));
}
#ifdef PETE_ALLOW_SCALAR_SHIFT
#endif // PETE_ALLOW_SCALAR_SHIFT

template<class T1, class T2, class T3>
inline typename MakeReturn<TrinaryNode<FnWhere,
                                       typename CreateLeaf<Expression<T1>>::Leaf_t,
                                       typename CreateLeaf<T2>::Leaf_t,
                                       typename CreateLeaf<T3>::Leaf_t>>::Expression_t
    where(const Expression<T1>& c, const T2& t, const T3& f)
{
  typedef TrinaryNode<FnWhere, typename CreateLeaf<Expression<T1>>::Leaf_t, typename CreateLeaf<T2>::Leaf_t,
                      typename CreateLeaf<T3>::Leaf_t>
      Tree_t;
  return MakeReturn<Tree_t>::make(
      Tree_t(CreateLeaf<Expression<T1>>::make(c), CreateLeaf<T2>::make(t), CreateLeaf<T3>::make(f)));
}
#endif // PETE_EXPRESSION_OPERATORS

template<class T1, class C1, class RHS>
inline Matrix<T1, C1>& assign(Matrix<T1, C1>& lhs, const RHS& rhs)
{
  using Leaf_t = typename CreateLeaf<RHS>::Leaf_t;
  evaluate(lhs, OpAssign(), MakeReturn<Leaf_t>::make(CreateLeaf<RHS>::make(rhs)));
  return lhs;
}

template<class T1, class C1, class RHS>
inline Matrix<T1, C1>& operator+=(Matrix<T1, C1>& lhs, const RHS& rhs)
{
  using Leaf_t = typename CreateLeaf<RHS>::Leaf_t;
  evaluate(lhs, OpAddAssign(), MakeReturn<Leaf_t>::make(CreateLeaf<RHS>::make(rhs)));
  return lhs;
}

template<class T1, class C1, class RHS>
inline Matrix<T1, C1>& operator-=(Matrix<T1, C1>& lhs, const RHS& rhs)
{
  using Leaf_t = typename CreateLeaf<RHS>::Leaf_t;
  evaluate(lhs, OpSubtractAssign(), MakeReturn<Leaf_t>::make(CreateLeaf<RHS>::make(rhs)));
  return lhs;
}

template<class T1, class C1, class RHS>
inline Matrix<T1, C1>& operator*=(Matrix<T1, C1>& lhs, const RHS& rhs)
{
  using Leaf_t = typename CreateLeaf<RHS>::Leaf_t;
  evaluate(lhs, OpMultiplyAssign(), MakeReturn<Leaf_t>::make(CreateLeaf<RHS>::make(rhs)));
  return lhs;
}

template<class T1, class C1, class RHS>
inline Matrix<T1, C1>& operator%=(Matrix<T1, C1>& lhs, const RHS& rhs)
{
  using Leaf_t = typename CreateLeaf<RHS>::Leaf_t;
  evaluate(lhs, OpModAssign(), MakeReturn<Leaf_t>::make(CreateLeaf<RHS>::make(rhs)));
  return lhs;
}

template<class T1, class C1, class RHS>
inline Matrix<T1, C1>& operator|=(Matrix<T1, C1>& lhs, const RHS& rhs)
{
  using Leaf_t = typename CreateLeaf<RHS>::Leaf_t;
  evaluate(lhs, OpBitwiseOrAssign(), MakeReturn<Leaf_t>::make(CreateLeaf<RHS>::make(rhs)));
  return lhs;
}

template<class T1, class C1, class RHS>
inline Matrix<T1, C1>& operator&=(Matrix<T1, C1>& lhs, const RHS& rhs)
{
  using Leaf_t = typename CreateLeaf<RHS>::Leaf_t;
  evaluate(lhs, OpBitwiseAndAssign(), MakeReturn<Leaf_t>::make(CreateLeaf<RHS>::make(rhs)));
  return lhs;
}

template<class T1, class C1, class RHS>
inline Matrix<T1, C1>& operator^=(Matrix<T1, C1>& lhs, const RHS& rhs)
{
  using Leaf_t = typename CreateLeaf<RHS>::Leaf_t;
  evaluate(lhs, OpBitwiseXorAssign(), MakeReturn<Leaf_t>::make(CreateLeaf<RHS>::make(rhs)));
  return lhs;
}

} // namespace qmcplusplus

#endif // OHMMS_MATRIXOPERATOR_H

/*
MatrixDefs.in
classes
-----
  ARG   = "class T[n], class C[n]"
  CLASS = "Matrix<T[n], C[n]>"


unaryOps
-----
  TAG = "OpUnaryMinus"
  FUNCTION = "operator-"
  EXPR = "return (-a);"
-----
  TAG = "OpUnaryPlus"
  FUNCTION = "operator+"
  EXPR = "return (+a);"
-----
  TAG = "OpBitwiseNot"
  FUNCTION = "operator~"
  EXPR = "return (~a);"
-----
  TAG = "OpIdentity"
  FUNCTION = "PETE_identity"
  EXPR = "return (a);"

unaryCastOps
-----
  TAG = "OpCast"
  FUNCTION = "peteCast"
  EXPR = "return T1(a);"

binaryOps
-----
  TAG = "OpAdd"
  FUNCTION = "operator+"
  EXPR = "return (a + b);"
-----
  TAG = "OpSubtract"
  FUNCTION = "operator-"
  EXPR = "return (a - b);"
-----
  TAG = "OpMultiply"
  FUNCTION = "operator*"
  EXPR = "return (a * b);"
-----
  TAG = "OpMod"
  FUNCTION = "operator%"
  EXPR = "return (a % b);"
-----
  TAG = "OpBitwiseAnd"
  FUNCTION = "operator&"
  EXPR = "return (a & b);"
-----
  TAG = "OpBitwiseOr"
  FUNCTION = "operator|"
  EXPR = "return (a | b);"
-----
  TAG = "OpBitwiseXor"
  FUNCTION = "operator^"
  EXPR = "return (a ^ b);"

assignOp
-----
  TAG = "OpAssign"
  FUNCTION = "assign"
  EXPR = "return (const_cast<T1 &>(a) = b);"

binaryAssignOps
-----
  TAG = "OpAddAssign"
  FUNCTION = "operator+="
  EXPR = "(const_cast<T1 &>(a) += b); return const_cast<T1 &>(a);"
-----
  TAG = "OpSubtractAssign"
  FUNCTION = "operator-="
  EXPR = "(const_cast<T1 &>(a) -= b); return const_cast<T1 &>(a);"
-----
  TAG = "OpMultiplyAssign"
  FUNCTION = "operator*="
  EXPR = "(const_cast<T1 &>(a) *= b); return const_cast<T1 &>(a);"
-----
  TAG = "OpModAssign"
  FUNCTION = "operator%="
  EXPR = "(const_cast<T1 &>(a) %= b); return const_cast<T1 &>(a);"
-----
  TAG = "OpBitwiseOrAssign"
  FUNCTION = "operator|="
  EXPR = "(const_cast<T1 &>(a) |= b); return const_cast<T1 &>(a);"
-----
  TAG = "OpBitwiseAndAssign"
  FUNCTION = "operator&="
  EXPR = "(const_cast<T1 &>(a) &= b); return const_cast<T1 &>(a);"
-----
  TAG = "OpBitwiseXorAssign"
  FUNCTION = "operator^="
  EXPR = "(const_cast<T1 &>(a) ^= b); return const_cast<T1 &>(a);"

trinaryOps
-----
  TAG = "FnWhere"
  FUNCTION = "where"
  EXPR = "if (a) return b; else return c;"
*/
