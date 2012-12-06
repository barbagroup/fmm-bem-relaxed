#pragma once
/** @file Vec.hpp
 * @brief A wrapper class to provide simple vector operations for
 * primitive types and classes that only require op[].
 */

#include <iostream>
#include <type_traits>
#include <cstdarg>
#include <math.h>
#include <assert.h>

//! Predeclaration of SmallVec
template <unsigned DIM, typename DATA>
class SmallVec;

/** @class Vec
 * A general class that endows types with vector operations. If the base type
 * is a primitive type, the Vec stores information as an array of that type.
 * If the base type is not a primitive type, then it must implement the bracket
 * operator (lvalue and rvalue)
 *
 * All of the following result in nearly the same public Vec interface
 * Vec<3, double> v1;
 * Vec<3, double[3]> v2;
 * Vec<3, MyClassThatImplementsOpBracketWithDoubles> v3;
 */
template <unsigned DIM, typename DATA>
using Vec = typename std::conditional<std::is_fundamental<DATA>::value,
                                      SmallVec<DIM, DATA[DIM]>,
                                      SmallVec<DIM, DATA>>::type;

#define for_i for(unsigned i=0; i!=dimension; ++i)

/** @class SmallVec
 * @brief Class representing ND points and vectors.
 *
 * SmallVec contains methods that support use of the underlying type
 * as points in ND space. The underlying data is stored as type POINT,
 * which must only provide the operators:
 * value_type& operator[](unsigned i)
 * const value_type& operator[](unsigned i) const
 */
template <unsigned DIM, typename POINT>
class SmallVec {
  static_assert(DIM > 0, "SmallVec DIM must be greater than 0");

  POINT a;

 public:
  static constexpr unsigned dimension = DIM;
  typedef POINT point_type;
  typedef typename std::remove_reference<decltype(a[0])>::type value_type;

  // CONSTRUCTORS

  SmallVec() {
    for_i a[i] = value_type(0);
  }
  SmallVec(const SmallVec& b) {
    for_i a[i] = b[i];
  }
  explicit SmallVec(const point_type& b) {
    for_i a[i] = b[i];
  }
  explicit SmallVec(value_type b) {
    for_i a[i] = b;
  }
  SmallVec(value_type b0, value_type b1) {
    static_assert(dimension == 2, "Calling 2-value Vec constructor for Non-2D Vec");
    a[0] = b0; a[1] = b1;
  }
  SmallVec(value_type b0, value_type b1, value_type b2) {
    static_assert(dimension == 3, "Calling 3-value Vec constructor for Non-3D Vec");
    a[0] = b0; a[1] = b1; a[2] = b2;
  }
  SmallVec(value_type b0, value_type b1, value_type b2, value_type b3) {
    static_assert(dimension == 4, "Calling 4-value Vec constructor for Non-4D Vec");
    a[0] = b0; a[1] = b1; a[2] = b2; a[3] = b3;
  }

  inline bool operator==(const SmallVec& b) const {
    for_i if (a[i] != b[i]) return false;
    return true;
  }

  // MODIFIERS

  /** Return a negated version of @a p. */
  inline SmallVec operator-() const {
    return SmallVec(-a[0], -a[1], -a[2]);
  }
  /** Scalar assignment */
  inline SmallVec& operator=(const value_type b) {
    for_i a[i] = b;
    return *this;
  }
  /** Add scalar @a b to this SmallVec */
  inline SmallVec& operator+=(const value_type b) {
    for_i a[i] += b;
    return *this;
  }
  /** Subtract scalar @a b from this SmallVec */
  inline SmallVec& operator-=(const value_type b) {
    for_i a[i] -= b;
    return *this;
  }
  /** Scale this SmallVec up by scalar @a b */
  inline SmallVec& operator*=(const value_type b) {
    for_i a[i] *= b;
    return *this;
  }
  /** Scale this SmallVec down by scalar @a b */
  inline SmallVec& operator/=(const value_type b) {
    for_i a[i] /= b;
    return *this;
  }
  /** Add SmallVec @a b to this SmallVec */
  inline SmallVec& operator+=(const SmallVec& b) {
    for_i a[i] += b[i];
    return *this;
  }
  /** Subtract SmallVec @a b from this SmallVec */
  inline SmallVec& operator-=(const SmallVec& b) {
    for_i a[i] -= b[i];
    return *this;
  }
  /** Scale this SmallVec up by factors in @a b */
  inline SmallVec& operator*=(const SmallVec& b) {
    for_i a[i] *= b[i];
    return *this;
  }
  /** Scale this SmallVec down by factors in @a b */
  inline SmallVec& operator/=(const SmallVec& b) {
    for_i a[i] /= b[i];
    return *this;
  }

  // ACCESSORS

  /** Access the @a i th (lvalue) element of this SmallVec
   * @pre i < dimension */
  inline value_type& operator[](unsigned i) {
    //assert(i < dimension);
    return a[i];
  }
  /** Access the @a i th (rvalue) element of this SmallVec
   * @pre i < dimension */
  inline const value_type& operator[](unsigned i) const {
    //assert(i < dimension);
    return a[i];
  }

  /** Compute the dot product of this SmallVec with another SmallVec */
  inline value_type dot(const SmallVec& b) const {
    value_type d(0);
    for_i d += a[i]*b[i];
    return d;
  }
  /** Compute the dot product of this SmallVec with another SmallVec */
  inline friend value_type dot(const SmallVec& a, const SmallVec& b) {
    return a.dot(b);
  }
  /** Compute cross product of two SmallVecs */
  inline friend SmallVec cross(const SmallVec& a, const SmallVec& b) {
    static_assert(SmallVec::dimension == 3, "Cross product only defined for 3D");
    SmallVec r;
    r[0] = a[1]*b[2] - a[2]*b[1];
    r[1] = -(a[0]*b[2] - a[2]*b[0]);
    r[2] = a[0]*b[1] - a[1]*b[0];
    return r;
  }
  /** Compute the squared L2 norm of this SmallVec */
  inline friend value_type normSq(const SmallVec& b) {
    return b.dot(b);
  }
  /** Compute the L2 norm of this SmallVec */
  inline friend value_type norm(const SmallVec& b) {
    return sqrt(normSq(b));
  }
  /** Write a SmallVec to an output stream */
  inline friend std::ostream& operator<<(std::ostream& s, const SmallVec& a) {
    for_i s << a[i] << " ";
    return s;
  }
};

#undef for_i

// ARITHMETIC

/** Unary plus: Return @a p. ("+p" should work if "-p" works.) */
template <unsigned D, typename P>
inline SmallVec<D,P> operator+(const SmallVec<D,P>& p) {
  return p;
}
template <unsigned D, typename P>
inline SmallVec<D,P> operator+(SmallVec<D,P> a, const SmallVec<D,P>& b) {
  return a += b;
}
template <unsigned D, typename P>
inline SmallVec<D,P> operator-(SmallVec<D,P> a, const SmallVec<D,P>& b) {
  return a -= b;
}
template <unsigned D, typename P>
inline SmallVec<D,P> operator*(SmallVec<D,P> p, double d) {
  return p *= d;
}
template <unsigned D, typename P>
inline SmallVec<D,P> operator*(double b, SmallVec<D,P> a) {
  return a *= b;
}
template <unsigned D, typename P>
inline SmallVec<D,P> operator*(SmallVec<D,P> a, const SmallVec<D,P>& b) {
  return a *= b;
}
template <unsigned D, typename P>
inline SmallVec<D,P> operator/(SmallVec<D,P> p, double d) {
  return p /= d;
}
template <unsigned D, typename P>
inline SmallVec<D,P> operator/(SmallVec<D,P> a, const SmallVec<D,P>& b) {
  return a /= b;
}
