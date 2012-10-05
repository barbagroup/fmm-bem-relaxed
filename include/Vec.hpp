#pragma once
/** @file Vec.hpp
 * @brief A wrapper class to provide simple vector operators for
 * a class that only requires the operator[]
 */

#include <iostream>
#include <cstdarg>
#include <math.h>

#define for_i for(unsigned i=0; i!=dimension; ++i)

/** @class Vec
 * @brief Class representing ND points and vectors.
 *
 * Vec contains methods that supports use as points in ND space.
 */
template <unsigned DIM, typename POINT>
class Vec {
 private:
  POINT a;

 public:
  static constexpr unsigned dimension = DIM;
  typedef POINT point_type;
  typedef typename std::remove_reference<decltype(a[0])>::type value_type;

  // CONSTRUCTORS

  Vec() {
    for_i a[i] = value_type(0);
  }
  Vec(const Vec& b) {
    for_i a[i] = b[i];
  }
  explicit Vec(const point_type& b) {
    for_i a[i] = b[i];
  }
  explicit Vec(value_type b) {
    for_i a[i] = b;
  }
  Vec(value_type b0, value_type b1) {
    static_assert(dimension == 3, "Calling 2-value Vec constructor for Non-2D Vec");
    a[0] = b0; a[1] = b1;
  }
  Vec(value_type b0, value_type b1, value_type b2) {
    static_assert(dimension == 3, "Calling 3-value Vec constructor for Non-3D Vec");
    a[0] = b0; a[1] = b1; a[2] = b2;
  }
  Vec(value_type b0, value_type b1, value_type b2, value_type b3) {
    static_assert(dimension == 4, "Calling 4-value Vec constructor for Non-4D Vec");
    a[0] = b0; a[1] = b1; a[2] = b2; a[3] = b3;
  }

  // MODIFIERS

  /** Return a negated version of @a p. */
  inline Vec operator-() const {
    return Vec(-a[0], -a[1], -a[2]);
  }
  /** Scalar assignment */
  Vec& operator=(const value_type b) {
    for_i a[i] = b;
    return *this;
  }
  /** Add scalar @a b to this Vec */
  Vec& operator+=(const value_type b) {
    for_i a[i] += b;
    return *this;
  }
  /** Subtract scalar @a b from this Vec */
  Vec& operator-=(const value_type b) {
    for_i a[i] -= b;
    return *this;
  }
  /** Scale this Vec up by scalar @a b */
  Vec& operator*=(const value_type b) {
    for_i a[i] *= b;
    return *this;
  }
  /** Scale this Vec down by scalar @a b */
  Vec& operator/=(const value_type b) {
    for_i a[i] /= b;
    return *this;
  }
  /** Add Vec @a b to this Vec */
  Vec& operator+=(const Vec& b) {
    for_i a[i] += b[i];
    return *this;
  }
  /** Subtract Vec @a b from this Vec */
  Vec& operator-=(const Vec& b) {
    for_i a[i] -= b[i];
    return *this;
  }
  /** Scale this Vec up by factors in @a b */
  Vec& operator*=(const Vec& b) {
    for_i a[i] *= b[i];
    return *this;
  }
  /** Scale this Vec down by factors in @a b */
  Vec& operator/=(const Vec& b) {
    for_i a[i] /= b[i];
    return *this;
  }

  // ACCESSORS

  /** Access the @a i th (lvalue) element of this Vec
   * @pre i < dimension */
  value_type& operator[](int i) {
    return a[i];
  }
  /** Access the @a i th (rvalue) element of this Vec
   * @pre i < dimension */
  const value_type& operator[](int i) const {
    return a[i];
  }

  /** Compute the squared L2 norm of this Vec */
  friend value_type norm(const Vec& b) {
    value_type c(0);
    for_i c += b[i]*b[i];
    return c;
  }
  /** Write a Vec to an output stream */
  friend std::ostream& operator<<(std::ostream& s, const Vec& a) {
    for_i s << a[i] << " ";
    return s;
  }
};

#undef for_i

// ARITHMETIC

/** Unary plus: Return @a p. ("+p" should work if "-p" works.) */
template <unsigned D, typename P>
inline Vec<D,P> operator+(const Vec<D,P>& p) {
  return p;
}
template <unsigned D, typename P>
inline Vec<D,P> operator+(Vec<D,P> a, const Vec<D,P>& b) {
  return a += b;
}
template <unsigned D, typename P>
inline Vec<D,P> operator-(Vec<D,P> a, const Vec<D,P>& b) {
  return a -= b;
}
template <unsigned D, typename P>
inline Vec<D,P> operator*(Vec<D,P> p, double d) {
  return p *= d;
}
template <unsigned D, typename P>
inline Vec<D,P> operator*(Vec<D,P> a, const Vec<D,P>& b) {
  return a *= b;
}
template <unsigned D, typename P>
inline Vec<D,P> operator/(Vec<D,P> p, double d) {
  return p /= d;
}
template <unsigned D, typename P>
inline Vec<D,P> operator/(Vec<D,P> a, const Vec<D,P>& b) {
  return a /= b;
}
