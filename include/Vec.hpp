#pragma once

#include <iostream>
#include <math.h>
#include <stdarg.h>

// TODO: Remove!
#include "vec.h"

/** @file Vec.hpp
 * @brief A wrapper class to provide simple vector operators for
 * a class that only requires the [const] operator[]
 */

#define for_i for(unsigned i=0; i!=dimension; ++i)

/** @class Vec
 * @brief Class representing ND points and vectors.
 *
 * Vec contains methods that support use as points in ND space.
 */
template <unsigned dim, typename point>
class Vec;

/** @class Vec
 * @brief Class representing 3D points and vectors.
 *
 * Vec contains methods that support use as points in 3D space.
 */
template <typename point>
class Vec<3, point> {
 private:
  point a;

 public:
  static constexpr unsigned dimension = 3;
  typedef point point_type;
  typedef typename std::remove_reference<decltype(a[0])>::type value_type;

  // CONSTRUCTORS

  Vec() {                                                          // Default constructor
    for_i a[i] = value_type(0);
  }
  Vec(const point_type& b) {                                       // Copy constructor (vector)
    for_i a[i] = b[i];
  }
  Vec(const Vec& b) {                                              // Copy constructor (vector)
    for_i a[i] = b[i];
  }
  Vec(value_type b) {                                             // Copy constructor (scalar)
    for_i a[i] = b;
  }
  Vec(value_type b0, value_type b1, value_type b2) {
    a[0] = b0; a[1] = b1; a[2] = b2;
  }


  // TODO: Remove!
  template <int N, typename T>
  Vec(const vec<N,T>& v) {
    for_i a[i] = v[i];
  }

  // MODIFIERS

  const Vec& operator=(const value_type b) {                       // Scalar assignment
    for_i a[i] = b;
    return *this;
  }
  template <typename T>
  const Vec& operator+=(const T b) {                               // Scalar compound assignment (add)
    for_i a[i] += b;
    return *this;
  }
  template <typename T>
  const Vec& operator-=(const T b) {                               // Scalar compound assignment (subtract)
    for_i a[i] -= b;
    return *this;
  }
  template <typename T>
  const Vec& operator*=(const T b) {                               // Scalar compound assignment (multiply)
    for_i a[i] *= b;
    return *this;
  }
  template <typename T>
  const Vec& operator/=(const T b) {                               // Scalar compound assignment (divide)
    for_i a[i] /= b;
    return *this;
  }
  const Vec& operator=(const Vec& b) {                             // Vector assignment
    for_i a[i] = b[i];
    return *this;
  }
  const Vec& operator+=(const Vec& b) {                            // Vector compound assignment (add)
    for_i a[i] += b[i];
    return *this;
  }
  const Vec& operator-=(const Vec& b) {                            // Vector compound assignment (subtract)
    for_i a[i] -= b[i];
    return *this;
  }
  const Vec& operator*=(const Vec& b) {                            // Vector compound assignment (multiply)
    for_i a[i] *= b[i];
    return *this;
  }
  const Vec& operator/=(const Vec& b) {                            // Vector compound assignment (divide)
    for_i a[i] /= b[i];
    return *this;
  }
  template <typename T>
  Vec operator+(const T b) const {                                 // Scalar arithmetic (add)
    Vec c;
    for_i c[i] = a[i] + b;
    return c;
  }
  template <typename T>
  Vec operator-(const T b) const {                                 // Scalar arithmetic (subtract)
    Vec c;
    for_i c[i] = a[i] - b;
    return c;
  }
  template <typename T>
  Vec operator*(const T b) const {                                 // Scalar arithmetic (multiply)
    Vec c;
    for_i c[i] = a[i] * b;
    return c;
  }
  template <typename T>
  Vec operator/(const T b) const {                                 // Scalar arithmetic (divide)
    Vec c;
    for_i c[i] = a[i] / b;
    return c;
  }
  Vec operator+(const Vec& b) const {                              // Vector arithmetic (add)
    Vec c;
    for_i c[i] = a[i] + b[i];
    return c;
  }
  Vec operator-(const Vec& b) const {                              // Vector arithmetic (subtract)
    Vec c;
    for_i c[i] = a[i] - b[i];
    return c;
  }
  Vec operator*(const Vec& b) const {                              // Vector arithmetic (multiply)
    Vec c;
    for_i c[i] = a[i] * b[i];
    return c;
  }
  Vec operator/(const Vec& b) const {                              // Vector arithmetic (divide)
    Vec c;
    for_i c[i] = a[i] / b[i];
    return c;
  }
  value_type& operator[](int i) {                                  // Indexing (lvalue)
    return a[i];
  }
  const value_type& operator[](int i) const {                      // Indexing (rvalue)
    return a[i];
  }
  //operator       point_type()       { return a; }                  // Type-casting (lvalue)
  //operator const point_type() const { return a; }                  // Type-casting (rvalue)

  friend std::ostream& operator<<(std::ostream& s, const Vec& a) { // Component-wise output stream
    for_i s << a[i] << " ";
    return s;
  }
  friend auto norm(const Vec& b) -> decltype(a[0]) {                                    // L2 norm squared
    value_type c = value_type(0);
    for_i c += b[i]*b[i];
    return c;
  }
};

#undef for_i
