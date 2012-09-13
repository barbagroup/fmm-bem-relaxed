#pragma once

#include <iostream>
#include <math.h>

/** @file Point.hpp
 * @brief Define the Point class for 3D points. */

/** @class Point
 * @brief Class representing 3D points and vectors.
 *
 * Point contains methods that support use as points in 3D space, and
 * use as 3-dimensional vectors (that can, for example, be cross-producted).
 *
 * Its x, y, and z components are publicly accessible under those names. They
 * can also be accessed as coordinate[0], coordinate[1], and coordinate[2],
 * respectively.
 */
struct Point {
  union {
    struct {
      double x;
      double y;
      double z;
    };
    double coordinate[3];
  };

  // CONSTRUCTORS

  /** Construct a zero-valued point. */
  Point()
    : x(0), y(0), z(0) {
  }
  /** Construct a point with the given coordinates. */
  Point(double _x, double _y, double _z)
    : x(_x), y(_y), z(_z) {
  }


  // ACCESSORS

  /** Return this point's magnitude */
  double magnitude() const {
    return sqrt(squared_mag());
  }
  /** Return the square of this point's magnitude. */
  double squared_mag() const {
    return x*x + y*y + z*z;
  }


  // COMPARATORS

  /** Test whether two Points are equal.
   *
   * Points are equal if they have the same coordinates. */
  friend inline bool operator==(const Point& a, const Point& b) {
    return a.x == b.x && a.y == b.y && a.z == b.z;
  }

  /** Test whether @a a is less than @a b, lexically by coordinates.
   *
   * @a a < @a b if @a a.x < @a b.x, or (@a a.x == @a b.x && @a a.y < @a b.y),
   * or (@a a.x == @a b.x && @a a.y == @a b.y && @a a.z < @a b.z).
   *
   * This ordering is not meaningful geometrically, but can be useful when
   * storing Point objects in containers like map<>. */
  friend inline bool operator<(const Point& a, const Point& b) {
    return a.x < b.x || (a.x == b.x && a.y < b.y)
      || (a.x == b.x && a.y == b.y && a.z < b.z);
  }

  /** Other comparisons are defined in terms of < and ==. */
  friend inline bool operator!=(const Point& a, const Point& b) {
    return !(a == b);
  }
  friend inline bool operator<=(const Point& a, const Point& b) {
    return !(b < a);
  }
  friend inline bool operator>(const Point& a, const Point& b) {
    return b < a;
  }
  friend inline bool operator>=(const Point& a, const Point& b) {
    return !(a < b);
  }


  // MODIFIERS

  /** Add @a p to this point. */
  Point& operator+=(const Point& p) {
    x += p.x;
    y += p.y;
    z += p.z;
    return *this;
  }
  /** Return a negated version of @a p. */
  inline Point operator-() const {
    return Point(-x, -y, -z);
  }
  /** Subtract @a p from this point. */
  Point& operator-=(const Point& p) {
    return *this += -p;
  }
  /** Scale this point by a factor. */
  Point& operator*=(double d) {
    x *= d;
    y *= d;
    z *= d;
    return *this;
  }
  /** Scale this point by three factors. */
  Point& operator*=(const Point& p) {
    x *= p.x;
    y *= p.y;
    z *= p.z;
    return *this;
  }
  /** Scale this point down by a factor. */
  Point& operator/=(double d) {
    return *this *= (1 / d);
  }
  /** Scale this point down by three factors. */
  Point& operator/=(const Point& p) {
    x /= p.x;
    y /= p.y;
    z /= p.z;
    return *this;
  }
  /** Normalize this point: scale it so that its magnitude() is 1.
   * @post If the old magnitude() != 0, then the new magnitude() == 1.
   * @return this point
   *
   * A zero-valued point cannot be normalized and is left unchanged. */
  Point& normalize() {
    double l = squared_mag();
    if (l != 0)
      *this /= sqrt(l);
    return *this;
  }
  /** Return the cross product of this point with @a p. */
  Point cross(const Point& p) const {
    return Point(y * p.z - z * p.y, z * p.x - x * p.z, x * p.y - y * p.x);
  }


  /** Return the dot product of this point with @a p. */
  double dot(const Point& p) const {
    return x * p.x + y * p.y + z * p.z;
  }
};

/** Read a Point from an input stream
 *
 * The Point is serialized as three whitespace-separated real numbers
 * in "x y z" order.
 */
inline std::istream& operator>>(std::istream& s, Point& p) {
  return (s >> p.x >> p.y >> p.z);
}

/** Write a Point to an output stream.
 *
 * The Point is written as three double-precision numbers
 * in "x y z" order, separated by spaces.
 */
inline std::ostream& operator<<(std::ostream& s, const Point& p) {
  return (s << p.x << ' ' << p.y << ' ' << p.z);
}

// ARITHMETIC
/** Unary plus: Return @a p. ("+p" should work if "-p" works.) */
inline Point operator+(const Point &p) {
  return p;
}
inline Point operator+(Point a, const Point& b) {
  return a += b;
}
inline Point operator-(Point a, const Point& b) {
  return a -= b;
}
inline Point operator*(Point p, double d) {
  return p *= d;
}
inline Point operator*(Point a, const Point& b) {
  return a *= b;
}
inline Point operator/(Point p, double d) {
  return p /= d;
}
inline Point operator/(Point a, const Point& b) {
  return a /= b;
}
