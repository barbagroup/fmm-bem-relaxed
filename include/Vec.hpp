#pragma once
/** @file Vec.hpp
 * @brief A wrapper class to provide simple vector operations for
 * primitive types and classes that only require op[].
 */

#include <boost/numeric/ublas/vector.hpp>
namespace ublas = boost::numeric::ublas;

template <unsigned N, typename T>
class Vec
    : public ublas::vector<T, ublas::bounded_array<T,N>>
{
  typedef ublas::vector<T, ublas::bounded_array<T,N>> super_type;
public:
  static constexpr unsigned dimension = N;

  /** Constructors */
  Vec()
      : super_type(N,T()) {
  }
  explicit Vec(const T& t)
      : super_type(N,t) {
  }
  explicit Vec(const T& t0, const T& t1)
      : super_type(N) {
    static_assert(N == 2, "Vec<2,T> constructor parameter number mismatch!");
    this->data()[0] = t0;
    this->data()[1] = t1;
  }
  explicit Vec(const T& t0, const T& t1, const T& t2)
      : super_type(N) {
    static_assert(N == 3, "Vec<3,T> constructor parameter number mismatch!");
    this->data()[0] = t0;
    this->data()[1] = t1;
    this->data()[2] = t2;
  }
  explicit Vec(const T& t0, const T& t1, const T& t2, const T& t3)
      : super_type(N) {
    static_assert(N == 4, "Vec<4,T> constructor parameter number mismatch!");
    this->data()[0] = t0;
    this->data()[1] = t1;
    this->data()[2] = t2;
    this->data()[3] = t3;
  }
  template <class E>
  Vec(const ublas::vector_expression<E>& e)
      : super_type(e) {
  }

  /** Assignment operators */
  inline Vec& operator=(const super_type& v) {
    super_type::operator=(v);
    return *this;
  }
  template <class E>
  inline Vec& operator=(const ublas::vector_expression<E>& v) {
    super_type::operator=(v);
    return *this;
  }
};

// OPERATORS
#include <algorithm>
#include <iostream>

/** Compute the dot product */
template <typename E1, typename E2>
inline auto dot(const ublas::vector_expression<E1>& a,
                const ublas::vector_expression<E2>& b)
    -> decltype(ublas::inner_prod(a,b)) {
  return ublas::inner_prod(a,b);
}
/** Compute the squared L2 norm */
template <typename E>
inline auto normSq(const ublas::vector_expression<E>& a)
    -> decltype(ublas::inner_prod(a,a)) {
  return ublas::inner_prod(a,a);
}
/** Compute the L2 norm */
template <typename E>
inline auto norm(const ublas::vector_expression<E>& a)
    -> decltype(ublas::norm_2(a)) {
  return ublas::norm_2(a);
}
/** Equality comparison */
template <typename T, typename A>
inline bool operator==(const ublas::vector<T,A>& a,
                       const ublas::vector<T,A>& b) {
  return std::equal(a.begin(), a.end(), b.begin());
}
/** Send to output stream */
template <unsigned N, typename T>
std::ostream& operator<<(std::ostream& s,
                         const Vec<N,T>& v) {
  s << "(";
  std::copy(v.begin(), v.end(), std::ostream_iterator<double>(s, ", "));
  return s << "\b\b)";
}

/////////////////////////////////
// Vector Expression Operators //
/////////////////////////////////


/** Scalar addition */
template <class T1, class E2>
BOOST_UBLAS_INLINE
typename boost::enable_if<
  boost::is_convertible<T1,typename E2::value_type>,
  typename ublas::vector_binary_scalar1_traits<
    const T1, E2, ublas::scalar_plus<T1, typename E2::value_type>
    >::result_type
  >::type
operator+(const T1& e1,
          const ublas::vector_expression<E2>& e2) {
  typedef typename ublas::vector_binary_scalar1_traits<
    const T1, E2, ublas::scalar_plus<T1, typename E2::value_type>
    >::expression_type expression_type;
  return expression_type(e1, e2());
}

/** Scalar addition */
template <class E1, class T2>
BOOST_UBLAS_INLINE
typename boost::enable_if<
  boost::is_convertible<typename E1::value_type,T2>,
  typename ublas::vector_binary_scalar2_traits<
    E1, const T2, ublas::scalar_plus<typename E1::value_type,T2>
    >::result_type
  >::type
operator+(const ublas::vector_expression<E1>& e1,
          const T2& e2) {
  typedef typename ublas::vector_binary_scalar2_traits<
    E1, const T2, ublas::scalar_plus<typename E1::value_type,T2>
    >::expression_type expression_type;
  return expression_type(e1(), e2);
}

/** Scalar subtraction */
template <class T1, class E2>
BOOST_UBLAS_INLINE
typename boost::enable_if<
  boost::is_convertible<T1,typename E2::value_type>,
  typename ublas::vector_binary_scalar1_traits<
    const T1, E2, ublas::scalar_minus<T1, typename E2::value_type>
    >::result_type
  >::type
operator-(const T1& e1,
          const ublas::vector_expression<E2>& e2) {
  typedef typename ublas::vector_binary_scalar1_traits<
    const T1, E2, ublas::scalar_minus<T1, typename E2::value_type>
    >::expression_type expression_type;
  return expression_type(e1, e2());
}

/** Scalar subtraction */
template <class E1, class T2>
BOOST_UBLAS_INLINE
typename boost::enable_if<
  boost::is_convertible<typename E1::value_type,T2>,
  typename ublas::vector_binary_scalar2_traits<
    E1, const T2, ublas::scalar_minus<typename E1::value_type,T2>
    >::result_type
  >::type
operator-(const ublas::vector_expression<E1>& e1,
          const T2& e2) {
  typedef typename ublas::vector_binary_scalar2_traits<
    E1, const T2, ublas::scalar_minus<typename E1::value_type,T2>
    >::expression_type expression_type;
  return expression_type(e1(), e2);
}

/** Elementwise division */
template <class E1, class E2>
BOOST_UBLAS_INLINE
typename boost::enable_if<
  boost::is_convertible<typename E1::value_type,typename E2::value_type>,
  typename ublas::vector_binary_traits<
    E1, E2, ublas::scalar_divides<typename E1::value_type,
                                  typename E2::value_type>
    >::result_type
  >::type
operator/(const ublas::vector_expression<E1>& e1,
          const ublas::vector_expression<E2>& e2) {
  typedef typename ublas::vector_binary_traits<
    E1, E2, ublas::scalar_divides<typename E1::value_type,
                                  typename E2::value_type>
    >::expression_type expression_type;
  return expression_type(e1(), e2());
}

/** Elementwise multiplication */
template <class E1, class E2>
BOOST_UBLAS_INLINE
typename boost::enable_if<
  boost::is_convertible<typename E1::value_type,typename E2::value_type>,
  typename ublas::vector_binary_traits<
    E1, E2, ublas::scalar_multiplies<typename E1::value_type,
                                     typename E2::value_type>
    >::result_type
  >::type
operator*(const ublas::vector_expression<E1>& e1,
          const ublas::vector_expression<E2>& e2) {
  typedef typename ublas::vector_binary_traits<
    E1, E2, ublas::scalar_multiplies<typename E1::value_type,
                                     typename E2::value_type>
    >::expression_type expression_type;
  return expression_type(e1(), e2());
}
