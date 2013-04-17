#pragma once
/** @file Vec.hpp
 * @brief A custom constant size vector class based on a Boost uBlas vector.
 *
 * This provides convenient constructors to initialize coordinates,
 * imports common operators such as norms and inner products, and implements
 * additional expression templates such as scalar addition/subtration and
 * elementwise multiplication/division.
 */

#include <boost/numeric/ublas/vector.hpp>
namespace ublas = boost::numeric::ublas;

#include <type_traits>


template <unsigned N, typename T>
class Vec
    : public ublas::vector<T, ublas::bounded_array<T,N>>
{
  typedef ublas::vector<T, ublas::bounded_array<T,N>> super_type;

  /** Template unrolling for assigning an argument pack to this->data() */
  template <unsigned I>
  inline void insert() {}
  template <unsigned I, typename A, typename ...Rest>
  inline void insert(const A& a, Rest... r) {
    this->data()[I] = a;
    insert<I+1>(r...);
  }

  /** Template for determining if *all* types in a pack
   * are convertible to type @a To */
  template <typename To, typename ...>
  struct all_convertible;
  template <typename To>
  struct all_convertible<To>
      : std::true_type {};
  template <typename To, typename From, typename ...Rest>
  struct all_convertible<To, From, Rest...>
      : std::integral_constant<bool,
                               std::is_convertible<From,To>::value &&
                               all_convertible<To,Rest...>::value> {};
public:
  static constexpr unsigned dimension = N;

  /** Constructors */
  Vec()
      : super_type(N,T()) {}
  Vec(const Vec& v)
      : super_type(v) {}
  template <class A>
  Vec(const ublas::vector<T,A>& v)
      : super_type (v) {}
  template <class E>
  Vec(const ublas::vector_expression<E>& e)
      : super_type(e) {}

  /** Construct with value for all coordinates */
  explicit Vec(const T& t)
      : super_type(N,t) {
  }
  /** Construct with values for each coordinate */
  template <typename ...Arg,
            typename std::enable_if<sizeof...(Arg) == N,int>::type=0,
            typename std::enable_if<all_convertible<T,Arg...>::value,int>::type=0>
  explicit Vec(Arg ...args)
      : super_type(N) {
    insert<0>(args...);
  }

  /** Assignment operators */
  BOOST_UBLAS_INLINE
  Vec& operator=(const Vec& v) {
    super_type::operator=(v);
    return *this;
  }
  template <class A>
  BOOST_UBLAS_INLINE
  Vec& operator=(const ublas::vector<T,A>& v) {
    super_type::operator=(v);
    return *this;
  }
  template <class E>
  BOOST_UBLAS_INLINE
  Vec& operator=(const ublas::vector_expression<E>& v) {
    super_type::operator=(v);
    return *this;
  }
  /** Container assignment without temporary */
  template <class C>
  BOOST_UBLAS_INLINE
  Vec& operator=(const ublas::vector_container<C>& v) {
    super_type::operator=(v);
    return *this;
  }
};

// OPERATORS
#include <algorithm>
#include <iostream>

/** Equality comparison (weak) */
template <unsigned N, typename T>
BOOST_UBLAS_INLINE
bool operator==(const Vec<N,T>& a,
                const Vec<N,T>& b) {
  return std::equal(a.begin(), a.end(), b.begin());
}
/** Send to output stream */
template <unsigned N, typename T>
std::ostream& operator<<(std::ostream& s,
                         const Vec<N,T>& v) {
  s << "(";
  std::copy(v.begin(), v.end(), std::ostream_iterator<T>(s, ", "));
  return s << "\b\b)";
}
/** Inner product */
using ublas::inner_prod;
/** L1 norm */
using ublas::norm_1;
/** L2 norm */
using ublas::norm_2;
/** L_inf norm */
using ublas::norm_inf;

// TEMP: BACKWARD COMPATABILITY

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
