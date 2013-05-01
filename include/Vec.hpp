#pragma once
/** @file Vec.hpp
 * @brief A custom fixed size vector class based on a Boost uBlas vector.
 *
 * This provides convenient constructors to initialize coordinates,
 * imports common operators such as norms and inner products, and implements
 * additional expression templates such as scalar addition/subtration and
 * elementwise multiplication/division.
 */

#define BOOST_UBLAS_NDEBUG
#include <boost/numeric/ublas/vector.hpp>
namespace ublas = boost::numeric::ublas;

#include <type_traits>

// Fixed array
template <class T, std::size_t N, class ALLOC = std::allocator<T>>
class fixed_array
    : public ublas::storage_array<fixed_array<T, N, ALLOC>>
{
  typedef fixed_array<T, N, ALLOC> self_type;
 public:
  // No allocator_type as ALLOC is not used for allocation
  typedef typename ALLOC::size_type size_type;
  typedef typename ALLOC::difference_type difference_type;
  typedef T value_type;
  typedef const T& const_reference;
  typedef T& reference;
  typedef const T* const_pointer;
  typedef T* pointer;
  typedef const_pointer const_iterator;
  typedef pointer iterator;

  // Construction and destruction
  BOOST_UBLAS_INLINE
  fixed_array() {
  }
  explicit BOOST_UBLAS_INLINE
  fixed_array(size_type size) {
    (void) size;
    BOOST_UBLAS_CHECK(size == N, ublas::bad_size());
    // data_ (an array) elements are already default constructed
  }
  BOOST_UBLAS_INLINE
  fixed_array(size_type size, const value_type& init) {
    (void) size;
    BOOST_UBLAS_CHECK(size == N, ublas::bad_size());
    // ISSUE elements should be value constructed here, but we must fill instead as already default constructed
    std::fill(begin(), end(), init) ;
  }
  BOOST_UBLAS_INLINE
  fixed_array(const fixed_array &c) {
    // ISSUE elements should be copy constructed here, but we must copy instead as already default constructed
    std::copy(c.begin(), c.end(), begin());
  }

  // Resizing
  BOOST_UBLAS_INLINE
  void resize(size_type size) {
    (void) size;
    BOOST_UBLAS_CHECK(size == N, ublas::bad_size());
  }
  BOOST_UBLAS_INLINE
  void resize(size_type size, value_type init) {
    (void) size;
    BOOST_UBLAS_CHECK(size == N, ublas::bad_size());
  }

  // Random Access Container
  BOOST_UBLAS_INLINE
  size_type max_size() const {
    return N;
  }

  BOOST_UBLAS_INLINE
  bool empty() const {
    return N == 0;
  }

  BOOST_UBLAS_INLINE
  size_type size() const {
    return N;
  }

  // Element access
  BOOST_UBLAS_INLINE
  const_reference operator[](size_type i) const {
    BOOST_UBLAS_CHECK (i < N, ublas::bad_index ());
    return data_[i];
  }
  BOOST_UBLAS_INLINE
  reference operator[](size_type i) {
    BOOST_UBLAS_CHECK (i < N, ublas::bad_index ());
    return data_[i];
  }

  // Assignment
  BOOST_UBLAS_INLINE
  fixed_array& operator=(const fixed_array& a) {
    if (this != &a) {
      std::copy(a.data_, a.data_ + N, data_);
    }
    return *this;
  }
  BOOST_UBLAS_INLINE
  fixed_array& assign_temporary(fixed_array& a) {
    *this = a;
    return *this;
  }

  // Swapping
  BOOST_UBLAS_INLINE
  void swap(fixed_array& a) {
    if (this != &a) {
      std::swap_ranges(data_, data_ + N, a.data_);
    }
  }
  BOOST_UBLAS_INLINE
  friend void swap(fixed_array& a1, fixed_array& a2) {
    a1.swap(a2);
  }

  BOOST_UBLAS_INLINE
  const_iterator begin() const {
    return data_;
  }
  BOOST_UBLAS_INLINE
  const_iterator end() const {
    return data_ + N;
  }

  BOOST_UBLAS_INLINE
  iterator begin() {
    return data_;
  }
  BOOST_UBLAS_INLINE
  iterator end() {
    return data_ + N;
  }

  // Reverse iterators
  typedef std::reverse_iterator<const_iterator> const_reverse_iterator;
  typedef std::reverse_iterator<iterator> reverse_iterator;

  BOOST_UBLAS_INLINE
  const_reverse_iterator rbegin() const {
    return const_reverse_iterator(end());
  }
  BOOST_UBLAS_INLINE
  const_reverse_iterator rend() const {
    return const_reverse_iterator(begin());
  }
  BOOST_UBLAS_INLINE
  reverse_iterator rbegin() {
    return reverse_iterator(end());
  }
  BOOST_UBLAS_INLINE
  reverse_iterator rend() {
    return reverse_iterator(begin());
  }

 private:
  value_type data_[N];
};


// --------------------
// Fixed vector class
// --------------------

/// @brief A dense vector of values of type @c T, of fixed size @f$N@f$.
/// A dense vector of values of type @c T, of fixed size @f$N@f$.
/// Elements are constructed by the storage type @c fixed_array
template<class T, std::size_t N>
class fixed_vector
    : public ublas::vector<T, fixed_array<T,N>>
{
  typedef ublas::vector<T, fixed_array<T,N>> vector_type;

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
  typedef typename vector_type::size_type size_type;
  static const size_type max_size = N;
  static constexpr size_type dimension = N;  // TEMP

  // Construction and destruction
  BOOST_UBLAS_INLINE
  fixed_vector()
      : vector_type(N,T()) {}
  BOOST_UBLAS_INLINE
  explicit fixed_vector(size_type size)
      : vector_type(N,T()) {
    (void) size;
    BOOST_UBLAS_CHECK(size == N, ublas::bad_size());
  }
  BOOST_UBLAS_INLINE
  fixed_vector(const fixed_vector& v)
      : vector_type(v) {}
  template<class A2>              // Allow vector<T,fixed_array<N> construction
  BOOST_UBLAS_INLINE
  fixed_vector(const ublas::vector<T, A2>& v)
      : vector_type(v) {}
  template<class AE>
  BOOST_UBLAS_INLINE
  fixed_vector(const ublas::vector_expression<AE>& ae)
      : vector_type(ae) {
  }
  BOOST_UBLAS_INLINE
  ~fixed_vector() {}

  /** Construct with value for all coordinates */
  //BOOST_UBLAS_INLINE
  //explicit fixed_vector(const T& t)
  //    : vector_type(N,t) {}
  /** Construct with values for each coordinate */
  template <typename ...Arg,
            typename std::enable_if<sizeof...(Arg) == N,int>::type=0,
            typename std::enable_if<all_convertible<T,Arg...>::value,int>::type=0>
  BOOST_UBLAS_INLINE
  explicit fixed_vector(Arg ...args)
      : vector_type(N) {
    insert<0>(args...);
  }

  // Assignment

  /*! @note "pass by value" the key idea to enable move semantics */
  BOOST_UBLAS_INLINE
  fixed_vector& operator=(const fixed_vector& v) {
    vector_type::operator=(v);
    return *this;
  }
  template<class A2>         // Generic vector assignment
  BOOST_UBLAS_INLINE
  fixed_vector& operator=(const ublas::vector<T, A2>& v) {
    vector_type::operator=(v);
    return *this;
  }
  template<class C>          // Container assignment without temporary
  BOOST_UBLAS_INLINE
  fixed_vector& operator=(const ublas::vector_container<C>& v) {
    vector_type::operator=(v);
    return *this;
  }
  template<class AE>
  BOOST_UBLAS_INLINE
  fixed_vector& operator=(const ublas::vector_expression<AE>& ae) {
    vector_type::operator=(ae);
    return *this;
  }
};


template <std::size_t N, typename T>
using Vec = fixed_vector<T,N>;


// OPERATORS
#include <algorithm>
#include <iostream>

/** Equality comparison (weak) */
template <std::size_t N, typename T>
BOOST_UBLAS_INLINE
bool operator==(const Vec<N,T>& a,
                const Vec<N,T>& b) {
  return std::equal(a.begin(), a.end(), b.begin());
}
template <std::size_t N, typename T>
BOOST_UBLAS_INLINE
bool operator!=(const Vec<N,T>& a,
                const Vec<N,T>& b) {
  return !(a == b);
}
/** Send to output stream */
template <std::size_t N, typename T>
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
