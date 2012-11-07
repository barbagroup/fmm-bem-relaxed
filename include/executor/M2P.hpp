#pragma once
/** @file M2P.hpp
 * @brief Dispatch methods for the M2P stage
 *
 */

struct M2P
{
  template <bool> struct UseVectorM2P {};

  /** Use *Substitution Failure Is Not An Error* and variadic templates
   * to determine if a Kernel has a M2P method with a certain signature
   */
  template <typename K, typename... Args>
  struct HasM2P {
    template <class A, void (A::*)(Args...) const> struct SFINAE {};
    template <class A> static std::true_type  sfinae(SFINAE<A, &A::M2P>*);
    template <class A> static std::false_type sfinae(...);
    static constexpr bool value = decltype(sfinae<K>(0))::value;
  };

  /** M2P evaluation.
   * The Kernel provides a vector M2P accumulator.
   */
  template <typename Kernel, typename PointIter, typename ResultIter>
  inline static void eval(Kernel& K,
                          const typename Kernel::multipole_type& M,
                          const typename Kernel::point_type& center,
                          PointIter p_begin, PointIter p_end, ResultIter r_begin,
                          UseVectorM2P<true>) {
    K.M2P(M, center, p_begin, p_end, r_begin);
  }

  /** M2P evaluation.
   * The Kernel provides a scalar M2P accumulator. Use it for each target point.
   */
  template <typename Kernel, typename PointIter, typename ResultIter>
  inline static void eval(Kernel& K,
                          const typename Kernel::multipole_type& M,
                          const typename Kernel::point_type& center,
                          PointIter p_begin, PointIter p_end, ResultIter r_begin,
                          UseVectorM2P<false>) {
    for ( ; p_begin != p_end; ++p_begin, ++r_begin)
      K.M2P(M, center, *p_begin, *r_begin);
  }

  /** M2P evaluation dispath.
   * Detects if Kernel has a vectorized M2P and uses it if available.
   */
  template <typename Kernel, typename PointIter, typename ResultIter>
  inline static void eval(Kernel& K,
                          const typename Kernel::multipole_type& M,
                          const typename Kernel::point_type& center,
                          PointIter p_begin, PointIter p_end, ResultIter r_begin) {
    typedef HasM2P<Kernel,
                   const typename Kernel::multipole_type&,
                   const typename Kernel::point_type&,
                   PointIter, PointIter, ResultIter> HasVectorM2P;

    M2P::eval(K, M, center, p_begin, p_end, r_begin,
              UseVectorM2P<HasVectorM2P::value>());
  }

 public:

  /** Unwrap the data from BoxContext and dispatch to the M2P evaluator
   */
  template <typename Kernel, typename BoxContext, typename Box>
  inline static void eval(Kernel& K,
			  BoxContext& bc,
			  const Box& source,
			  const Box& target) {
#ifdef DEBUG
    printf("M2P: %d to %d\n", source.index(), target.index());
#endif

    M2P::eval(K,
              bc.multipole_expansion(source), bc.center(source),
              bc.point_begin(target), bc.point_end(target),
              bc.result_begin(target));
  }
};
