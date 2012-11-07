#pragma once
/** @file P2M.hpp
 * @brief Dispatch methods for the P2M stage
 *
 */

struct P2M
{
  template <bool> struct UseVectorP2M {};

  /** Use *Substitution Failure Is Not An Error* and variadic templates
   * to determine if a Kernel has a P2M method with a certain signature
   */
  template <typename K, typename... Args>
  struct HasP2M {
    template <class A, void (A::*)(Args...) const> struct SFINAE {};
    template <class A> static constexpr void sfinae(SFINAE<A, &A::P2M>*);
    template <class A> static constexpr char sfinae(...);
    static constexpr bool value = std::is_void<decltype(sfinae<K>(0))>::value;
  };

  /** P2M evaluation.
   * The Kernel provides a vector P2M accumulator.
   */
  template <typename Kernel, typename SourceIter, typename ChargeIter>
  inline static void eval(Kernel& K,
                          SourceIter s_begin, SourceIter s_end,
			  ChargeIter c_begin,
                          const typename Kernel::point_type& center,
                          typename Kernel::multipole_type& M,
                          UseVectorP2M<true>) {
    K.P2M(s_begin, s_end, c_begin, center, M);
  }

  /** P2M evaluation.
   * The Kernel provides a scalar P2M accumulator. Use it for each source point.
   */
  template <typename Kernel, typename SourceIter, typename ChargeIter>
  inline static void eval(Kernel& K,
                          SourceIter s_begin, SourceIter s_end,
			  ChargeIter c_begin,
                          const typename Kernel::point_type& center,
                          typename Kernel::multipole_type& M,
                          UseVectorP2M<false>) {
    for ( ; s_begin != s_end; ++s_begin, ++c_begin)
      K.P2M(*s_begin, *c_begin, center, M);
  }

  /** P2M evaluation dispath.
   * Detects if Kernel has a vectorized P2M and uses it if available.
   */
  template <typename Kernel, typename SourceIter, typename ChargeIter>
  inline static void eval(Kernel& K,
                          SourceIter s_begin, SourceIter s_end,
			  ChargeIter c_begin,
                          const typename Kernel::point_type& center,
                          typename Kernel::multipole_type& M) {
    typedef HasP2M<Kernel,
                   SourceIter, SourceIter, ChargeIter,
                   const typename Kernel::point_type&,
                   typename Kernel::multipole_type&> HasVectorP2M;

    P2M::eval(K, s_begin, s_end, c_begin, center, M,
              UseVectorP2M<HasVectorP2M::value>());
  }

 public:

  /** Unwrap the data from the BoxContext and dispatch to P2M evaluator
   */
  template <typename Kernel, typename BoxContext, typename Box>
  inline static void eval(Kernel& K,
			  BoxContext& bc,
			  const Box& box)
  {
#ifdef DEBUG
    printf("P2M: %d\n", box.index());
#endif

    P2M::eval(K,
              bc.source_begin(box), bc.source_end(box), bc.charge_begin(box),
              bc.center(box), bc.multipole_expansion(box));
  }
};
