#pragma once
/** @file L2P.hpp
 * @brief Dispatch methods for the L2P stage
 *
 */

struct L2P
{
  template <bool> struct UseVectorL2P {};

  /** Use *Substitution Failure Is Not An Error* and variadic templates
   * to determine if a Kernel has a L2P method with a certain signature
   */
  template <typename K, typename... Args>
  struct HasL2P {
    template <class A, void (A::*)(Args...) const> struct SFINAE {};
    template <class A> static std::true_type  sfinae(SFINAE<A, &A::L2P>*);
    template <class A> static std::false_type sfinae(...);
    static constexpr bool value = decltype(sfinae<K>(0))::value;
  };

  /** L2P evaluation.
   * The Kernel provides a vector L2P accumulator.
   */
  template <typename Kernel, typename TargetIter, typename ResultIter>
  inline static void eval(Kernel& K,
                          const typename Kernel::local_type& L,
                          const typename Kernel::point_type& center,
                          TargetIter t_begin, TargetIter t_end,
			  ResultIter r_begin,
                          UseVectorL2P<true>) {
    K.L2P(L, center, t_begin, t_end, r_begin);
  }

  /** L2P evaluation.
   * The Kernel provides a scalar L2P accumulator. Use it for each target point.
   */
  template <typename Kernel, typename TargetIter, typename ResultIter>
  inline static void eval(Kernel& K,
                          const typename Kernel::local_type& L,
                          const typename Kernel::point_type& center,
                          TargetIter t_begin, TargetIter t_end,
			  ResultIter r_begin,
                          UseVectorL2P<false>) {
    for ( ; t_begin != t_end; ++t_begin, ++r_begin)
      K.L2P(L, center, *t_begin, *r_begin);
  }

  /** L2P evaluation dispath.
   * Detects if Kernel has a vectorized L2P and uses it if available.
   */
  template <typename Kernel, typename TargetIter, typename ResultIter>
  inline static void eval(Kernel& K,
                          const typename Kernel::local_type& L,
                          const typename Kernel::point_type& center,
                          TargetIter t_begin, TargetIter t_end,
			  ResultIter r_begin) {
    typedef HasL2P<Kernel,
                   const typename Kernel::local_type&,
                   const typename Kernel::point_type&,
                   TargetIter, TargetIter, ResultIter> HasVectorL2P;

    L2P::eval(K, L, center, t_begin, t_end, r_begin,
              UseVectorL2P<HasVectorL2P::value>());
  }

 public:

  /** Unwrap the data from BoxContext and dispatch to the L2P evaluator
   */
  template <typename Kernel, typename BoxContext, typename Box>
  inline static void eval(Kernel& K,
			  BoxContext& bc,
			  const Box& box)
  {
#ifdef DEBUG
    printf("L2P: %d\n", box.index());
#endif

    L2P::eval(K,
              bc.local_expansion(box), bc.center(box),
              bc.target_begin(box), bc.target_end(box),
              bc.result_begin(box));
  }
};
