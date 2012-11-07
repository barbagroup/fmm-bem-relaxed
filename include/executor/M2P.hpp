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
    template <class A> static constexpr void sfinae(SFINAE<A, &A::M2P>*);
    template <class A> static constexpr char sfinae(...);
    static constexpr bool value = std::is_void<decltype(sfinae<K>(0))>::value;
  };

  /** M2P evaluation.
   * The Kernel provides a vector M2P accumulator.
   */
  template <typename Kernel, typename TargetIter, typename ResultIter>
  inline static void eval(Kernel& K,
                          const typename Kernel::multipole_type& M,
                          const typename Kernel::point_type& center,
                          TargetIter t_begin, TargetIter t_end,
			  ResultIter r_begin,
                          UseVectorM2P<true>) {
    K.M2P(M, center, t_begin, t_end, r_begin);
  }

  /** M2P evaluation.
   * The Kernel provides a scalar M2P accumulator. Use it for each target point.
   */
  template <typename Kernel, typename TargetIter, typename ResultIter>
  inline static void eval(Kernel& K,
                          const typename Kernel::multipole_type& M,
                          const typename Kernel::point_type& center,
                          TargetIter t_begin, TargetIter t_end,
			  ResultIter r_begin,
                          UseVectorM2P<false>) {
    for ( ; t_begin != t_end; ++t_begin, ++r_begin)
      K.M2P(M, center, *t_begin, *r_begin);
  }

  /** M2P evaluation dispath.
   * Detects if Kernel has a vectorized M2P and uses it if available.
   */
  template <typename Kernel, typename TargetIter, typename ResultIter>
  inline static void eval(Kernel& K,
                          const typename Kernel::multipole_type& M,
                          const typename Kernel::point_type& center,
                          TargetIter t_begin, TargetIter t_end,
			  ResultIter r_begin) {
    typedef HasM2P<Kernel,
                   const typename Kernel::multipole_type&,
                   const typename Kernel::point_type&,
                   TargetIter, TargetIter, ResultIter> HasVectorM2P;

    M2P::eval(K, M, center, t_begin, t_end, r_begin,
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
              bc.target_begin(target), bc.target_end(target),
              bc.result_begin(target));
  }
};
