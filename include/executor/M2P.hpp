#pragma once
/** @file M2P.hpp
 * @brief Dispatch methods for the M2P stage
 *
 */

#include "KernelTraits.hpp"


struct M2P
{
  template <bool> struct UseVectorM2P {};
  template <bool> struct UseM2P {};

  /** The Kernel does not provide an M2P
   */
  template <typename Kernel, typename TargetIter, typename ResultIter>
  inline static void eval(Kernel&,
                          const typename Kernel::multipole_type&,
                          const typename Kernel::point_type&,
                          TargetIter, TargetIter,
			  ResultIter,
			  UseM2P<false>,
                          UseVectorM2P<false>) {
    std::cerr << "Error! Kernel has no M2P method to call.\n";
    exit(1);
  }

  /** M2P evaluation.
   * The Kernel provides a vector M2P accumulator.
   */
  template <typename Kernel, typename TargetIter, typename ResultIter,
	    bool A>
  inline static void eval(Kernel& K,
                          const typename Kernel::multipole_type& M,
                          const typename Kernel::point_type& center,
                          TargetIter t_begin, TargetIter t_end,
			  ResultIter r_begin,
			  UseM2P<A>,
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
			  UseM2P<true>,
                          UseVectorM2P<false>) {
    for ( ; t_begin != t_end; ++t_begin, ++r_begin)
      K.M2P(M, center, *t_begin, *r_begin);
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
              bc.result_begin(target),
	      UseM2P<ExpansionTraits<Kernel>::has_M2P>(),
              UseVectorM2P<ExpansionTraits<Kernel>::has_vector_M2P>());
  }
};
