#pragma once
/** @file M2P.hpp
 * @brief Dispatch methods for the M2P stage
 *
 */

#include "KernelTraits.hpp"

class M2P
{
  /** If no other M2P dispatcher matches */
  template <typename... Args>
  inline static void eval(Args...) {
    std::cerr << "Kernel does not have a correct M2P!\n";
    exit(1);
  }

  /** M2P evaluation.
   * The Kernel provides a vector M2P accumulator.
   */
  template <typename Kernel, typename TargetIter, typename ResultIter>
  inline static
  typename std::enable_if<ExpansionTraits<Kernel>::has_vector_M2P>::type
  eval(const Kernel& K,
       const typename Kernel::multipole_type& M,
       const typename Kernel::point_type& center,
       TargetIter t_begin, TargetIter t_end,
       ResultIter r_begin) {
    K.M2P(M, center, t_begin, t_end, r_begin);
  }

  /** M2P evaluation.
   * The Kernel provides a scalar M2P accumulator. Use it for each target point.
   */
  template <typename Kernel, typename TargetIter, typename ResultIter>
  inline static
  typename std::enable_if<ExpansionTraits<Kernel>::has_M2P &
                          !ExpansionTraits<Kernel>::has_vector_M2P>::type
  eval(const Kernel& K,
       const typename Kernel::multipole_type& M,
       const typename Kernel::point_type& center,
       TargetIter t_begin, TargetIter t_end,
       ResultIter r_begin) {
    for ( ; t_begin != t_end; ++t_begin, ++r_begin)
      K.M2P(M, center, *t_begin, *r_begin);
  }

 public:

  /** Unwrap the data from BoxContext and dispatch to the M2P evaluator
   */
  template <typename Kernel, typename Context>
  inline static void eval(const Kernel& K,
                          Context& bc,
                          const typename Context::box_type& source,
                          const typename Context::box_type& target) {
#ifdef DEBUG
    printf("M2P: %d to %d\n", source.index(), target.index());
#endif

    M2P::eval(K,
              bc.multipole_expansion(source), bc.center(source),
              bc.target_begin(target), bc.target_end(target),
              bc.result_begin(target));
  }
};
