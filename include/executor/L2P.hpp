#pragma once
/** @file L2P.hpp
 * @brief Dispatch methods for the L2P stage
 *
 */

#include "KernelTraits.hpp"
#include <type_traits>

struct L2P
{
  /** If no other L2P dispatcher matches */
  template <typename... Args>
  inline static void eval(Args...) {
    std::cerr << "Kernel does not have a correct L2P!\n";
    exit(1);
  }

  /** L2P evaluation.
   * The Kernel provides a vector L2P accumulator.
   */
  template <typename Kernel, typename TargetIter, typename ResultIter>
  inline static
  typename std::enable_if<ExpansionTraits<Kernel>::has_vector_L2P>::type
  eval(const Kernel& K,
       const typename Kernel::local_type& L,
       const typename Kernel::point_type& center,
       TargetIter t_begin, TargetIter t_end, ResultIter r_begin) {
    K.L2P(L, center, t_begin, t_end, r_begin);
  }

  /** L2P evaluation.
   * The Kernel provides a scalar L2P accumulator. Use it for each target point.
   */
  template <typename Kernel, typename TargetIter, typename ResultIter>
  inline static
  typename std::enable_if<ExpansionTraits<Kernel>::has_L2P &
                          !ExpansionTraits<Kernel>::has_vector_L2P>::type
  eval(const Kernel& K,
       const typename Kernel::local_type& L,
       const typename Kernel::point_type& center,
       TargetIter t_begin, TargetIter t_end, ResultIter r_begin) {
    for ( ; t_begin != t_end; ++t_begin, ++r_begin)
      K.L2P(L, center, *t_begin, *r_begin);
  }

 public:

  /** Unwrap the data from BoxContext and dispatch to the L2P evaluator
   */
  template <typename Kernel, typename Context>
  inline static void eval(const Kernel& K,
                          Context& bc,
                          const typename Context::box_type& box)
  {
#ifdef DEBUG
    printf("L2P: %d\n", box.index());
#endif

    L2P::eval(K,
              bc.local_expansion(box), bc.center(box),
              bc.target_begin(box), bc.target_end(box), bc.result_begin(box));
  }
};
