#pragma once
/** @file P2M.hpp
 * @brief Dispatch methods for the P2M stage
 *
 */

#include "KernelTraits.hpp"
#include <type_traits>

class P2M
{
  /** If no other P2M dispatcher matches */
  template <typename... Args>
  inline static void eval(Args...) {
    std::cerr << "Kernel does not have a correct P2M!\n";
    exit(1);
  }

  /** P2M evaluation.
   * The Kernel provides a vector P2M accumulator.
   */
  template <typename Kernel, typename SourceIter, typename ChargeIter>
  inline static
  typename std::enable_if<ExpansionTraits<Kernel>::has_vector_P2M>::type
  eval(const Kernel& K,
       SourceIter s_begin, SourceIter s_end,
       ChargeIter c_begin,
       const typename Kernel::point_type& center,
       typename Kernel::multipole_type& M) {
    K.P2M(s_begin, s_end, c_begin, center, M);
  }

  /** P2M evaluation.
   * The Kernel provides a scalar P2M accumulator. Use it for each source point.
   */
  template <typename Kernel, typename SourceIter, typename ChargeIter>
  inline static
  typename std::enable_if<ExpansionTraits<Kernel>::has_P2M &
                          !ExpansionTraits<Kernel>::has_vector_P2M>::type
  eval(const Kernel& K,
       SourceIter s_begin, SourceIter s_end,
       ChargeIter c_begin,
       const typename Kernel::point_type& center,
       typename Kernel::multipole_type& M) {
    for ( ; s_begin != s_end; ++s_begin, ++c_begin)
      K.P2M(*s_begin, *c_begin, center, M);
  }

 public:

  /** Unwrap the data from the BoxContext and dispatch to P2M evaluator
   */
  template <typename Kernel, typename Context>
  inline static void eval(const Kernel& K,
                          Context& bc,
                          const typename Context::box_type& box)
  {
#ifdef DEBUG
    printf("P2M: %d\n", box.index());
    printf("  Bodies %d-%d\n", box.body_begin()->index(), (--box.body_end())->index());
    printf("  Bodies %d-%d\n", box.body_begin()->number(), (--box.body_end())->number());
#endif

    P2M::eval(K,
              bc.source_begin(box), bc.source_end(box), bc.charge_begin(box),
              bc.center(box), bc.multipole_expansion(box));
  }
};
