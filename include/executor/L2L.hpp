#pragma once
/** @file L2L.hpp
 * @brief Dispatch methods for the L2L stage
 *
 */

#include "KernelTraits.hpp"
#include <type_traits>

struct L2L
{
  /** If no other L2L dispatcher matches */
  template <typename... Args>
  inline static void eval(Args...) {
    std::cerr << "Kernel does not have a correct L2L!\n";
    exit(1);
  }

  template <typename Kernel>
  inline static
  typename std::enable_if<ExpansionTraits<Kernel>::has_L2L>::type
  eval(const Kernel& K,
       const typename Kernel::local_type& source,
       typename Kernel::local_type& target,
       const typename Kernel::point_type& translation) {
    K.L2L(source, target, translation);
  }

 public:

  template <typename Kernel, typename Context>
  inline static void eval(const Kernel& K,
                          Context& bc,
                          const typename Context::box_type& source,
                          const typename Context::box_type& target)
  {
#ifdef DEBUG
    printf("M2M: %d to %d\n", source.index(), target.index());
#endif

    typename Kernel::point_type r = bc.center(target) - bc.center(source);
    L2L::eval(K,
              bc.local_expansion(source),
              bc.local_expansion(target),
              r);
  }
};
