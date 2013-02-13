#pragma once
/** @file M2L.hpp
 * @brief Dispatch methods for the M2L stage
 *
 */

#include "KernelTraits.hpp"
#include <type_traits>

class M2L
{
  /** If no other M2L dispatcher matches */
  template <typename... Args>
  inline static void eval(Args...) {
    std::cerr << "Kernel does not have a correct M2L!\n";
    exit(1);
  }

  template <typename Kernel>
  inline static
  typename std::enable_if<ExpansionTraits<Kernel>::has_M2L>::type
  eval(const Kernel& K,
       const typename Kernel::multipole_type& source,
       typename Kernel::local_type& target,
       const typename Kernel::point_type& translation,
       unsigned p) {
    K.M2L(source, target, translation, p);
  }

 public:

  template <typename Kernel, typename Context>
  inline static void eval(Kernel& K,
                          Context& bc,
                          const typename Context::box_type& source,
                          const typename Context::box_type& target,
                          unsigned p)
  {
#ifdef DEBUG
    printf("M2L: %d to %d\n", source.index(), target.index());
#endif

    M2L::eval(K,
              bc.multipole_expansion(source), bc.local_expansion(target),
              bc.center(target) - bc.center(source), p);
  }
};
