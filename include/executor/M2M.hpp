#pragma once
/** @file M2M.hpp
 * @brief Dispatch methods for the M2M stage
 *
 */

#include "KernelTraits.hpp"
#include <type_traits>

class M2M
{
  /** If no other M2M dispatcher matches */
  template <typename... Args>
  inline static void eval(Args...) {
    std::cerr << "Kernel does not have a correct M2M!\n";
    exit(1);
  }

  template <typename Kernel>
  inline static
  typename std::enable_if<ExpansionTraits<Kernel>::has_M2M>::type
  eval(const Kernel& K,
       const typename Kernel::multipole_type& source,
       typename Kernel::multipole_type& target,
       const typename Kernel::point_type& translation) {
    K.M2M(source, target, translation);
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
    M2M::eval(K,
              bc.multipole_expansion(source),
              bc.multipole_expansion(target),
              r);
  }
};
