#pragma once
/** @file INITL.hpp
 * @brief Dispatch methods for the initializing a local expansion
 *
 */

#include "KernelTraits.hpp"
#include <type_traits>

struct INITL
{
  /** If no init_local dispatcher matches */
  template <typename... Args>
  inline static void eval(Args...) {
    // Do nothing
  }

  template <typename Kernel>
  inline static
  typename std::enable_if<ExpansionTraits<Kernel>::has_init_local>::type
  eval(const Kernel& K,
       typename Kernel::local_type& L,
       const typename Kernel::point_type& extents,
       unsigned level) {
    K.init_local(L, extents, level);
  }

  template <typename Kernel, typename Context>
  inline static void eval(const Kernel& K,
                          Context& bc,
                          const typename Context::box_type& b)
  {
#ifdef DEBUG
    printf("initL: %d\n", b.index());
#endif

    INITL::eval(K, bc.local_expansion(b), b.extents(), b.level());
  }
};
