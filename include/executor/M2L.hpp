#pragma once
/** @file M2L.hpp
 * @brief Dispatch methods for the M2L stage
 *
 */

struct M2L
{
  template <typename Kernel, typename BoxContext, typename Box>
  inline static void eval(Kernel& K,
			  BoxContext& bc,
			  const Box& source,
			  const Box& target)
  {
#ifdef DEBUG
    printf("M2L: %d to %d\n", source.index(), target.index());
#endif

    K.M2L(bc.multipole_expansion(source),
          bc.local_expansion(target),
          bc.center(target) - bc.center(source));
  }
};
