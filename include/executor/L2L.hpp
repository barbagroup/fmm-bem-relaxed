#pragma once
/** @file L2L.hpp
 * @brief Dispatch methods for the L2L stage
 *
 */

struct L2L
{
  template <typename Kernel, typename BoxContext, typename Box>
  inline static void eval(Kernel& K,
			  BoxContext& bc,
			  const Box& source,
			  const Box& target)
  {
#ifdef DEBUG
    printf("L2L: %d to %d\n", source.index(), target.index());
#endif

    K.L2L(bc.local_expansion(source),
          bc.local_expansion(target),
          bc.center(target) - bc.center(source));
  }
};
