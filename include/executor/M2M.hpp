#pragma once
/** @file M2M.hpp
 * @brief Dispatch methods for the M2M stage
 *
 */

struct M2M
{
  template <typename Kernel, typename BoxContext, typename Box>
  inline static void eval(Kernel& K,
			  BoxContext& bc,
			  const Box& source,
			  const Box& target)
  {
#ifdef DEBUG
    printf("M2M: %d to %d\n", source.index(), target.index());
#endif

    K.M2M(bc.multipole_expansion(source),
          bc.multipole_expansion(target),
          bc.center(target) - bc.center(source));
  }
};
