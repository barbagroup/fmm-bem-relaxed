#pragma once
/** @file M2P.hpp
 * @brief Dispatch methods for the M2P stage
 *
 */

struct M2P
{
  template <typename Kernel, typename BoxContext, typename Box>
  inline static void eval(Kernel& K,
			  BoxContext& bc,
			  const Box& source,
			  const Box& target)
  {
#ifdef DEBUG
    printf("M2P: %d to %d\n", source.index(), target.index());
#endif

    K.M2P(bc.multipole_expansion(source),
          source.center(),
          bc.point_begin(target), bc.point_end(target),
          bc.result_begin(target));
  }
};
