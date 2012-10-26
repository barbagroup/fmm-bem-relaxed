#pragma once
/** @file L2P.hpp
 * @brief Dispatch methods for the L2P stage
 *
 */

struct L2P
{
  template <typename Kernel, typename BoxContext, typename Box>
  inline static void eval(Kernel& K,
			  BoxContext& bc,
			  const Box& box)
  {
#ifdef DEBUG
    printf("L2P: %d\n", box.index());
#endif

    K.L2P(bc.local_expansion(box), box.center(),
          bc.point_begin(box), bc.point_end(box),
          bc.result_begin(box));
  }
};
