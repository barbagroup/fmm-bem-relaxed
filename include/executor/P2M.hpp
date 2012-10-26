#pragma once
/** @file P2M.hpp
 * @brief Dispatch methods for the P2M stage
 *
 */

struct P2M
{
  template <typename Kernel, typename BoxContext, typename Box>
  inline static void eval(Kernel& K,
			  BoxContext& bc,
			  const Box& box)
  {
#ifdef DEBUG
    printf("P2M: %d\n", box.index());
#endif

    K.P2M(bc.point_begin(box), bc.point_end(box), bc.charge_begin(box),
          box.center(),
          bc.multipole_expansion(box));
  }
};
