#pragma once
/** @file P2P.hpp
 * @brief Dispatch methods for the P2P stage
 *
 */

struct P2P
{
  struct ONE_SIDED {};
  struct TWO_SIDED {};

  template <typename Kernel, typename BoxContext, typename BOX>
  inline static void eval(Kernel& K,
			  BoxContext& bc,
			  const BOX& source,
			  const BOX& target,
			  ONE_SIDED)
  {
#ifdef DEBUG
    printf("P2P: %d to %d\n", source.index(), target.index());
#endif

    Direct::matvec(K,
		   bc.point_begin(source), bc.point_end(source),
		   bc.charge_begin(source),
		   bc.point_begin(target), bc.point_end(target),
		   bc.result_begin(target));
  }

  // TODO: TWO_SIDED
};
