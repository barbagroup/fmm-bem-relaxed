#pragma once
/** @file P2P.hpp
 * @brief Dispatch methods for the P2P stage
 *
 */

struct P2P
{
  struct ONE_SIDED {};
  struct TWO_SIDED {};

  /** One sided P2P
   */
  template <typename Kernel, typename BoxContext, typename Box>
  inline static void eval(const Kernel& K,
			  BoxContext& bc,
			  const Box& source,
			  const Box& target,
			  ONE_SIDED)
  {
#ifdef DEBUG
    printf("P2P: %d to %d\n", source.index(), target.index());
#endif

    Direct::matvec(K,
		   bc.source_begin(source), bc.source_end(source),
		   bc.charge_begin(source),
		   bc.target_begin(target), bc.target_end(target),
		   bc.result_begin(target));
  }

  /** Two sided P2P
   */
  template <typename Kernel, typename BoxContext, typename Box>
  inline static void eval(const Kernel& K,
			  BoxContext& bc,
			  const Box& source,
			  const Box& target,
			  TWO_SIDED)
  {
#ifdef DEBUG
    printf("P2P: %d to %d\n", source.index(), target.index());
#endif

    Direct::matvec(K,
		   bc.source_begin(source), bc.source_end(source),
                   bc.charge_begin(source),
		   bc.target_begin(target), bc.target_end(target),
                   bc.charge_begin(target),
                   bc.result_begin(source), bc.result_begin(target));
  }
};
