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
          target.center() - source.center());
  }

  template <typename Kernel, typename BoxContext, typename Box>
  inline static void eval(Kernel& K,
			  BoxContext& bc,
			  const Box& box)
  {
    // For all the children, L2L
    auto c_end = box.child_end();
    for (auto cit = box.child_begin(); cit != c_end; ++cit)
      L2L::eval(K, bc, box, *cit);
  }
};
