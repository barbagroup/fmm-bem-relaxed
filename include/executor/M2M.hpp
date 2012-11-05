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

  template <typename Kernel, typename BoxContext, typename Box>
  inline static void eval(Kernel& K,
			  BoxContext& bc,
			  const Box& box)
  {
    // For all the children, M2M
    auto c_end = box.child_end();
    for (auto cit = box.child_begin(); cit != c_end; ++cit)
      M2M::eval(K, bc, *cit, box);
  }
};
