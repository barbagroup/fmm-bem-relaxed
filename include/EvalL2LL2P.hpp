#pragma once

#include "Evaluator.hpp"



template <typename Tree, typename Kernel, typename Options>
class EvalL2LL2P : public Evaluator<EvalL2LL2P<Tree,Kernel,Options>>
{
  const Tree& tree;
  const Kernel& K;
public:

  EvalL2LL2P(const Tree& t, const Kernel& k, const Options& options)
      : tree(t), K(k) {
    // any precomputation here
    (void) options;
  }

  template <typename BoxContext, typename BOX>
      void evalL2P(BoxContext& bc,
                   const BOX& box) const {
    auto t_begin = bc.point_begin(box);
    auto t_end   = bc.point_end(box);
    auto r_begin = bc.result_begin(box);

    printf("L2P: %d\n",box.index());
    K.L2P(bc.local_expansion(box), box.center(),
          t_begin, t_end,
          r_begin);
  }

  template <typename BoxContext, typename BOX>
      void evalL2L(BoxContext& bc,
                   const BOX& box,
                   const BOX& cbox) const {
    printf("L2L: %d to %d\n", box.index(), cbox.index());
    K.L2L(bc.local_expansion(box),
          bc.local_expansion(cbox),
          cbox.center() - box.center());
  }

  template <typename BoxContext>
      void execute(BoxContext& bc) const {
    // For the highest level down to the lowest level
    for (unsigned l = 1; l < tree.levels(); ++l) {
      // For all boxes at this level
      auto b_end = tree.box_end(l);
      for (auto bit = tree.box_begin(l); bit != b_end; ++bit) {
        auto box = *bit;

        // Initialize box data
        if (box.is_leaf()) {
          // If leaf, make L2P calls
          evalL2P(bc, box);
        } else {
          // If not leaf, make L2L calls

          // For all the children, L2L
          auto c_end = box.child_end();
          for (auto cit = box.child_begin(); cit != c_end; ++cit) {
            auto cbox = *cit;
            evalL2L(bc, box, cbox);
          }
        }
      }
    }
  }
};

template <typename Tree, typename Kernel, typename Options>
EvalL2LL2P<Tree,Kernel,Options>* make_L2LL2P(const Tree& tree,
					     const Kernel& kernel,
					     const Options& opts) {
  return new EvalL2LL2P<Tree,Kernel,Options>(tree, kernel, opts);
}
