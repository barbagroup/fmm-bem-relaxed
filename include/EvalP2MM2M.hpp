#pragma once

#include "Evaluator.hpp"


// FMM / treecode upward sweep
template <typename Tree, typename Kernel, typename Options>
class EvalP2MM2M : public Evaluator<EvalP2MM2M<Tree,Kernel,Options>>
{
  const Tree& tree;
  const Kernel& K;

  template <typename BoxContext, typename BOX>
      void evalP2M(BoxContext& bc,
                   const BOX& box) const {
    // Point iterator pair
    auto p_begin = bc.point_begin(box);
    auto p_end   = bc.point_end(box);
    // Charge iterator
    auto c_begin = bc.charge_begin(box);

#ifdef DEBUG
    printf("P2M: box: %d\n", box.index());
#endif
    K.P2M(p_begin, p_end, c_begin,
          box.center(),
          bc.multipole_expansion(box));
  }

  template <typename BoxContext, typename BOX>
      void evalM2M(BoxContext& bc,
                   const BOX& cbox,
                   const BOX& box) const {
#ifdef DEBUG
    printf("M2M: %d to %d\n", cbox.index(), box.index());
#endif
    K.M2M(bc.multipole_expansion(cbox),
          bc.multipole_expansion(box),
          box.center() - cbox.center());
  }

public:
  //! constructor
  EvalP2MM2M(const Tree& t, const Kernel& k, const Options& options)
      : tree(t), K(k) {
    (void) options;
  };

  template <typename BoxContext>
      void execute(BoxContext& bc) const {
#ifdef DEBUG
    unsigned lowest_level = tree.levels();
    printf("lowest level in tree: %d\n",(int)lowest_level);
#endif

    // For the lowest level up to the highest level
    for (unsigned l = tree.levels()-1; l != 0; --l) {
      // For all boxes at this level
      auto b_end = tree.box_end(l);
      for (auto bit = tree.box_begin(l); bit != b_end; ++bit) {
        auto box = *bit;

        // Initialize box data
        double box_size = box.side_length();
        // K.init_multipole(M[idx], box_size);
        K.init_multipole(bc.multipole_expansion(box), box_size);
        // K.init_local(L[idx], box_size);
        K.init_local(bc.local_expansion(box), box_size);

        if (box.is_leaf()) {
          // If leaf, make P2M calls
          evalP2M(bc, box);

        } else {
          // If not leaf, make M2M calls

          // For all the children, M2M
          auto c_end = box.child_end();
          for (auto cit = box.child_begin(); cit != c_end; ++cit) {
            auto cbox = *cit;
            evalM2M(bc, cbox, box);
          }
        }
      }
    }
  }
};

template <typename Tree, typename Kernel, typename Options>
EvalP2MM2M<Tree,Kernel,Options>* make_P2MM2M(const Tree& tree,
					     const Kernel& kernel,
					     const Options& opts) {
  return new EvalP2MM2M<Tree,Kernel,Options>(tree, kernel, opts);
}

