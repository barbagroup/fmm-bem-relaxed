#pragma once

#include "Evaluator.hpp"

#include "P2M.hpp"
#include "M2M.hpp"

template <typename Tree, typename Kernel>
class EvalUpward : public Evaluator<EvalUpward<Tree,Kernel>>
{
  const Tree& tree;
  const Kernel& K;

public:
  EvalUpward(const Tree& t, const Kernel& k)
  : tree(t), K(k) {
  };

  template <typename BoxContext>
  void execute(BoxContext& bc) const {
    // For the lowest level up to the highest level
    for (unsigned l = tree.levels()-1; l != 0; --l) {
      // For all boxes at this level
      auto b_end = tree.box_end(l);
      for (auto bit = tree.box_begin(l); bit != b_end; ++bit) {
        auto box = *bit;

        // Initialize box data
        double box_size = box.side_length();
	// TODO: initialize on-demand
        K.init_multipole(bc.multipole_expansion(box), box_size);
        K.init_local(bc.local_expansion(box), box_size);

        if (box.is_leaf()) {
          // If leaf, make P2M calls
	  P2M::eval(K, bc, box);
        } else {
          // If not leaf, then for all the children M2M
	  auto c_end = box.child_end();
	  for (auto cit = box.child_begin(); cit != c_end; ++cit)
	    M2M::eval(K, bc, *cit, box);
        }
      }
    }
  }
};

template <typename Tree, typename Kernel, typename Options>
EvalUpward<Tree,Kernel>* make_upward(const Tree& tree,
				     const Kernel& kernel,
				     const Options&) {
  return new EvalUpward<Tree,Kernel>(tree, kernel);
}

