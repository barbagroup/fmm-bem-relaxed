#pragma once

#include "Evaluator.hpp"

#include "INITM.hpp"
#include "INITL.hpp"

#include "P2M.hpp"
#include "M2M.hpp"

struct EvalUpward : public Evaluator<EvalUpward>
{
  template <typename BoxContext>
  void execute(BoxContext& bc) const {
    auto tree = bc.source_tree();
    auto K = bc.kernel();

    // For the lowest level up to the highest level
    for (unsigned l = tree.levels()-1; l != 0; --l) {
      // For all boxes at this level
      auto b_end = tree.box_end(l);
      for (auto bit = tree.box_begin(l); bit != b_end; ++bit) {
        auto box = *bit;

        // Initialize box data
        double box_size = box.side_length();
	// TODO: initialize on-demand
	INITM::eval(K, bc, box, box_size);
	INITL::eval(K, bc, box, box_size);

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
