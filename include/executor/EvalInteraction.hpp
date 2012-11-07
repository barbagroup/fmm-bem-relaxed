#pragma once

#include "Evaluator.hpp"
#include "Direct.hpp"

#include "M2L.hpp"
#include "M2P.hpp"
#include "P2P.hpp"

#include <functional>

template <typename Tree, typename Kernel, FMMOptions::EvalType TYPE>
class EvalInteraction : public Evaluator<EvalInteraction<Tree,Kernel,TYPE>>
{
  const Tree& tree;
  const Kernel& K;

  std::function<bool(typename Tree::box_type,
		     typename Tree::box_type)> acceptMultipole;

 public:

  template <typename MAC>
  EvalInteraction(const Tree& t, const Kernel& k, const MAC& mac)
  : tree(t), K(k), acceptMultipole(mac) {
    // any precomputation here
  }

  template <typename BoxContext, typename BOX, typename Q>
  void interact(BoxContext& bc, const BOX& b1, const BOX& b2, Q& pairQ) const {
    if (acceptMultipole(b1, b2)) {
      // These boxes satisfy the multipole acceptance criteria
      if (TYPE == FMMOptions::FMM)
	      M2L::eval(K, bc, b1, b2);
      else if (TYPE == FMMOptions::TREECODE)
	      M2P::eval(K, bc, b1, b2);
    } else if(b1.is_leaf() && b2.is_leaf()) {
      P2P::eval(K, bc, b2, b1, P2P::ONE_SIDED());
    } else {
      pairQ.push_back(std::make_pair(b1,b2));
    }
  }

  template <typename BoxContext>
  void execute(BoxContext& bc) const {
    typedef typename Tree::box_type Box;
    typedef typename std::pair<Box, Box> box_pair;
    std::deque<box_pair> pairQ;

    if(tree.root().is_leaf())
      return P2P::eval(K, bc, tree.root(), tree.root(), P2P::ONE_SIDED());

    // Queue based tree traversal for P2P, M2P, and/or M2L operations
    pairQ.push_back(box_pair(tree.root(), tree.root()));

    while (!pairQ.empty()) {
      auto b1 = pairQ.front().first;
      auto b2 = pairQ.front().second;
      pairQ.pop_front();

      if (b2.is_leaf() || (!b1.is_leaf() && b1.side_length() > b2.side_length())) {
        // Split the first box into children and interact
        auto c_end = b1.child_end();
        for (auto cit = b1.child_begin(); cit != c_end; ++cit)
          interact(bc, *cit, b2, pairQ);
      } else {
        // Split the second box into children and interact
        auto c_end = b2.child_end();
        for (auto cit = b2.child_begin(); cit != c_end; ++cit)
          interact(bc, b1, *cit, pairQ);
      }
    }
  }
};

template <typename Tree, typename Kernel, typename Options>
EvalInteraction<Tree,Kernel,FMMOptions::FMM>*
make_fmm_inter(const Tree& tree,
	       const Kernel& K,
	       const Options& opts) {
  return new EvalInteraction<Tree,Kernel,FMMOptions::FMM>(tree,K,opts.MAC);
}

template <typename Tree, typename Kernel, typename Options>
EvalInteraction<Tree,Kernel,FMMOptions::TREECODE>*
make_tree_inter(const Tree& tree,
		const Kernel& K,
		const Options& opts) {
  return new EvalInteraction<Tree,Kernel,FMMOptions::TREECODE>(tree,K,opts.MAC);
}
