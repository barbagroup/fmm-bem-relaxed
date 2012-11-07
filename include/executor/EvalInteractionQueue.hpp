#pragma once

#include "Evaluator.hpp"
#include "Direct.hpp"

#include "M2L.hpp"
#include "M2P.hpp"
#include "P2P.hpp"

#include <functional>

template <typename Tree, typename Kernel, FMMOptions::EvalType TYPE>
class EvalInteractionQueue : public Evaluator<EvalInteractionQueue<Tree,Kernel,TYPE>>
{
  const Tree& tree;
  const Kernel& K;

  //! type of box
  typedef typename Tree::box_type box_type;
  //! Pair of boxees
  typedef std::pair<box_type, box_type> box_pair;
  //! List for P2P interactions
  mutable std::vector<box_pair> P2P_list;
  //! List for Long-range (M2P / M2L) interactions
  mutable std::vector<box_pair> LR_list;

  std::function<bool(typename Tree::box_type,
		     typename Tree::box_type)> acceptMultipole;

 public:

  template <typename MAC>
  EvalInteractionQueue(const Tree& t, const Kernel& k, const MAC& mac)
  : tree(t), K(k), P2P_list(0), LR_list(0), acceptMultipole(mac) {
    // any precomputation here
  }

  template <typename BoxContext>
  void eval_LR_list(BoxContext& bc) const
  {
    for (auto it=LR_list.begin(); it!=LR_list.end(); ++it) {
      // evaluate this pair using M2L / M2P
      if (TYPE == FMMOptions::FMM) {
        M2L::eval(K,bc,it->first,it->second);
      }
      else if (TYPE == FMMOptions::TREECODE) {
        M2P::eval(K,bc,it->first,it->second);
      }
    }
  }

  template <typename BoxContext>
  void eval_P2P_list(BoxContext& bc) const 
  {
    for (auto it=P2P_list.begin(); it!=P2P_list.end(); ++it) {
      // evaluate this pair using P2P
      P2P::eval(K,bc,it->first,it->second,P2P::ONE_SIDED());
    }
  }

  template <typename BoxContext, typename BOX, typename Q>
  void interact(BoxContext& bc, const BOX& b1, const BOX& b2, Q& pairQ) const {
    (void) bc; // quiet warning for now
    if (acceptMultipole(b1, b2)) {
      // These boxes satisfy the multipole acceptance criteria
      if (TYPE == FMMOptions::FMM) {
	      // M2L::eval(K, bc, b1, b2);
        LR_list.push_back(std::make_pair(b1,b2));
      }
      else if (TYPE == FMMOptions::TREECODE) {
	      // M2P::eval(K, bc, b1, b2);
        LR_list.push_back(std::make_pair(b1,b2));
      }
    } else if(b1.is_leaf() && b2.is_leaf()) {
      // P2P::eval(K, bc, b2, b1, P2P::ONE_SIDED());
      P2P_list.push_back(std::make_pair(b2,b1));
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

    // now evaluate lists
    eval_LR_list(bc);
    eval_P2P_list(bc);
  }
};

template <typename Tree, typename Kernel, typename Options>
EvalInteractionQueue<Tree,Kernel,FMMOptions::FMM>*
make_fmm_inter_queue(const Tree& tree,
	       const Kernel& K,
	       const Options& opts) {
  return new EvalInteractionQueue<Tree,Kernel,FMMOptions::FMM>(tree,K,opts.MAC);
}

template <typename Tree, typename Kernel, typename Options>
EvalInteraction<Tree,Kernel,FMMOptions::TREECODE>*
make_tree_inter_queue(const Tree& tree,
		const Kernel& K,
		const Options& opts) {
  return new EvalInteractionQueue<Tree,Kernel,FMMOptions::TREECODE>(tree,K,opts.MAC);
}
