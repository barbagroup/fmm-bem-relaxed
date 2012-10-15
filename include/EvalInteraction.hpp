#pragma once

#include "Evaluator.hpp"

#include <functional>

template <typename Tree, typename Kernel, int EXECCODE>
class EvalInteraction : public Evaluator<EvalInteraction<Tree,Kernel,EXECCODE>>
{
  const Tree& tree;
  const Kernel& K;

  std::function<bool(typename Tree::box_type,
		     typename Tree::box_type)> acceptMultipole;

  /** One-sided P2P!!
   */
  template <typename BoxContext, typename BOX>
  void evalP2P(BoxContext& bc,
	       const BOX& b1,
	       const BOX& b2) const {
    // Point iters
    auto p1_begin = bc.point_begin(b1);
    auto p1_end   = bc.point_end(b1);
    auto p2_begin = bc.point_begin(b2);
    auto p2_end   = bc.point_end(b2);

    // Charge iters
    auto c1_begin = bc.charge_begin(b1);

    // Result iters
    auto r2_begin = bc.result_begin(b2);

    printf("P2P: %d to %d\n",b1.index(),b2.index());
    K.P2P(p1_begin, p1_end, c1_begin,
          p2_begin, p2_end,
          r2_begin);
  }

  template <typename BoxContext, typename BOX>
  void evalM2P(BoxContext& bc,
	       const BOX& b1,
	       const BOX& b2) const {
    // Target point iters
    auto t_begin = bc.point_begin(b2);
    auto t_end   = bc.point_end(b2);

    // Target result iters
    auto r_begin = bc.result_begin(b2);

    printf("M2P: %d to %d\n", b1.index(), b2.index());

    K.M2P(bc.multipole_expansion(b1),
          b1.center(),
          t_begin, t_end,
          r_begin);
  }

  template <typename BoxContext, typename BOX>
  void evalM2L(BoxContext& bc,
	       const BOX& b1,
	       const BOX& b2) const {
    printf("M2L: %d to %d\n", b2.index(), b1.index());

    K.M2L(bc.multipole_expansion(b1),
          bc.local_expansion(b2),
          b2.center() - b1.center());
  }

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
      if (EXECCODE == 0)
	evalM2L(bc, b1, b2);
      else if (EXECCODE == 1)
	evalM2P(bc, b1, b2);
    } else if(b1.is_leaf() && b2.is_leaf()) {
      evalP2P(bc, b2, b1);
    } else {
      pairQ.push_back(std::make_pair(b1,b2));
    }
  }

  template <typename BoxContext>
  void execute(BoxContext& bc) const {
    typedef typename Tree::box_type Box;
    typedef typename std::pair<Box, Box> box_pair;
    std::deque<box_pair> pairQ;

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


/*
template <typename Tree, typename Kernel, typename Options>
EvalInteraction<Tree,Kernel,int>* make_inter(const Tree& tree,
					     const Kernel& K,
					     const Options& opts) {
  if (opts.evaluator == FMMOptions::FMM)
    return new EvalInteraction<Tree,Kernel,0>(tree,K,opts);
  else
    return new EvalInteraction<Tree,Kernel,1>(tree,K,opts);
}
*/

