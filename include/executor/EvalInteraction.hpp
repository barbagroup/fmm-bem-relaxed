#pragma once

#include "EvaluatorBase.hpp"

#include "M2L.hpp"
#include "M2P.hpp"
#include "P2P.hpp"

#include <deque>

template <typename Context, bool IS_FMM>
class EvalInteraction : public EvaluatorBase<Context>
{
  //! type of box
  typedef typename Context::box_type box_type;
  //! Pair of boxees
  typedef std::pair<box_type, box_type> box_pair;

 public:
  void execute(Context& bc) const {
    // Queue based tree traversal for P2P, M2P, and/or M2L operations
    std::deque<box_pair> pairQ;
    pairQ.push_back(box_pair(bc.source_tree().root(),
                             bc.target_tree().root()));

    while (!pairQ.empty()) {
      auto b1 = pairQ.front().first;
      auto b2 = pairQ.front().second;
      pairQ.pop_front();

      if (b1.is_leaf()) {
	      if (b2.is_leaf()) {
		      // Both are leaves, P2P
		      P2P::eval(bc.kernel(), bc, b1, b2, P2P::ONE_SIDED());
	      } else {
		      // b2 is not a leaf, Split the second box into children and interact
		      auto c_end = b2.child_end();
		      for (auto cit = b2.child_begin(); cit != c_end; ++cit)
			      interact(bc, b1, *cit, pairQ);
	      }
      } else if (b2.is_leaf()) {
	      // b1 is not a leaf, Split the first box into children and interact
	      auto c_end = b1.child_end();
	      for (auto cit = b1.child_begin(); cit != c_end; ++cit)
		      interact(bc, *cit, b2, pairQ);
      } else {
	      // Neither are leaves, Split the larger into children and interact
	      if (b1.side_length() > b2.side_length()) {   // TODO: optimize?
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
  }

  template <typename BOX, typename Q>
  void interact(Context& bc,
                const BOX& b1, const BOX& b2,
                Q& pairQ) const {
    if (bc.accept_multipole(b1, b2)) {
      // These boxes satisfy the multipole acceptance criteria
      if (IS_FMM)
	      M2L::eval(bc.kernel(), bc, b1, b2);
      else
	      M2P::eval(bc.kernel(), bc, b1, b2);
    } else {
      pairQ.push_back(box_pair(b1,b2));
    }
  }
};




template <typename Context, typename Options>
EvaluatorBase<Context>* make_interact(Context&, Options& opts) {
  if (opts.evaluator == FMMOptions::FMM)
    return new EvalInteraction<Context, true>();
  if (opts.evaluator == FMMOptions::TREECODE)
    return new EvalInteraction<Context, false>();
  return nullptr;
}
