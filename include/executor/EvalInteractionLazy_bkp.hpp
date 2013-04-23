#pragma once

#include "EvaluatorBase.hpp"

#include "INITM.hpp"
#include "INITL.hpp"
#include "P2M.hpp"
#include "M2M.hpp"
#include "M2L.hpp"
#include "M2P.hpp"
#include "P2P.hpp"
#include "L2P.hpp"
#include "L2L.hpp"

#include <functional>
#include <set>
#include <unordered_set>


template <typename Context, bool IS_FMM>
class EvalInteractionLazy : public EvaluatorBase<Context>
{
  //! Type of box
  typedef typename Context::box_type box_type;
  //! Pair of boxes
  typedef std::pair<box_type, box_type> box_pair;
  //! List for P2P interactions    TODO: could further compress these...
  std::vector<box_pair> P2P_list;
  //! List for Long-range (M2P / M2L) interactions
  std::vector<box_pair> LR_list;

  //! List to recursively compute multipoles
  std::vector<box_type> multipole_cover;
  //! List to recursively initialize locals and perform downward pass
  std::vector<box_type> local_cover;


  template <typename Q>
  void interact(Context& bc,
                const box_type& b1, const box_type& b2,
                Q& pairQ) {
    if (bc.accept_multipole(b1, b2))
      // These boxes satisfy the multipole acceptance criteria
      LR_list.push_back(box_pair(b1,b2));
    else
      pairQ.push_back(box_pair(b1,b2));
  }

 public:

	/** Constructor
	 * Precompute the interaction lists, P2P_list and LR_list
	 */
	EvalInteractionLazy(Context& bc) {
    // Queue based tree traversal for P2P, M2P, and/or M2L operations
    std::deque<box_pair> pairQ;
    pairQ.push_back(box_pair(bc.source_tree().root(),
                             bc.target_tree().root()));

    while (!pairQ.empty()) {
      auto b1 = pairQ.front().first;
      auto b2 = pairQ.front().second;
      pairQ.pop_front();

      char code = (b1.is_leaf() << 1) | (b2.is_leaf() << 0);
      switch (code) {
        case 3: {    // Both boxes are leaves
          P2P_list.push_back(std::make_pair(b2,b1));
        } break;

        case 2: {    // Box 1 is a leaf
          // Split box 2 into children and interact
		      auto c_end = b2.child_end();
		      for (auto cit = b2.child_begin(); cit != c_end; ++cit)
            interact(bc, b1, *cit, pairQ);
        } break;

        case 1: {    // Box 2 is a leaf
          // Split box 1 into children and interact
          auto c_end = b1.child_end();
          for (auto cit = b1.child_begin(); cit != c_end; ++cit)
            interact(bc, *cit, b2, pairQ);
        } break;

        case 0: {    // Neither box is leaf
          // Split the larger of the two into children and interact
          if (b1.side_length() > b2.side_length()) {
            // Split box 1 into children and interact
            auto c_end = b1.child_end();
            for (auto cit = b1.child_begin(); cit != c_end; ++cit)
              interact(bc, *cit, b2, pairQ);
          } else {
            // Split box 2 into children and interact
            auto c_end = b2.child_end();
            for (auto cit = b2.child_begin(); cit != c_end; ++cit)
              interact(bc, b1, *cit, pairQ);
          }
        } break;
      }
    }

    // P2P and LR list are created.

    // Construct a covering set of multipole boxes

  }

	/** Execute this evaluator by applying the operators to the interaction lists
	 */
  void execute(Context& bc) const {
    // Evaluate queued P2P interactions
    eval_P2P_list(bc);
    // Evaluate queued long-range interactions
    eval_LR_list(bc);
  }

 private:

  /** Recursively resolve all needed multipole expansions */
  void resolve_multipole(Context& bc, const box_type& b) const
  {
    // Early exit if already initialised
    if (initialised_M.count(b.index())) return;

    // setup memory for the expansion
    INITM::eval(bc.kernel(), bc, b);

    if (b.is_leaf()) {
      // eval P2M
      // P2M::eval(bc.kernel(), bc, b);
      P2M_list.push_back(b);
    } else {
      // recursively call resolve_multipole on children
      for (auto it=b.child_begin(); it!=b.child_end(); ++it) {
        // resolve the lower multipole
        resolve_multipole(bc, *it);
        // now invoke M2M to get child multipoles
        M2M_list.push_back(std::make_pair(*it,b));
      }
    }
    // set this box as initialised
    initialised_M.insert(b.index());
  }

  /** Downward pass of L2L & L2P
   *  We can always assume that the called Local expansion has been initialised
   */
  void propagate_local(Context& bc, const box_type& b) const
  {
    if (b.is_leaf()) {
      // call L2P
      if (!L2P_set.count(b.index())) {
        L2P_list.push_back(b);
        L2P_set.insert(b.index());
      }
    } else {
      // loop over children and propagate
      for (auto cit=b.child_begin(); cit!=b.child_end(); ++cit) {
        if (!initialised_L.count(cit->index())) {
          // initialised the expansion if necessary
	        INITL::eval(bc.kernel(), bc, *cit);
          initialised_L.insert(cit->index());
        }

        // call L2L on parent -> child
        L2L_list.push_back(std::make_pair(b,*cit));
        // now recurse down the tree
        propagate_local(bc, *cit);
      }
    }
  }

  /** Evaluate long-range interactions
   *  Generate all necessary multipoles */
  // void eval_LR_list(Context& bc) const
  void resolve_LR_interactions(Context& bc) const
  {
    for (auto it=LR_list.begin(); it!=LR_list.end(); ++it) {
      // resolve all needed multipole expansions from lower levels of the tree
      resolve_multipole(bc, it->first);

      // If needed, initialise Local expansion for target box
      if (IS_FMM) {
        if (!initialised_L.count(it->second.index())) {
          // initialise Local expansion at target box if necessary
	        INITL::eval(bc.kernel(), bc, it->second);
          initialised_L.insert(it->second.index());
        }
      }
    }
  }

  void eval_L_list(Context& bc) const
  {
    for (auto it=L_list.begin(); it!=L_list.end(); ++it) {
      // propagate this local expansion down the tree
      propagate_local(bc, *it);
    }
  }

  void eval_P2P_list(Context& bc) const {
    for (box_pair b2b : P2P_list)
      P2P::eval(bc.kernel(), bc, b2b.first, b2b.second, P2P::ONE_SIDED());
  }

  void init_multipole(const box_type& b) const {
    if (M_init[b.index()]) return;

    // Initialize this multipole
    INITM::eval(bc.kernel(), bc, b);

    if (b.is_leaf()) {
      P2M::eval(bc.kernel(), bc, b);
    } else {
      // Initialize and accumulate the children
      auto c_end = b.child_end();
      for (auto ci = b.child_begin(); ci != c_end; ++ci) {
        init_multipole(*ci);
        M2M::eval(bc.kernel(), bc, *ci, b);
      }
    }

    M_init[b.index()] = true;
  }

  void eval_LR_list(Context& bc) const {
    for (auto b2b = LR_list) {
      init_multipole(b2b.first);
      if (IS_FMM)
        M2L::eval(bc.kernel(), bc, b2b.first, b2b.second);
      else
        M2P::eval(bc.kernel(), bc, b2b.first, b2b.second);
    }
  }
};


template <typename Context, typename Options>
EvaluatorBase<Context>* make_lazy_eval(Context& c, Options& opts) {
  if (opts.evaluator == FMMOptions::FMM)
	  return new EvalInteractionLazy<Context, true>(c);
  if (opts.evaluator == FMMOptions::TREECODE)
	  return new EvalInteractionLazy<Context, false>(c);
  return nullptr;
}
