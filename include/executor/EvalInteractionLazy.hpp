#pragma once

#include "Evaluator.hpp"
#include "Direct.hpp"

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

template <typename Tree, FMMOptions::EvalType TYPE>
class EvalInteractionLazy : public Evaluator<EvalInteractionLazy<Tree,TYPE>>
{
  //! type of box
  typedef typename Tree::box_type box_type;
  //! Pair of boxees
  typedef std::pair<box_type, box_type> box_pair;
  //! List for P2P interactions
  mutable std::vector<box_pair> P2P_list;
  //! List for Long-range (M2P / M2L) interactions
  mutable std::vector<box_pair> LR_list;
  //! Set of unsigned integers
  typedef std::unordered_set<int> Set;
  //! Set for Needed local expansions
  mutable std::set<box_type> L_list;
  //! Set of initialised multipole expansions
  mutable Set initialised_M;
  //! Set of initialised local expansions
  mutable Set initialised_L;

 public:

  template <typename Kernel, typename STree, typename TTree, typename Options>
  void create(const Kernel&, STree&, TTree&, Options&) {
    //acceptMultipole = opts.MAC();
  }

  template <typename BoxContext>
  void execute(BoxContext& bc) const {
    typedef typename Tree::box_type Box;
    typedef typename std::pair<Box, Box> box_pair;
    std::deque<box_pair> pairQ;

    auto stree = bc.source_tree();
    auto ttree = bc.target_tree();

    if(stree.root().is_leaf() || ttree.root().is_leaf())
      return P2P::eval(bc.kernel(), bc,
                       stree.root(), ttree.root(), P2P::ONE_SIDED());

    // Queue based tree traversal for P2P, M2P, and/or M2L operations
    pairQ.push_back(box_pair(stree.root(), ttree.root()));

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

    // Evaluate queued long-range interactions
    eval_LR_list(bc);
    // Evaluate queued P2P interactions
    eval_P2P_list(bc);
    // Evaluate queued L2L / L2P interactions
    eval_L_list(bc);
  }

private:

  /** Recursively resolve all needed multipole expansions */
  template <typename BoxContext, typename BOX>
  void resolve_multipole(BoxContext& bc, BOX b) const
  {
    // Early exit if already initialised
    if (initialised_M.count(b.index())) return;

    // setup memory for the expansion
    bc.kernel().init_multipole(bc.multipole_expansion(b),b.side_length());

    if (b.is_leaf()) {
      // eval P2M
      P2M::eval(bc.kernel(), bc, b);
    }
    else {
      // recursively call resolve_multipole on children
      for (auto it=b.child_begin(); it!=b.child_end(); ++it) {
        // resolve the lower multipole
        resolve_multipole(bc, *it);
        // now invoke M2M to get child multipoles
        M2M::eval(bc.kernel(), bc, *it, b);
      }
    }
    // set this box as initialised
    initialised_M.insert(b.index());
  }

  /** Downward pass of L2L & L2P
   *  We can always assume that the called Local expansion has been initialised */
  template <typename BoxContext>
  void propagate_local(BoxContext& bc, box_type b) const
  {
    if (b.is_leaf()) {
      // call L2P
      L2P::eval(bc.kernel(), bc, b);
    } else {
      // loop over children and propagate
      for (auto cit=b.child_begin(); cit!=b.child_end(); ++cit) {
        if (!initialised_L.count(cit->index())) {
          // initialised the expansion if necessary
          bc.kernel().init_local(bc.local_expansion(*cit),0);
          initialised_L.insert(cit->index());
        }

        // call L2L on parent -> child
        L2L::eval(bc.kernel(), bc, b, *cit);
        // now recurse down the tree
        propagate_local(bc, *cit);
      }
    }
  }

  template <typename BoxContext>
  void eval_LR_list(BoxContext& bc) const
  {
    for (auto it=LR_list.begin(); it!=LR_list.end(); ++it) {
      // resolve all needed multipole expansions from lower levels of the tree
      resolve_multipole(bc, it->first);

      // evaluate this pair using M2L / M2P
      if (TYPE == FMMOptions::FMM) {
        if (!initialised_L.count(it->second.index())) {
          // initialise Local expansion at target box if necessary
          bc.kernel().init_local(bc.local_expansion(it->second),0);
          initialised_L.insert(it->second.index());
        }

        // Perform the translation
        M2L::eval(bc.kernel(), bc, it->first, it->second);
      }
      else if (TYPE == FMMOptions::TREECODE) {
        M2P::eval(bc.kernel(), bc, it->first, it->second);
      }
    }
  }

  template <typename BoxContext>
  void eval_L_list(BoxContext& bc) const
  {
    for (auto it=L_list.begin(); it!=L_list.end(); ++it) {
      // propagate this local expansion down the tree
      propagate_local(bc, *it);
    }
  }

  template <typename BoxContext>
  void eval_P2P_list(BoxContext& bc) const
  {
    for (auto it=P2P_list.begin(); it!=P2P_list.end(); ++it) {
      // evaluate this pair using P2P
      P2P::eval(bc.kernel(), bc, it->first, it->second, P2P::ONE_SIDED());
    }
  }

  template <typename BoxContext, typename BOX, typename Q>
  void interact(BoxContext& bc, const BOX& b1, const BOX& b2, Q& pairQ) const {
    if (bc.accept_multipole(b1, b2)) {
      // These boxes satisfy the multipole acceptance criteria
      if (TYPE == FMMOptions::FMM) {
        LR_list.push_back(std::make_pair(b1,b2));
        L_list.insert(b2); // push_back(b2);
      }
      else if (TYPE == FMMOptions::TREECODE) {
        LR_list.push_back(std::make_pair(b1,b2));
      }
    } else if(b1.is_leaf() && b2.is_leaf()) {
      P2P_list.push_back(std::make_pair(b2,b1));
    } else {
      pairQ.push_back(std::make_pair(b1,b2));
    }
  }
};
