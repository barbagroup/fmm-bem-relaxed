#pragma once

#include "Evaluator.hpp"
#include "Direct.hpp"

#include "M2L.hpp"
#include "M2P.hpp"
#include "P2P.hpp"

#include <functional>

template <typename Tree, FMMOptions::EvalType TYPE>
class EvalInteractionQueue : public Evaluator<EvalInteractionQueue<Tree,TYPE>>
{
  //! type of box
  typedef typename Tree::box_type box_type;
  //! Pair of boxees
  typedef std::pair<box_type, box_type> box_pair;
  //! List for P2P interactions
  mutable std::vector<box_pair> P2P_list;
  //! List for Long-range (M2P / M2L) interactions
  mutable std::vector<box_pair> LR_list;

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

    // now evaluate lists
    eval_LR_list(bc);
    eval_P2P_list(bc);
  }

private:

  template <typename BoxContext>
  void eval_LR_list(BoxContext& bc) const
  {
    for (auto it=LR_list.begin(); it!=LR_list.end(); ++it) {
      // evaluate this pair using M2L / M2P
      if (TYPE == FMMOptions::FMM) {
        M2L::eval(bc.kernel(), bc, it->first, it->second);
      }
      else if (TYPE == FMMOptions::TREECODE) {
        M2P::eval(bc.kernel(), bc, it->first, it->second);
      }
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
};

template <typename Tree, typename Kernel, typename Options>
EvalInteractionQueue<Tree,FMMOptions::FMM>*
make_fmm_inter_queue(const Tree&,
                     const Kernel&,
                     const Options&) {
  return new EvalInteractionQueue<Tree,FMMOptions::FMM>();
}

template <typename Tree, typename Kernel, typename Options>
EvalInteractionQueue<Tree,FMMOptions::TREECODE>*
make_tree_inter_queue(const Tree&,
                      const Kernel&,
                      const Options&) {
  return new EvalInteractionQueue<Tree,FMMOptions::TREECODE>();
}
