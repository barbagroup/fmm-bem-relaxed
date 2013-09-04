#pragma once

#include "EvaluatorBase.hpp"

#include "P2M.hpp"
#include "M2M.hpp"
#include "M2L.hpp"
#include "M2P.hpp"
#include "P2P.hpp"
#include "L2P.hpp"
#include "L2L.hpp"

#include "timing.hpp"

#include <functional>
#include <set>
#include <unordered_set>
#include <list>


template <typename Context, bool IS_FMM>
class EvalInteractionLazy : public EvaluatorBase<Context>
{
  //! Type of box
  typedef typename Context::box_type box_type;
  //! Pair of boxees
  typedef std::pair<box_type, box_type> box_pair;
  typedef std::pair<int, int> int_pair;
  //! List for P2P interactions    TODO: could further compress these...
  mutable std::vector<std::vector<int>> P2P_lists;
  //! List for P2M calls
  mutable std::vector<int> P2M_list;
  //! List for M2M calls
  mutable std::vector<int_pair> M2M_list;
  //! List for Long-range (M2P / M2L) interactions
  mutable std::vector<int_pair> LR_list;
  //! List for L2L calls
  mutable std::vector<int_pair> L2L_list;
  //! List for L2P calls
  mutable std::vector<int> L2P_list;
  //! Set of unsigned integers
  typedef std::unordered_set<int> Set;
  //! Set for Needed local expansions
  mutable std::vector<box_type> L_list;
  //! Set of initialised multipole expansions
  mutable Set initialised_M;
  //! Set of initialised local expansions
  mutable Set initialised_L;
  //! Set of already added L2P operations
  mutable std::set<unsigned> L2P_set;
  //! keep track of # of sparse matrix entries needed
  mutable std::vector<unsigned> mat_entries;

 public:

	/** Constructor
	 * Precompute the interaction lists, P2P_list and LR_list
	 */
	EvalInteractionLazy(Context& bc) : mat_entries(bc.source_tree().bodies()) {
    // Queue based tree traversal for P2P, M2P, and/or M2L operations
    // initialise P2P lists
    auto num_leaves = bc.source_tree().boxes();
    P2P_lists.resize(num_leaves);

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
		      // P2P_list.push_back(std::make_pair(b2.index(),b1.index()));
          // printf("P2P pushing: %d x %d\n",(int)b2.index(), (int)b1.index());
		      P2P_lists[b1.index()].push_back(b2.index());
	      } else {
		      // Split the second box into children and interact
		      auto c_end = b2.child_end();
		      for (auto cit = b2.child_begin(); cit != c_end; ++cit)
			      interact(bc, b1, *cit, pairQ);
	      }
      } else if (b2.is_leaf()) {
	      // Split the first box into children and interact
	      auto c_end = b1.child_end();
	      for (auto cit = b1.child_begin(); cit != c_end; ++cit)
		      interact(bc, *cit, b2, pairQ);
      } else {
	      // Split the larger of the two into children and interact
	      if (b1.side_length() > b2.side_length()) {
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
    // run through interaction lists and generate all call lists
    resolve_LR_interactions(bc);

    /* print out the # of P2P interactions & estimated sparse matrix size */
    /*
      int mat_elements = std::accumulate(mat_entries.begin(),mat_entries.end(),0);
      double elements_row = (double)mat_elements / mat_entries.size();
      int mat_size = mat_elements*(sizeof(double)+sizeof(int)) + (mat_entries.size()+1)*sizeof(int);
      double mat_storage = (double)mat_size / 1024. / 1024.;
      printf("Sparse matrix entries: %d (%.3lg): %.4eMB\n",mat_elements,elements_row,mat_storage);
    */
	}

	/** Execute this evaluator by applying the operators to the interaction lists
   *  Note this is implicitly cached as lists generated in the constructor
	 */
  void execute(Context& bc) const {
    // Reset/Initialise all multipole & local expansions
    auto it_end = bc.source_tree().box_end();
    for (auto it = bc.source_tree().box_begin(); it != it_end; ++it)
      INITM::eval(bc.kernel(), bc, *it);
    if (IS_FMM) {
      auto it_end = bc.target_tree().box_end();
      for (auto it = bc.target_tree().box_begin(); it != it_end; ++it)
        INITL::eval(bc.kernel(), bc, *it);
    }
    // Generate all Multipole coefficients
    eval_P2M_list(bc);
    // Evaluate all M2M operations
    eval_M2M_list(bc);
    // Evaluate queued long-range interactions
    double tic, toc, m2l_time = 0., p2p_time = 0.;
    tic = get_time();
    eval_LR_list(bc);
    toc = get_time();
    m2l_time = toc-tic;
    // Evaluate L2L operations
    eval_L2L_list(bc);
    // Evaluate L2P operations
    eval_L2P_list(bc);
    // Evaluate queued P2P interactions
    tic = get_time();
    eval_P2P_lists(bc);
    toc = get_time();
    p2p_time = toc-tic;

    printf("P2P: %.4gs, M2L (%d): %.4gs\n",p2p_time,(int)LR_list.size(),m2l_time);
  }

 private:

  /** Recursively resolve all needed multipole expansions */
  void resolve_multipole(Context& bc, const box_type& b) const
  {
    // Early exit if already initialised
    if (initialised_M.count(b.index())) return;

    if (b.is_leaf()) {
      // eval P2M
      // P2M::eval(bc.kernel(), bc, b);
      P2M_list.push_back(b.index());
    }
    else {
      // recursively call resolve_multipole on children
      for (auto it=b.child_begin(); it!=b.child_end(); ++it) {
        // resolve the lower multipole
        resolve_multipole(bc, *it);
        // now invoke M2M to get child multipoles
        M2M_list.push_back(std::make_pair(it->index(),b.index()));
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
        L2P_list.push_back(b.index());
        L2P_set.insert(b.index());
      }
    } else {
      // loop over children and propagate
      for (auto cit=b.child_begin(); cit!=b.child_end(); ++cit) {

        if (!initialised_L.count(cit->index())) {
          initialised_L.insert(cit->index());
          // call L2L on parent -> child
          L2L_list.push_back(std::make_pair(b.index(),cit->index()));
          // now recurse down the tree
          propagate_local(bc, *cit);
        }
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
      resolve_multipole(bc, bc.source_tree().box(it->first));

      // if using an FMM, mark this Local expansion as initialised and propagate down the tree
      if (IS_FMM && !initialised_L.count(it->second)) {
        initialised_L.insert(it->second);
        propagate_local(bc, bc.target_tree().box(it->second));
      }
    }
  }

  template <typename Q>
  void interact(Context& bc,
                const box_type& b1, const box_type& b2,
                Q& pairQ) const {
    if (bc.accept_multipole(b1, b2)) {
      // These boxes satisfy the multipole acceptance criteria
      // LR_list.push_back(box_pair(b1,b2));
      LR_list.push_back(std::make_pair(b1.index(),b2.index()));
      if (IS_FMM)
        L_list.push_back(b2);
    } else {
      pairQ.push_back(box_pair(b1,b2));
    }
  }

  void eval_P2P_lists(Context& bc) const
  {
    unsigned j;
#pragma omp parallel for private(j)
    for (unsigned i=0; i<P2P_lists.size(); i++) {
      // evaluate this pair using P2P
      for (j=0; j<P2P_lists[i].size(); j++) {
        P2P::eval(bc.kernel(), bc,
                  bc.source_tree().box(P2P_lists[i][j]),
                  bc.target_tree().box(i),
                  P2P::ONE_SIDED());
      }
    }
  }

  void eval_P2M_list(Context& bc) const
  {
#pragma omp parallel for
    for (unsigned i=0; i<P2M_list.size(); i++) {
      P2M::eval(bc.kernel(), bc, bc.source_tree().box(P2M_list[i]));
    }
  }

  void eval_M2M_list(Context& bc) const
  {
    for (unsigned i=0; i<M2M_list.size(); i++) {
      M2M::eval(bc.kernel(), bc, bc.source_tree().box(M2M_list[i].first), bc.source_tree().box(M2M_list[i].second));
    }
  }

  void eval_LR_list(Context& bc) const
  {
#pragma omp parallel for
    for (unsigned i=0; i<LR_list.size(); i++) {
      if (IS_FMM) {
        M2L::eval(bc.kernel(), bc,
                  bc.source_tree().box(LR_list[i].first),
                  bc.target_tree().box(LR_list[i].second));
      } else {
        M2P::eval(bc.kernel(), bc,
                  bc.source_tree().box(LR_list[i].first),
                  bc.target_tree().box(LR_list[i].second));
      }
    }
  }

  void eval_L2L_list(Context& bc) const
  {
    for (unsigned i=0; i<L2L_list.size(); i++) {
      L2L::eval(bc.kernel(), bc,
                bc.target_tree().box(L2L_list[i].first),
                bc.target_tree().box(L2L_list[i].second));
    }
  }

  void eval_L2P_list(Context& bc) const
  {
#pragma omp parallel for
    for (unsigned i=0; i<L2P_list.size(); i++) {
      L2P::eval(bc.kernel(), bc, bc.target_tree().box(L2P_list[i]));
    }
  }
};


template <typename Context, typename Options>
EvaluatorBase<Context>* make_lazy_eval(Context& c, Options& opts) {
  if (opts.evaluator == FMMOptions::FMM) {
	  return new EvalInteractionLazy<Context, true>(c);
  } else if (opts.evaluator == FMMOptions::TREECODE) {
	  return new EvalInteractionLazy<Context, false>(c);
  }
  return nullptr;
}
