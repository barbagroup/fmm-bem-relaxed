#pragma once

#include "EvaluatorBase.hpp"
#include "EvalP2P.hpp"
#include "Matvec.hpp"

#include <deque>

/** Only evaluate local (direct) portion of the tree
 */
template <typename Context>
class EvalLocalSparse
    : public EvaluatorBase<Context> {
  //! Type of box
  typedef typename Context::box_type box_type;
  //! Pair of boxees
  typedef std::pair<box_type, box_type> box_pair;
  //! kernel type
  typedef typename Context::kernel_type kernel_type;
  //! kernel value type
  typedef typename kernel_type::kernel_value_type kernel_value_type;
  //! kernel charge type
  typedef typename kernel_type::charge_type charge_type;
  //! kernel result type
  typedef typename kernel_type::result_type result_type;

  //! matrix pair
  typedef std::pair<int, kernel_value_type> matrix_pair;
  //! sparse matrix
  ublas::compressed_matrix<kernel_value_type> A;

 public:
  // constructor -- create matrix
  EvalLocalSparse(Context& bc) {
    // Local P2P evaluator to construct the interaction matrix
    P2P_Lazy<Context> p2p_lazy(bc);

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
          p2p_lazy.insert(b1, b2);
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

    A = p2p_lazy.to_matrix();
  } // end constructor

  void execute(Context& bc) const {
    // printf("EvalLocalSparse::execute(Context&)\n");
    // do shit here
    auto root = bc.source_tree().root();

    //double tic, toc;
    //static double t_charge = 0., t_result = 0.;

    //tic = get_time();
    typedef typename Context::charge_type charge_type;
    ublas::vector<charge_type> charges(bc.source_tree().bodies());
    std::copy(bc.charge_begin(root), bc.charge_end(root), charges.begin());
    //toc = get_time();
    //t_charge += toc-tic;

    // call the matvec
    typedef typename Context::result_type result_type;
    // ublas::vector<result_type> results = ublas::prod(A, charges);
    // ublas::vector<result_type> results = Matvec(A, charges);
    ublas::vector<result_type> results = Matvec<ublas::compressed_matrix<kernel_value_type>,
                                                ublas::vector<charge_type>,
                                                ublas::vector<result_type>>(A,charges);

    // copy results back into iterator
    //tic = get_time();
    std::transform(results.begin(), results.end(),
                   bc.result_begin(root), bc.result_begin(root),
                   std::plus<result_type>());
    //toc = get_time();
    //t_result += toc-tic;

    //printf("charge copy: %.3es, result copy: %.3es\n",t_charge, t_result);
  }

  template <typename BOX, typename Q>
  void interact(Context& bc,
                const BOX& b1, const BOX& b2,
                Q& pairQ) const {
    // ignore accepted multipoles
    if (!bc.accept_multipole(b1, b2)) {
      pairQ.push_back(box_pair(b1,b2));
    }
  }
};


template <typename Context, typename Options>
EvaluatorBase<Context>* make_sparse_local_eval(Context& bc, Options& opts) {
  (void) opts;
  return new EvalLocalSparse<Context>(bc);
}
