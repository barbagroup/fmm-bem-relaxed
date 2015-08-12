#pragma once

#include "EvaluatorBase.hpp"
#include "EvalP2P.hpp"
#include "Matvec.hpp"

#include <deque>

/** Only evaluate local (direct) portion of the tree
 */
template <typename Context>
class EvalDiagonalSparse
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
  EvalDiagonalSparse(Context& bc) {
    // Local P2P evaluator to construct the interaction matrix
    P2P_Lazy<Context> p2p_lazy(bc);

    // get the source root
    auto& tree = bc.source_tree();

    // loop through all leaf boxes, pushing back self-interactions
    for (auto bi = tree.box_begin(); bi!=tree.box_end(); ++bi)
    {
      if (bi->is_leaf()) {
        p2p_lazy.insert(*bi,*bi);
      }
    }

    A = p2p_lazy.to_matrix();
  } // end constructor

  void execute(Context& bc) const {
    auto root = bc.source_tree().root();
    typedef typename Context::charge_type charge_type;
    ublas::vector<charge_type> charges(bc.source_tree().bodies());
    std::copy(bc.charge_begin(root), bc.charge_end(root), charges.begin());

    // call the matvec
    typedef typename Context::result_type result_type;
    // ublas::vector<result_type> results = ublas::prod(A, charges);
    // ublas::vector<result_type> results = Matvec(A, charges);
    ublas::vector<result_type> results = Matvec<ublas::compressed_matrix<kernel_value_type>,
                                                ublas::vector<charge_type>,
                                                ublas::vector<result_type>>(A,charges);

    // copy results back into iterator
    std::transform(results.begin(), results.end(),
                   bc.result_begin(root), bc.result_begin(root),
                   std::plus<result_type>());
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
EvaluatorBase<Context>* make_sparse_diagonal_eval(Context& bc, Options& opts) {
  (void) opts;
  return new EvalDiagonalSparse<Context>(bc);
}
