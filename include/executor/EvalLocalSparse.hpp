#pragma once

#include "EvaluatorBase.hpp"
#include "P2P.hpp"
#include "SparseMatrix.hpp"
#include <deque>


/**
 * Only evaluate local (direct) portion of the tree
 */
template <typename Context>
class EvalLocalSparse : public EvaluatorBase<Context>
{
  //! type of box
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

  //! temp vector for charges
  mutable std::vector<charge_type> charges;
  //! temp vector for results
  mutable std::vector<result_type> results;
  //! matrix pair
  typedef std::pair<int, kernel_value_type> matrix_pair;
  //! sparse matrix
  SparseMatrix<int,kernel_value_type> A;


 public:
  // constructor -- create matrix
  EvalLocalSparse(Context& bc) {
    // source / target details
    auto num_sources = bc.source_tree().bodies();
    auto num_targets = bc.target_tree().bodies();
    // resize charge vector
    charges.resize(num_sources);
    results.resize(num_targets);
    // initialise temporary storage for the matrix
    std::vector<std::vector<matrix_pair>> temp_matrix;

    temp_matrix.resize(num_targets);
    //for (auto i=0u; i<num_targets; i++) {
      //printf("init temp_matrix[%d]: %d\n",(int)i,(int)num_sources);
      //temp_matrix[i].resize(num_sources);
    //}

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
          // modify this to generate this line of the sparse matrix
          evalP2P(bc, b1, b2, temp_matrix);
		      // P2P::eval(bc.kernel(), bc, b1, b2, P2P::ONE_SIDED());
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
    // now we have all values in temporary storage, get the total nnz
    auto nnz = 0u;
    for (auto i=0u; i<temp_matrix.size(); i++) {
      // printf("nnz += %d\n",(int)temp_matrix[i].size());
      nnz += temp_matrix[i].size();
    };
    // build matrix
    printf("Initialising matrix: %d x %d : %d\n",num_targets, num_sources, nnz);
    //A = SparseMatrix<int,kernel_value_type>(num_targets, num_sources, nnz);
    A.resize(num_targets, num_sources, nnz);
    printf("Matrix created: %.4eMB\n",(double)A.storage_size()/1024./1024.);

    // loop through temp values, streaming into sparse matrix
    int offset = 0;
    for (auto i=0u; i<temp_matrix.size(); i++) {

      for (auto j=0u; j<temp_matrix[i].size(); j++) {
        A.indices[offset+j] = temp_matrix[i][j].first;
        A.vals[offset+j] = temp_matrix[i][j].second;
      }

      offset += temp_matrix[i].size();
      A.offsets[i+1] = offset;
      // now clear this row of the temp matrix to save excess memory
      temp_matrix[i].resize(0);
    }
    // now we should have a completed matrix representation of the near-field

  } // end constructor

  void execute(Context& bc) const {
    // printf("EvalLocalSparse::execute(Context&)\n");
    // do shit here
    auto root = bc.source_tree().root();

    //double tic, toc;
    //static double t_charge = 0., t_result = 0.;

    //tic = get_time();
    int i=0;
    for (auto it=bc.charge_begin(root); it != bc.charge_end(root); ++it, i++) {
      charges[i] = *it; // *const_cast<charge_type*>(&(*it));
    }
    //toc = get_time();
    //t_charge += toc-tic;

    // call the matvec
    auto r = matvec(A, charges);
    // copy results back into iterator
    //tic = get_time();
    i=0;
    for (auto it=bc.result_begin(root); it!= bc.result_end(root); ++it, ++i) {
      *it = r[i];
    }
    //toc = get_time();
    //t_result += toc-tic;

    //printf("charge copy: %.3es, result copy: %.3es\n",t_charge, t_result);
  }

  template <typename BOX, typename Storage>
  void evalP2P(Context& bc, BOX& b1, BOX& b2, Storage& temp_matrix)
  {
    auto& K = bc.kernel();
    // targets
    auto t_it = bc.source_begin(b1);
    for (auto it=b1.body_begin(); it!=b1.body_end(); ++it, ++t_it) {
      // sources
      auto s_it = bc.source_begin(b2);
      for (auto jit=b2.body_begin(); jit!=b2.body_end(); ++jit, ++s_it) {
        // row number
        auto target_idx = it->index();
        // column number
        auto source_idx = jit->index();

        // value
        auto v = K(*t_it, *s_it);

        // insert into list / directly into matrix
        // printf("interacting bodies %d and %d: %.4lg\n", (int)target_idx, (int)source_idx, (double)v);
        temp_matrix[target_idx].push_back(std::make_pair(source_idx, v));
      }
    }
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
