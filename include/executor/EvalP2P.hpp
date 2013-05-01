#pragma once

#define BOOST_UBLAS_NDEBUG
#include <boost/numeric/ublas/matrix_sparse.hpp>
namespace ublas = boost::numeric::ublas;

#include <cmath>

#include "P2P.hpp"

/** A lazy P2P evaluator which saves a list of pairs of boxes
 * That are sent to the P2P dispatcher on demand.
 */
template <typename Context>
class P2P_Lazy
    : public EvaluatorBase<Context> {
  // The context of this evaluator
  Context& bc;

  typedef typename Context::kernel_type kernel_type;
  //! Kernel value type
  typedef typename Context::kernel_value_type kernel_value_type;

  //! Type of box
  typedef typename Context::box_type box_type;
  //! Box list for P2P interactions    TODO: could further compress these...
  typedef std::pair<box_type,box_type> box_pair;
  std::vector<box_pair> p2p_list;

 public:
  P2P_Lazy(Context& _bc)
      : bc(_bc) {
  }

  /** Insert a source-target box interaction to the interaction list */
  void insert(const box_type& box1, const box_type& box2) {
    p2p_list.push_back(std::make_pair(box1,box2));
  }

  /** Compute all interations in the interaction list */
  void execute(Context& bc) const {
    for (const box_pair& b2b : p2p_list)
      P2P::eval(bc.kernel(), bc, b2b.first, b2b.second, P2P::ONE_SIDED());
  }

  /** Convert the interaction list to an interaction matrix */
  ublas::compressed_matrix<kernel_value_type> to_matrix() {
    auto first_source = bc.source_begin(bc.source_tree().root());
    auto first_target = bc.target_begin(bc.target_tree().root());

    // Interaction list for each target body
    unsigned rows = bc.target_tree().bodies();
    unsigned cols = 0;
    unsigned nnz = 0;
    std::vector<std::vector<unsigned>> csr(rows);

    for (const box_pair& b2b : p2p_list) {
      const box_type& box1 = b2b.first;
      const box_type& box2 = b2b.second;

      auto source_end = bc.source_end(box1);
      auto target_end = bc.target_end(box2);
      for (auto t = bc.target_begin(box2); t != target_end; ++t) {
        // Row
        unsigned i = t - first_target;
        std::vector<unsigned>& csri = csr[i];

        for (auto s = bc.source_begin(box1); s != source_end; ++s) {
          // Column
          unsigned j = s - first_source;

          //assert(std::find(csri.begin(), csri.end(), j) == csri.end());
          csri.push_back(j);
          ++nnz;
          cols = std::max(cols, j);
        }
      }
    }
    ++cols;

    // The precomputed interaction matrix
    ublas::compressed_matrix<kernel_value_type> m(rows, cols, nnz);

    typedef typename kernel_type::source_type source_type;
    typedef typename kernel_type::target_type target_type;
    for (unsigned i = 0; i < csr.size(); ++i) {
      // Insert elements to compressed_matrix in-order for efficiency
      std::sort(csr[i].begin(), csr[i].end());
      const target_type& target = first_target[i];

      for (unsigned j : csr[i]) {
        const source_type& source = first_source[j];
        m.push_back(i, j, bc.kernel()(target, source));
      }
    }

    return m;
  }
};
