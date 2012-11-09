#pragma once

#include "ExecutorBase.hpp"
#include "Evaluator.hpp"

#include "TransformIterator.hpp"

#include <type_traits>

/** @class Executor
 * @brief A very general Executor class. This provides a context to any tree
 * that provides the following interface:
 * Tree
 *   Tree(Iter, Iter, Options)          // Tree constructor
 * Tree::box_type
 *   int index() const                  // Each box has a unique index
 *   point_type center() const          // Each box has a center point
 *   body_iterator body_begin() const
 *   body_iterator body_end() const     // Iterator pair to the box's bodies
 * Tree::body_type
 *   int index() const                  // Each body has a unique index
 *   Kernel::point_type point() const   // A body has a point_type location
 * This class assumes nothing else about the tree.
 */
template <typename Kernel, typename Tree, typename Eval>
class MinimalExecutor : public ExecutorBase<Kernel>
{
  //! Tree box type
  typedef typename Tree::box_type box_type;
  //! Tree body type
  typedef typename Tree::body_type body_type;
  //! Kernel charge type
  typedef typename Kernel::multipole_type multipole_type;
  //! Kernel result type
  typedef typename Kernel::local_type local_type;
  //! Kernel point type
  typedef typename Kernel::point_type point_type;
  //! Kernel source type
  typedef typename Kernel::source_type source_type;
  //! Kernel target type
  typedef typename Kernel::target_type target_type;
  //! Kernel charge type
  typedef typename Kernel::charge_type charge_type;
  //! Kernel result type
  typedef typename Kernel::result_type result_type;

private:
  //! The tree of sources
  Tree source_tree_;
  //! Evaluator algorithm to apply
  Eval eval_;
  //! Multipole expansions corresponding to Box indices in Tree
  std::vector<multipole_type> M_;
  //! Local expansions corresponding to Box indices in the Tree
  std::vector<local_type> L_;
  //! The sources associated with bodies in the source_tree
  std::vector<source_type> s_;
  //! The charges associated with bodies in the source_tree
  std::vector<charge_type> c_;
  //! The targets associated with bodies in the target_tree
  std::vector<target_type> t_;
  //! The results associated with bodies in the target_tree
  std::vector<result_type> r_;


  template <typename T>
  struct body_map {
    typedef typename std::remove_const<T>::type vT;
    typedef typename std::vector<vT>::iterator iterator;
    typedef typename std::vector<vT>::const_iterator const_iterator;
    static constexpr bool is_const = std::is_const<T>::value;
    typename std::conditional<is_const, const_iterator, iterator>::type begin_;
    T& operator()(const body_type& b) const {
      return begin_[b.number()];
    }
  };

  typedef body_map<const source_type> body2source_t;
  typedef body_map<const charge_type> body2charge_t;
  typedef body_map<const target_type> body2target_t;
  typedef body_map<result_type>       body2result_t;

  body2source_t body2source;
  body2charge_t body2charge;
  body2target_t body2target;
  body2result_t body2result;

  typedef typename Tree::body_iterator                  body_iterator;
  typedef transform_iterator<body_iterator,body2source_t> source_iterator;
  typedef transform_iterator<body_iterator,body2charge_t> charge_iterator;
  typedef transform_iterator<body_iterator,body2target_t> target_iterator;
  typedef transform_iterator<body_iterator,body2result_t> result_iterator;

public:
  //! Constructor
  template <typename SourceIter, typename Options>
  MinimalExecutor(const Kernel& K,
		  SourceIter first, SourceIter last,
		  Options& opts)
    : source_tree_(first, last, opts),
      eval_(K, source_tree_, opts),
      M_(source_tree_.boxes()),
      L_(source_tree_.boxes()),
      s_(first, last),
      t_(first, last),  // TEMP HACK
      r_(t_.size()) {
    // TODO: Init L_ on the target_tree
  }

#if 0
  void set_sources(const std::vector<charge_type>& sources) {
  }

  void set_charges(const std::vector<charge_type>& charges) {

  }

  void set_targets(const std::vector<charge_type>& targets) {
  }

  void get_results() {
  }
#endif

  void execute(const std::vector<typename Kernel::charge_type>& charges,
	       std::vector<typename Kernel::result_type>& results) {
    // TEMP Hack
    body2source.begin_ = s_.begin();
    body2charge.begin_ = charges.begin();
    body2target.begin_ = t_.begin();
    body2result.begin_ = results.begin();
    eval_.execute(*this);
  }

  // Accessors to make this Executor into a BoxContext
  inline multipole_type& multipole_expansion(const box_type& b) {
    return M_[b.index()];
  }
  inline const multipole_type& multipole_expansion(const box_type& b) const {
    return M_[b.index()];
  }
  inline local_type& local_expansion(const box_type& b) {
    return L_[b.index()];
  }
  inline const multipole_type& local_expansion(const box_type& b) const {
    return L_[b.index()];
  }

  inline point_type center(const box_type& b) const {
    return b.center();
  }

  inline source_iterator source_begin(const box_type& b) {
    return make_transform_iterator(b.body_begin(), body2source);
  }
  inline source_iterator source_end(const box_type& b) {
    return make_transform_iterator(b.body_end(), body2source);
  }
  inline charge_iterator charge_begin(const box_type& b) {
    return make_transform_iterator(b.body_begin(), body2charge);
  }
  inline charge_iterator charge_end(const box_type& b) {
    return make_transform_iterator(b.body_end(), body2charge);
  }
  inline target_iterator target_begin(const box_type& b) {
    return make_transform_iterator(b.body_begin(), body2source);
  }
  inline target_iterator target_end(const box_type& b) {
    return make_transform_iterator(b.body_end(), body2source);
  }
  inline result_iterator result_begin(const box_type& b) {
    return make_transform_iterator(b.body_begin(), body2result);
  }
  inline result_iterator result_end(const box_type& b) {
    return make_transform_iterator(b.body_end(), body2result);
  }
};


//! Type hiding constructor for MinimalExecutor
template <typename Tree, typename Eval,
	  typename Kernel, typename SourceIter, typename Options>
ExecutorBase<Kernel>* make_minimal_executor(const Kernel& K,
					    SourceIter first, SourceIter last,
					    Options& opts) {
  return new MinimalExecutor<Kernel, Tree, Eval>(K, first, last, opts);
}



// TODO
  /*
    // This assumes the the body indices are consecutive within a box
    // The assumption is true for Octree,
    // but should be left for an "optimized" Evaluator that
    // explicitely makes this assumption
  inline charge_iterator charge_begin(const box_type& b) const {
    return charges_begin + b.body_begin()->index();
  }
  inline charge_iterator charge_end(const box_type& b) const {
    return charges_begin + b.body_end()->index();
  }
  */
  /*
    // This assumes the the body indices are consecutive within a box
    // The assumption is true for Octree,
    // but should be left for an "optimized" Evaluator that
    // explicitely makes this assumption
  inline result_iterator result_begin(const box_type& b) {
    return results_begin + b.body_begin()->index();
  }
  inline result_iterator result_end(const box_type& b) {
    return results_begin + b.body_end()->index();
  }
  */
