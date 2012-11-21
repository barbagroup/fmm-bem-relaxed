#pragma once

#include "ExecutorBase.hpp"
#include "Evaluator.hpp"

#include <TransformIterator.hpp>
#include <type_traits>

/** @class Executor
 * @brief A very general Executor class. This provides a context to any tree
 * that provides the following interface:
 * Tree
 *   Tree(Iter, Iter, Options)          // Tree constructor
 *   int boxes() const                  // The number of boxes
 * Tree::box_type
 *   int index() const                  // Each box has a unique index
 *   point_type center() const          // Each box has a center point
 *   body_iterator body_begin() const
 *   body_iterator body_end() const     // Iterator pair to the box's bodies
 * Tree::body_iterator
 * Tree::body_type
 *   int index() const                  // Each body has a unique index
 *   int number() const                 // Original index of this body TODO TEMP
 * This class assumes nothing else about the tree.
 */
template <typename Kernel, typename Tree, typename Eval>
class MinimalExecutor : public ExecutorBase<Kernel>
{
 public:
  //! Kernel type
  typedef Kernel kernel_type;
  //! Tree type
  typedef Tree tree_type;
  //! Evaluator type
  typedef Eval eval_type;
  //! Tree box type
  typedef typename tree_type::box_type box_type;
  //! Tree body type
  typedef typename tree_type::body_type body_type;
  //! Kernel charge type
  typedef typename kernel_type::multipole_type multipole_type;
  //! Kernel result type
  typedef typename kernel_type::local_type local_type;
  //! Kernel point type
  typedef typename kernel_type::point_type point_type;
  //! Kernel source type
  typedef typename kernel_type::source_type source_type;
  //! Kernel target type
  typedef typename kernel_type::target_type target_type;
  //! Kernel charge type
  typedef typename kernel_type::charge_type charge_type;
  //! Kernel result type
  typedef typename kernel_type::result_type result_type;

protected:
  //! The tree of sources
  Tree source_tree_;
  //! Multipole expansions corresponding to Box indices in Tree
  std::vector<multipole_type> M_;
  //! Local expansions corresponding to Box indices in the Tree
  std::vector<local_type> L_;
  //! The sources associated with bodies in the source_tree
  std::vector<source_type> s_;

  //! Evaluator algorithm to apply
  Eval eval_;
  //! Reference to the Kernel
  const kernel_type& K_;

  template <typename T>
  struct body_map {
    typedef typename std::remove_const<T>::type vT;
    typedef typename std::vector<vT>::iterator iter;
    typedef typename std::vector<vT>::const_iterator const_iter;
    static constexpr bool is_const = std::is_const<T>::value;
    typename std::conditional<is_const, const_iter, iter>::type begin_;
    T& operator()(const body_type& b) const {
      return begin_[b.number()];  // TODO TEMP: number to workaround permute
    }
  };

  body_map<const source_type> body2source;
  body_map<const charge_type> body2charge;
  body_map<const target_type> body2target;
  body_map<      result_type> body2result;

  typedef typename tree_type::body_iterator              body_iterator;
  template <typename UnaryF>
  using transform_body_iterator = transform_iterator<body_iterator, UnaryF>;

  typedef transform_body_iterator<decltype(body2source)> source_iterator;
  typedef transform_body_iterator<decltype(body2charge)> charge_iterator;
  typedef transform_body_iterator<decltype(body2target)> target_iterator;
  typedef transform_body_iterator<decltype(body2result)> result_iterator;

public:
  //! Constructor
  template <typename SourceIter, typename Options>
  MinimalExecutor(const kernel_type& K,
		  SourceIter first, SourceIter last,
		  Options& opts)
    : source_tree_(first, last, opts),
      M_(source_tree_.boxes()),
      L_(source_tree_.boxes()),
      s_(first, last),
      K_(K) {
    if (opts.print_tree())
      std::cout << source_tree_ << std::endl;
  }

  template <typename Options>
  void create_eval(Options& opts) {
    eval_.create(K_, source_tree_, source_tree_, opts);
  }

  virtual void execute(const std::vector<charge_type>& charges,
		       std::vector<result_type>& results) {
    // TEMP Hack
    body2source.begin_ = s_.begin();
    body2charge.begin_ = charges.begin();
    body2target.begin_ = s_.begin();
    body2result.begin_ = results.begin();
    eval_.execute(*this);
  }

  const kernel_type& kernel() const {
    return K_;
  }

  tree_type& source_tree() {
    return source_tree_;
  }

  tree_type& target_tree() {
    return source_tree_;
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
    return make_transform_iterator(b.body_begin(), body2target);
  }
  inline target_iterator target_end(const box_type& b) {
    return make_transform_iterator(b.body_end(), body2target);
  }
  inline result_iterator result_begin(const box_type& b) {
    return make_transform_iterator(b.body_begin(), body2result);
  }
  inline result_iterator result_end(const box_type& b) {
    return make_transform_iterator(b.body_end(), body2result);
  }
};


template <typename Kernel, typename Tree, typename Eval>
class MinimalDualExecutor : public MinimalExecutor<Kernel, Tree, Eval>
{
  typedef MinimalExecutor<Kernel, Tree, Eval> super_type;

  typedef typename super_type::kernel_type kernel_type;
  typedef typename super_type::target_type target_type;
  typedef typename super_type::charge_type charge_type;
  typedef typename super_type::result_type result_type;

  //! The tree of targets
  Tree target_tree_;

  //! The targets associated with bodies in the target_tree
  std::vector<target_type> t_;

public:
  //! Constructor
  template <typename SourceIter, typename TargetIter, typename Options>
  MinimalDualExecutor(const kernel_type& K,
		      SourceIter sfirst, SourceIter slast,
		      TargetIter tfirst, TargetIter tlast,
		      Options& opts)
    : super_type(K, sfirst, slast, opts),
      target_tree_(tfirst, tlast, opts),
      t_(tfirst, tlast) {
    if (opts.print_tree())
      std::cout << target_tree_ << std::endl;
  }

  template <typename Options>
  void create_eval(Options& opts) {
    this->eval_.create(this->K_, this->source_tree_, target_tree_, opts);
  }

  virtual void execute(const std::vector<charge_type>& charges,
		       std::vector<result_type>& results) {
    // TEMP Hack
    this->body2source.begin_ = this->s_.begin();
    this->body2charge.begin_ = charges.begin();
    this->body2target.begin_ = this->t_.begin();
    this->body2result.begin_ = results.begin();
    this->eval_.execute(*this);
  }
};


//! Type hiding constructor for MinimalExecutor
template <typename Tree, typename Eval,
	  typename Kernel, typename SourceIter, typename Options>
ExecutorBase<Kernel>* make_minimal_executor(const Kernel& K,
					    SourceIter first, SourceIter last,
					    Options& opts) {
  typedef MinimalExecutor<Kernel, Tree, Eval> Executor;
  Executor* exec = new Executor(K, first, last, opts);
  exec->create_eval(opts);
  return exec;
}

template <typename Tree, typename Eval,
	  typename Kernel, typename SourceIter, typename TargetIter,
	  typename Options>
ExecutorBase<Kernel>* make_minimal_executor(const Kernel& K,
					    SourceIter sfirst, SourceIter slast,
					    TargetIter tfirst, TargetIter tlast,
					    Options& opts) {
  typedef MinimalDualExecutor<Kernel, Tree, Eval> Executor;
  Executor* exec = new Executor(K,
				sfirst, slast,
				tfirst, tlast,
				opts);
  exec->create_eval(opts);
  return exec;
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
