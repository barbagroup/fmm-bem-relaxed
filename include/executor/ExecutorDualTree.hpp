#pragma once

#include "ExecutorBase.hpp"
#include "EvaluatorBase.hpp"

#include "tree/TreeContext.hpp"

#include <type_traits>
#include <functional>

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
template <typename Kernel, typename Tree>
class ExecutorDualTree : public ExecutorBase<Kernel>
{
 public:
  //! This type
  typedef ExecutorDualTree<Kernel,Tree> self_type;
  //! Kernel type
  typedef Kernel kernel_type;
  //! Tree type
  typedef Tree tree_type;
  //! Source tree type
  typedef tree_type source_tree_type;
  //! Target tree type
  typedef tree_type target_tree_type;
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
  //! Reference to the Kernel
  const kernel_type& K_;

  //! The tree of sources
  tree_type source_tree_;
  //! The tree of targets
  tree_type target_tree_;

  // TODO: Fix const correctness

  //! Multipole expansions corresponding to Box indices in Tree
  BoxMap<tree_type, std::vector<multipole_type>> M_;
  //! Local expansions corresponding to Box indices in Tree
  BoxMap<tree_type, std::vector<local_type>> L_;
  //! The sources associated with bodies in the source_tree
  BodyMap<tree_type, std::vector<source_type>> s_;
  //! The charges associated with bodies in the source_tree
  BodyMap<tree_type, typename std::vector<charge_type>::const_iterator> c_;
  //! The targets associated with bodies in the source_tree
  BodyMap<tree_type, std::vector<target_type>> t_;
  //! The results associated with bodies in the source_tree
  BodyMap<tree_type, typename std::vector<result_type>::iterator> r_;

  //! Evaluator algorithms to apply
  EvaluatorCollection<self_type> evals_;

  //! Multipole acceptance
  std::function<bool(const box_type&, const box_type&)> acceptMultipole;

  //! is treecode evaluator being used
  bool isTreecode;

 public:
  //! Constructor
  template <typename SourceIter, typename TargetIter, typename Options>
  ExecutorDualTree(const kernel_type& K,
                   SourceIter sfirst, SourceIter slast,
                   TargetIter tfirst, TargetIter tlast,
                   Options& opts)
      : K_(K),
        source_tree_(sfirst, slast, opts),
        target_tree_(tfirst, tlast, opts),
        M_(std::vector<multipole_type>(source_tree_.boxes())),
        L_(std::vector<local_type>(opts.evaluator == FMMOptions::TREECODE ?
                                   0 : source_tree_.boxes())),
        s_(std::vector<source_type>(sfirst, slast)),
    t_(std::vector<target_type>(tfirst, tlast)),
    acceptMultipole(opts.MAC()),
    isTreecode(opts.evaluator == FMMOptions::TREECODE ? true : false) {
  }

  void insert(EvaluatorBase<self_type>* eval) {
    evals_.insert(eval);
  }

  virtual void execute(const std::vector<charge_type>& charges,
                       std::vector<result_type>& results) {
    c_ = charges.begin();
    r_ = results.begin();
    evals_.execute(*this);
  }

  bool accept_multipole(const box_type& source, const box_type& target) const {
    return acceptMultipole(source, target);
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

  // Re-initialise all expansions
  void reset_expansions()
  {
    for (auto it=this->source_tree().box_begin(); it!=this->source_tree().box_end(); ++it) {
      INITM::eval(this->kernel(), *this, *it);
      if (!isTreecode) {
        INITL::eval(this->kernel(), *this, *it);
      }
    }
  }

  // Accessors to make this Executor into a BoxContext
  inline multipole_type& multipole_expansion(const box_type& box) {
    return M_(box);
  }
  inline const multipole_type& multipole_expansion(const box_type& box) const {
    return M_(box);
  }
  inline local_type& local_expansion(const box_type& box) {
    return L_(box);
  }
  inline const local_type& local_expansion(const box_type& box) const {
    return L_(box);
  }

  inline point_type center(const box_type& b) const {
    return b.center();
  }
  inline double box_size(const box_type& b) const {
    return b.side_length();
  }

  typedef typename decltype(s_)::body_value_const_iterator source_iterator;
  inline source_iterator source_begin(const box_type& b) const {
    return s_.begin(b);
  }
  inline source_iterator source_end(const box_type& b) const {
    return s_.end(b);
  }
  typedef typename decltype(c_)::body_value_const_iterator charge_iterator;
  inline charge_iterator charge_begin(const box_type& b) const {
    return c_.begin(b);
  }
  inline charge_iterator charge_end(const box_type& b) const {
    return c_.end(b);
  }
  typedef typename decltype(t_)::body_value_const_iterator target_iterator;
  inline target_iterator target_begin(const box_type& b) const {
    return t_.begin(b);
  }
  inline target_iterator target_end(const box_type& b) const {
    return t_.end(b);
  }
  typedef typename decltype(r_)::body_value_iterator result_iterator;
  inline result_iterator result_begin(const box_type& b) {
    return r_.begin(b);
  }
  inline result_iterator result_end(const box_type& b) {
    return r_.end(b);
  }
};


template <typename Tree, typename Kernel,
          typename SourceIter, typename TargetIter,
          typename Options>
ExecutorDualTree<Kernel,Tree>* make_executor(const Kernel& K,
                                             SourceIter sfirst, SourceIter slast,
                                             TargetIter tfirst, TargetIter tlast,
                                             Options& opts) {
  return new ExecutorDualTree<Kernel,Tree>(K,
                                           sfirst, slast,
                                           tfirst, tlast,
                                           opts);
}
