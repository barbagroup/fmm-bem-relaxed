#pragma once

#include "ExecutorBase.hpp"
#include "EvaluatorBase.hpp"

//#include "tree/TreeContext.hpp"
#include "executor/INITM.hpp"
#include "executor/INITL.hpp"

#include <type_traits>
#include <functional>

#include <boost/iterator/transform_iterator.hpp>
using boost::transform_iterator;
using boost::make_transform_iterator;

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

  //! Tree type
  typedef Tree tree_type;
  //! Source tree type
  typedef tree_type source_tree_type;
  //! Target tree type
  typedef tree_type target_tree_type;
  //! Tree box type
  typedef typename tree_type::box_type box_type;
  //! Tree box iterator
  typedef typename tree_type::box_iterator box_iterator;
  //! Tree body type
  typedef typename tree_type::body_type body_type;
  //! Tree body iterator
  typedef typename tree_type::body_iterator body_iterator;

  //! Kernel type
  typedef Kernel kernel_type;
  //! Kernel value type
  typedef typename Kernel::kernel_value_type kernel_value_type;
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
  // Transform a body to a value
  template <typename Indexable>
  struct BodyTransformer {
    typedef typename Indexable::reference result_type;
    typedef body_type argument_type;
    BodyTransformer(const Indexable& value) : value_(value) {}
    result_type operator()(const body_type& body) const {
      return value_[body.number()]; // TODO: TEMP to avoid permutation
    }
   private:
    Indexable value_;
  };

  template <typename Indexable>
  BodyTransformer<Indexable> body_transformer(Indexable v) const {
    return BodyTransformer<Indexable>(v);
  }

  template <typename Indexable>
  using body_transform = transform_iterator<BodyTransformer<Indexable>,
                                            body_iterator>;

  //! Reference to the Kernel
  const kernel_type& K_;

  //! The tree of sources
  tree_type source_tree_;
  //! The tree of targets
  tree_type target_tree_;
  //! Multipole acceptance
  std::function<bool(const box_type&, const box_type&)> acceptMultipole;

  //! Multipole expansions corresponding to Box indices in Tree
  typedef std::vector<multipole_type> multipole_container;
  multipole_container M_;
  //! Local expansions corresponding to Box indices in Tree
  typedef std::vector<local_type> local_container;
  local_container L_;

  //! The sources associated with bodies in the source_tree
  typedef const std::vector<source_type> source_container;
  typedef typename source_container::const_iterator source_iterator;
  source_container sources;
  source_iterator s_;
  //! The targets associated with bodies in the target_tree
  typedef const std::vector<target_type> target_container;
  typedef typename target_container::const_iterator target_iterator;
  target_container targets;
  target_iterator t_;

  //! Iterator to the start of the charge vector
  typedef std::vector<charge_type> charge_container;
  typedef typename charge_container::const_iterator charge_iterator;
  charge_iterator c_;
  //! Iterator to the start of the result vector
  typedef std::vector<result_type> result_container;
  typedef typename result_container::iterator result_iterator;
  result_iterator r_;

  //! Evaluator algorithms to apply
  EvaluatorCollection<self_type> evals_;

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
        acceptMultipole(opts.MAC()),
        M_(source_tree_.boxes()),
        L_((opts.evaluator == FMMOptions::TREECODE ? 0 : target_tree_.boxes())),
        sources(sfirst, slast),
        targets(tfirst, tlast) {
    s_ = sources.begin();
    t_ = targets.begin();
  }

  void insert(EvaluatorBase<self_type>* eval) {
    evals_.insert(eval);
  }

  virtual void execute(const std::vector<charge_type>& charges,
                       std::vector<result_type>& results) {
    s_ = sources.begin();
    c_ = charges.begin();
    t_ = targets.begin();
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
    return target_tree_;
  }

  // Accessors to make this Executor into a BoxContext
  inline multipole_type& multipole_expansion(const box_type& box) {
    return M_[box.index()];
  }
  inline const multipole_type& multipole_expansion(const box_type& box) const {
    return M_[box.index()];
  }
  inline local_type& local_expansion(const box_type& box) {
    return L_[box.index()];
  }
  inline const local_type& local_expansion(const box_type& box) const {
    return L_[box.index()];
  }

  inline point_type center(const box_type& b) const {
    return b.center();
  }

  typedef body_transform<source_iterator> body_source_iterator;
  inline body_source_iterator source_begin(const box_type& b) const {
    return make_transform_iterator(b.body_begin(), body_transformer(s_));
  }
  inline body_source_iterator source_end(const box_type& b) const {
    return make_transform_iterator(b.body_end(), body_transformer(s_));
  }

  typedef body_transform<charge_iterator> body_charge_iterator;
  inline body_charge_iterator charge_begin(const box_type& b) const {
    return make_transform_iterator(b.body_begin(), body_transformer(c_));
  }
  inline body_charge_iterator charge_end(const box_type& b) const {
    return make_transform_iterator(b.body_end(), body_transformer(c_));
  }

  typedef body_transform<target_iterator> body_target_iterator;
  inline body_target_iterator target_begin(const box_type& b) const {
    return make_transform_iterator(b.body_begin(), body_transformer(t_));
  }
  inline body_target_iterator target_end(const box_type& b) const {
    return make_transform_iterator(b.body_end(), body_transformer(t_));
  }

  typedef body_transform<result_iterator> body_result_iterator;
  inline body_result_iterator result_begin(const box_type& b) {
    return make_transform_iterator(b.body_begin(), body_transformer(r_));
  }
  inline body_result_iterator result_end(const box_type& b) {
    return make_transform_iterator(b.body_end(), body_transformer(r_));
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
