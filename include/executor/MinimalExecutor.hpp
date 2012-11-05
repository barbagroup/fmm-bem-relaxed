#pragma once

#include "ExecutorBase.hpp"
#include "Evaluator.hpp"
#include <TransformIterator.hpp>

/** @class Executor
 * @brief A very general Executor class. This provides a context to any tree
 * that provides the following interface:
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
template <typename Tree, typename Kernel, typename E>
class MinimalExecutor : public ExecutorBase<Tree,Kernel>
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
  //! Kernel charge type
  typedef typename Kernel::charge_type charge_type;
  //! Kernel result type
  typedef typename Kernel::result_type result_type;

private:
  //! Multipole expansions corresponding to Box indices in Tree
  std::vector<multipole_type> M_;
  //! Local expansions corresponding to Box indices in the Tree
  std::vector<local_type> L_;
  //! Evaluator algorithm to apply
  const E* eval_;

  //! Anonymous struct to convert body_type to point_type (TODO: Use lambdas?)
  struct {
    const point_type& operator()(const body_type& b) const {
      return b.point();
    }
  } body2point;

  //! Charges corresponding to bodies in the Tree
  typedef typename std::vector<charge_type>::const_iterator charge_iterator;
  //! Anonymous struct to convert body_type to charge_type (TODO: Use lambdas?)
  struct {
    charge_iterator charges_begin;
    const charge_type& operator()(const body_type& b) const {
      return charges_begin[b.index()];
    }
  } body2charge;

  //! Results corresponding to bodies in the Tree
  typedef typename std::vector<result_type>::iterator result_iterator;
  //! Anonymous struct to convert body_type to result_type (TODO: Use lambdas?)
  struct {
    result_iterator results_begin;
    result_type& operator()(const body_type& b) const {
      return results_begin[b.index()];
    }
  } body2result;

public:
  //! Constructor
  MinimalExecutor(const Tree& tree, const Kernel& K, const Evaluator<E>* eval)
      : M_(tree.boxes()), L_(tree.boxes()),
	eval_(static_cast<const E*>(eval)) {
    (void) K;
  }
  // Destructor
  ~MinimalExecutor() {
    delete eval_;
  }

  // Virtual function to run this Executor
  void execute(const std::vector<typename Kernel::charge_type>& charges,
               std::vector<typename Kernel::result_type>& results) {
    body2charge.charges_begin = charges.begin();
    body2result.results_begin = results.begin();
    eval_->execute(*this);
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
  inline auto point_begin(const box_type& b)
      -> decltype(make_transform_iterator(b.body_begin(), body2point)) {
    return make_transform_iterator(b.body_begin(), body2point);
  }
  inline auto point_end(const box_type& b)
      -> decltype(make_transform_iterator(b.body_end(), body2point)) {
    return make_transform_iterator(b.body_end(), body2point);
  }
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
  inline auto charge_begin(const box_type& b)
      -> decltype(make_transform_iterator(b.body_begin(), body2charge)) {
    return make_transform_iterator(b.body_begin(), body2charge);
  }
  inline auto charge_end(const box_type& b)
      -> decltype(make_transform_iterator(b.body_end(), body2charge)) {
    return make_transform_iterator(b.body_end(), body2charge);
  }
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
  inline auto result_begin(const box_type& b)
      -> decltype(make_transform_iterator(b.body_begin(), body2result)) {
    return make_transform_iterator(b.body_begin(), body2result);
  }
  inline auto result_end(const box_type& b)
      -> decltype(make_transform_iterator(b.body_end(), body2result)) {
    return make_transform_iterator(b.body_end(), body2result);
  }
};


//! Type hiding constructor for MinimalExecutor
template <typename Tree, typename Kernel, typename E>
ExecutorBase<Tree,Kernel>* make_minimal_executor(const Tree& tree,
						 const Kernel& K,
						 const Evaluator<E>* eval) {
  return new MinimalExecutor<Tree, Kernel, E>(tree, K, eval);
}
