#pragma once

#include "Evaluator.hpp"

#include <TransformIterator.hpp>


template <typename Tree, typename Kernel>
class BoxContext {
  //! Tree type
  typedef Tree tree_type;
  //! Tree box type
  typedef typename Tree::box_type box_type;
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

  //! Multipole expansions corresponding to Box indices in Tree
  std::vector<multipole_type> M;
  //! Local expansions corresponding to Box indices in the Tree
  std::vector<local_type> L;
  //! Charges corresponding to bodies in the Tree
  std::vector<charge_type> charges;
  //! Results correspondint to bodies in the Tree
  std::vector<result_type> results;

  struct {
    const point_type& operator()(const typename Tree::body_type& b) const {
      return b.point();
    }
  } body2point;

 public:
  BoxContext(const Tree& tree)
      : M(tree.boxes()), L(tree.boxes()),
        charges(tree.bodies()), results(tree.bodies()) {
  }

  // Accessors
  inline multipole_type& multipole_expansion(const box_type& b) {
    return M[b.index()];
  }
  inline const multipole_type& multipole_expansion(const box_type& b) const {
    return M[b.index()];
  }
  inline local_type& local_expansion(const box_type& b) {
    return L[b.index()];
  }
  inline const multipole_type& local_expansion(const box_type& b) const {
    return L[b.index()];
  }

  inline auto point_begin(const box_type& b)
      -> decltype(make_transform_iterator(b.body_begin(), body2point)) {
    return make_transform_iterator(b.body_begin(), body2point);
  }
  inline auto point_end(const box_type& b)
      -> decltype(make_transform_iterator(b.body_end(), body2point)) {
    return make_transform_iterator(b.body_end(), body2point);
  }

  inline typename std::vector<charge_type>::const_iterator charge_begin(const box_type& b) const {
    return charges.begin() + b.body_begin()->index();
  }
  inline typename std::vector<charge_type>::const_iterator charge_end(const box_type& b) const {
    return charges.begin() + b.body_end()->index();
  }
  inline typename std::vector<result_type>::iterator result_begin(const box_type& b) {
    return results.begin() + b.body_begin()->index();
  }
  inline typename std::vector<result_type>::iterator result_end(const box_type& b) {
    return results.begin() + b.body_end()->index();
  }

  // TODO: What to do with these...
  template <typename Vector>
  inline void set_charges(Vector& c) {
    charges = c;
  }
  template <typename Vector>
  inline void set_results(Vector& r) {
    results = r;
  }
  inline std::vector<result_type> get_results() {
    return results;
  }
};



template <typename Tree, typename Kernel>
struct ExecutorBase {
  ExecutorBase() {}
  virtual ~ExecutorBase() {};
  // TODO: improve
  virtual void execute(const std::vector<typename Kernel::charge_type>& charges,
                       std::vector<typename Kernel::result_type>& results) = 0;
};


template <typename Tree, typename Kernel, typename E>
class Executor : public ExecutorBase<Tree,Kernel>
{
  BoxContext<Tree,Kernel> bc_;
  const E* eval_;

 public:
  Executor(const Tree& tree, const Kernel& K, const Evaluator<E>* eval)
      : bc_(tree), eval_(static_cast<const E*>(eval)) {
    (void) K;
  }
  ~Executor() {
    delete eval_;
  }

  void execute(const std::vector<typename Kernel::charge_type>& charges,
               std::vector<typename Kernel::result_type>& results) {
    bc_.set_charges(charges);
    bc_.set_results(results);
    eval_->execute(bc_);
    results = bc_.get_results();
  }
};


template <typename Tree, typename Kernel, typename E>
ExecutorBase<Tree,Kernel>* make_executor(const Tree& tree,
                                         const Kernel& K,
                                         const Evaluator<E>* eval) {
  return new Executor<Tree, Kernel, E>(tree, K, eval);
}
