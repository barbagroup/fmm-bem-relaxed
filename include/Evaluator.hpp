#pragma once

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




template <typename DerivedType>
struct Evaluator {
};

template <typename E1, typename E2>
class EvalPair : public Evaluator<EvalPair<E1,E2>> {
  const E1* e1_;
  const E2* e2_;
public:
  EvalPair(const Evaluator<E1>* e1, const Evaluator<E2>* e2)
      : e1_(static_cast<const E1*>(e1)), e2_(static_cast<const E2*>(e2)) {
  }
  ~EvalPair() {
    delete e1_;
    delete e2_;
  }
  template <typename BoxContext>
      inline void execute(BoxContext& bc) const {
    e1_->execute(bc);
    e2_->execute(bc);
  }
};


template <typename E1, typename E2>
EvalPair<E1,E2>* make_evaluator(const Evaluator<E1>* e1,
                                const Evaluator<E2>* e2) {
  return new EvalPair<E1,E2>(e1, e2);
}

template <typename E1, typename E2, typename E3>
EvalPair<E1,EvalPair<E2,E3>>* make_evaluator(const Evaluator<E1>* e1,
                                             const Evaluator<E2>* e2,
                                             const Evaluator<E3>* e3) {
  return make_evaluator(e1, make_evaluator(e2, e3));
}
