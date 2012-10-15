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


// Example evaluator
struct NullEval : public Evaluator<NullEval> {
  template <typename BoxContext>
  void execute(BoxContext&) const {}
};

// FMM / treecode upward sweep
template <typename Tree, typename Kernel, typename Options>
class EvalUpward : public Evaluator<EvalUpward<Tree,Kernel,Options>>
{
  const Tree& tree;
  const Kernel& K;

public:
  typedef typename Kernel::charge_type charge_type;
  typedef typename Kernel::result_type result_type;
  typedef typename Kernel::point_type  point_type;

  typedef typename Kernel::multipole_type multipole_type;
  typedef typename Kernel::local_type  local_type;

  //! constructor
  EvalUpward(const Tree& t, const Kernel& k, const Options& options)
      : tree(t), K(k) {
    (void) options;
  };

  template <typename BoxContext, typename BOX>
      void evalP2M(BoxContext& bc,
                   const BOX& box) const {
    auto p_begin = bc.point_begin(box);
    auto p_end   = bc.point_end(box);
    auto c_begin = bc.charge_begin(box);

    printf("P2M: box: %d\n", box.index());
    K.P2M(p_begin, p_end, c_begin,
          box.center(),
          bc.multipole_expansion(box));
  }

  template <typename BoxContext, typename BOX>
      void evalM2M(BoxContext& bc,
                   const BOX& cbox,
                   const BOX& box) const {
    printf("M2M: %d to %d\n", cbox.index(), box.index());
    K.M2M(bc.multipole_expansion(cbox),
          bc.multipole_expansion(box),
          box.center() - cbox.center());
  }

  template <typename BoxContext>
      void execute(BoxContext& bc) const {
    unsigned lowest_level = tree.levels();
    printf("lowest level in tree: %d\n",(int)lowest_level);

    // For the lowest level up to the highest level
    for (unsigned l = tree.levels()-1; l != 0; --l) {
      // For all boxes at this level
      auto b_end = tree.box_end(l);
      for (auto bit = tree.box_begin(l); bit != b_end; ++bit) {
        auto box = *bit;

        // Initialize box data
        double box_size = box.side_length();
        // K.init_multipole(M[idx], box_size);
        K.init_multipole(bc.multipole_expansion(box), box_size);
        // K.init_local(L[idx], box_size);
        K.init_local(bc.local_expansion(box), box_size);

        if (box.is_leaf()) {
          // If leaf, make P2M calls
          evalP2M(bc, box);

        } else {
          // If not leaf, make M2M calls

          // For all the children, M2M
          auto c_end = box.child_end();
          for (auto cit = box.child_begin(); cit != c_end; ++cit) {
            auto cbox = *cit;
            evalM2M(bc, cbox, box);
          }
        }
      }
    }
  }
};

template <typename Tree, typename Kernel, typename Options>
class EvalInteraction : public Evaluator<EvalInteraction<Tree,Kernel,Options>>
{
public:
  typedef typename Kernel::charge_type charge_type;
  typedef typename Kernel::result_type result_type;
  typedef typename Kernel::point_type  point_type;

  typedef typename Kernel::multipole_type multipole_type;
  typedef typename Kernel::local_type  local_type;

private:
  const Tree& tree;
  const Kernel& K;
  double THETA;

  /** One-sided P2P!!
   */
  template <typename BoxContext, typename BOX>
      void evalP2P(BoxContext& bc,
                   const BOX& b1,
                   const BOX& b2) const {
    // Point iters
    auto p1_begin = bc.point_begin(b1);
    auto p1_end   = bc.point_end(b1);
    auto p2_begin = bc.point_begin(b2);
    auto p2_end   = bc.point_end(b2);

    // Charge iters
    auto c1_begin = bc.charge_begin(b1);

    // Result iters
    auto r2_begin = bc.result_begin(b2);

    printf("P2P: %d to %d\n",b1.index(),b2.index());
    K.P2P(p1_begin, p1_end, c1_begin,
          p2_begin, p2_end,
          r2_begin);
  }

  template <typename BoxContext, typename BOX>
      void evalM2P(BoxContext& bc,
                   const BOX& b1,
                   const BOX& b2) const {
    // Target point iters
    auto t_begin = bc.point_begin(b2);
    auto t_end   = bc.point_end(b2);

    // Target result iters
    auto r_begin = bc.result_begin(b2);

    printf("M2P: %d to %d\n", b1.index(), b2.index());

    K.M2P(bc.multipole_expansion(b1),
          b1.center(),
          t_begin, t_end,
          r_begin);
  }

  template <typename BoxContext, typename BOX>
      void evalM2L(BoxContext& bc,
                   const BOX& b1,
                   const BOX& b2) const {
    printf("M2L: %d to %d\n", b2.index(), b1.index());

    K.M2L(bc.multipole_expansion(b1),
          bc.local_expansion(b2),
          b2.center() - b1.center());
  }

public:

  EvalInteraction(const Tree& t, const Kernel& k, const Options& options)
      : tree(t), K(k), THETA(options.THETA) {
    // any precomputation here
  }

  template <typename BoxContext, typename BOX, typename Q>
      void interact(BoxContext& bc, const BOX& b1, const BOX& b2, Q& pairQ) const {
    double r0_norm = norm(b1.center() - b2.center());
    if (r0_norm * THETA > b1.side_length()/2 + b2.side_length()/2) {
      // These boxes satisfy the multipole acceptance criteria
      evalM2L(bc, b1, b2);
    } else if(b1.is_leaf() && b2.is_leaf()) {
      evalP2P(bc, b2, b1);
    } else {
      pairQ.push_back(std::make_pair(b1,b2));
    }
  }

  template <typename BoxContext>
      void execute(BoxContext& bc) const {
    typedef typename Octree<point_type>::Box Box;
    typedef typename std::pair<Box, Box> box_pair;
    std::deque<box_pair> pairQ;

    // Queue based tree traversal for P2P, M2P, and/or M2L operations
    pairQ.push_back(box_pair(tree.root(), tree.root()));

    while (!pairQ.empty()) {
      auto b1 = pairQ.front().first;
      auto b2 = pairQ.front().second;
      pairQ.pop_front();

      if (b2.is_leaf() || (!b1.is_leaf() && b1.side_length() > b2.side_length())) {
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
};

template <typename Tree, typename Kernel, typename Options>
class EvalDownward : public Evaluator<EvalDownward<Tree,Kernel,Options>>
{
public:
  typedef typename Kernel::charge_type charge_type;
  typedef typename Kernel::result_type result_type;
  typedef typename Kernel::point_type  point_type;

  typedef typename Kernel::multipole_type multipole_type;
  typedef typename Kernel::local_type  local_type;

private:
  const Tree& tree;
  const Kernel& K;
public:

  EvalDownward(const Tree& t, const Kernel& k, const Options& options)
      : tree(t), K(k) {
    // any precomputation here
    (void) options;
  }

  template <typename BoxContext, typename BOX>
      void evalL2P(BoxContext& bc,
                   const BOX& box) const {
    auto t_begin = bc.point_begin(box);
    auto t_end   = bc.point_end(box);
    auto r_begin = bc.result_begin(box);

    printf("L2P: %d\n",box.index());
    K.L2P(bc.local_expansion(box), box.center(),
          t_begin, t_end,
          r_begin);
  }

  template <typename BoxContext, typename BOX>
      void evalL2L(BoxContext& bc,
                   const BOX& box,
                   const BOX& cbox) const {
    printf("L2L: %d to %d\n", box.index(), cbox.index());
    K.L2L(bc.local_expansion(box),
          bc.local_expansion(cbox),
          cbox.center() - box.center());
  }

  template <typename BoxContext>
      void execute(BoxContext& bc) const {
    // For the highest level down to the lowest level
    for (unsigned l = 1; l < tree.levels(); ++l) {
      // For all boxes at this level
      auto b_end = tree.box_end(l);
      for (auto bit = tree.box_begin(l); bit != b_end; ++bit) {
        auto box = *bit;

        // Initialize box data
        if (box.is_leaf()) {
          // If leaf, make L2P calls
          evalL2P(bc, box);
        } else {
          // If not leaf, make L2L calls

          // For all the children, L2L
          auto c_end = box.child_end();
          for (auto cit = box.child_begin(); cit != c_end; ++cit) {
            auto cbox = *cit;
            evalL2L(bc, box, cbox);
          }
        }
      }
    }
  }
};

