#pragma once

#include <TransformIterator.hpp>

template <typename Tree, typename Kernel>
class ExpansionContext {
  //! Tree box type
  typedef typename Tree::box_type box_type;
  //! Kernel charge type
  typedef typename Kernel::multipole_type multipole_type;
  //! Kernel result type
  typedef typename Kernel::local_type local_type;

  //! Multipole expansions corresponding to Box indices in Octree
  std::vector<multipole_type> M;
  //! Local expansions corresponding to Box indices in Octree
  std::vector<local_type> L;
 public:
  ExpansionContext(unsigned size)
      : M(size), L(size) {
  }
  // resize the expansions
  void resize(unsigned size)
  {
    printf("resizing expansions\n");
    M.resize(size);
    printf("M resized\n");
    L.resize(size);
    printf("L resized\n");
  }

  // Accessors
  multipole_type& multipole_expansion(const box_type& b) {
    return M[b.index()];
  }
  const multipole_type& multipole_expansion(const box_type& b) const {
    return M[b.index()];
  }
  local_type& local_expansion(const box_type& b) {
    return L[b.index()];
  }
  const multipole_type& local_expansion(const box_type& b) const {
    return L[b.index()];
  }
};

template <typename Tree, typename Kernel>
class ChargeContext {
  //! Tree box type
  typedef typename Tree::box_type box_type;
  //! Kernel charge type
  typedef typename Kernel::charge_type charge_type;
  //! Kernel result type
  typedef typename Kernel::result_type result_type;

  std::vector<charge_type> charges;
  std::vector<result_type> results;
 public:
  ChargeContext(unsigned size)
    : charges(size), results(size) {
  }
  // Accessors
  typename std::vector<charge_type>::const_iterator charge_begin(const box_type& b) const {
    return charges.begin() + b.body_begin()->index();
  }
  typename std::vector<charge_type>::const_iterator charge_end(const box_type& b) const {
    return charges.begin() + b.body_end()->index();
  }
  typename std::vector<result_type>::iterator result_begin(const box_type& b) {
    return results.begin() + b.body_begin()->index();
  }
  typename std::vector<result_type>::iterator result_end(const box_type& b) {
    return results.begin() + b.body_end()->index();
  }
  template <typename Vector>
  void set_charges(Vector& c) {
    charges = c;
  }
  template <typename Vector>
  void set_results(Vector& r) {
    results = r;
  }
};



template <typename DerivedType>
struct Evaluator {
  const DerivedType& derived() const {
    return static_cast<const DerivedType&>(*this);
  }
 private:
  Evaluator(const Evaluator<DerivedType>& other) {(void)other;};
  Evaluator<DerivedType>& operator=(const Evaluator<DerivedType>& other) {(void)other;};
};

template <typename E1, typename E2>
class EvalPair : public Evaluator<EvalPair<E1,E2>> {
  const E1& e1_;
  const E2& e2_;
public:
  EvalPair(const E1& e1, const E2& e2)
      : e1_(e1), e2_(e2) {
  }
  template <typename ExpansionContext, typename ChargeContext>
  void execute(ExpansionContext& ec, ChargeContext& cc) const {
    e1_.execute(ec, cc);
    e2_.execute(ec, cc);
  }
};

template <typename E1, typename E2>
EvalPair<E1,E2> make_evaluator(const Evaluator<E1>& e1,
                               const Evaluator<E2>& e2) {
  return EvalPair<E1,E2>(e1.derived(), e2.derived());
}

template <typename E1, typename E2, typename E3>
EvalPair<E1,EvalPair<E2,E3>> make_evaluator(const Evaluator<E1>& e1,
                                            const Evaluator<E2>& e2,
                                            const Evaluator<E3>& e3) {
  return make_evaluator(e1.derived(), make_evaluator(e2.derived(), e3.derived()));
}

template <typename Tree, typename Kernel>
struct ExecutorBase {
  virtual ~ExecutorBase() {};
  virtual void execute(ExpansionContext<Tree,Kernel>& ec,
                       ChargeContext<Tree,Kernel>& cc) const = 0;
  // TODO: remove
  virtual void execute(const std::vector<typename Kernel::charge_type>& charges,
                       std::vector<typename Kernel::result_type>& results) = 0;
};


template <typename Tree, typename Kernel, typename E>
class Executor : public ExecutorBase<Tree,Kernel>
{
  const E& eval_;
  ExpansionContext<Tree,Kernel> ec_;
  ChargeContext<Tree,Kernel> cc_;

 public:
  Executor(const Tree& tree, const Kernel& K, const E& eval)
      : eval_(eval), ec_(tree.boxes()), cc_(tree.bodies()) {
    (void) K;
  }
  ~Executor() {}

  void execute(ExpansionContext<Tree,Kernel>& ec,
               ChargeContext<Tree,Kernel>& cc) const {
    printf("execute(ec, cc) called\n");
    eval_.execute(ec, cc);
  }

  void execute(const std::vector<typename Kernel::charge_type>& charges,
               std::vector<typename Kernel::result_type>& results) {
      (void) charges;
      (void) results;
      printf("execute(charges, results) called\n");
      cc_.set_charges(charges);
      cc_.set_results(results);
      eval_.execute(ec_,cc_);
  }
};

template <typename Tree, typename Kernel, typename E>
ExecutorBase<Tree,Kernel>* make_executor(const Tree& tree, const Kernel& K,
                                         const Evaluator<E>& eval) {
  printf("ExecutorBase::make_executor: &tree: %p\n",&tree);
  printf("ExecutorBase::make_executor: tree.bodies(): %d\n",(int)tree.bodies());
  return new Executor<Tree, Kernel, decltype(eval.derived())>(tree, K, eval.derived());
}



// Example evaluator
struct NullEval : public Evaluator<NullEval> {
  template <typename ExpansionContext, typename ChargeContext>
  void execute(ExpansionContext&, ChargeContext&) const {}
};

// FMM / treecode upward sweep
template <typename Tree, typename Kernel, typename Options>
class EvalUpward : public Evaluator<EvalUpward<Tree,Kernel,Options>>
{
  const Tree& tree;
  const Kernel& K;

 private:
  EvalUpward(const EvalUpward<Tree,Kernel,Options>& other) {(void)other;};
  EvalUpward<Tree,Kernel,Options>& operator=(const EvalUpward<Tree,Kernel,Options>& other) {(void)other;};
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
    // Do any precomputation
    printf("EvalUpward::EvalUpward(): &{t,tree}: %p, %p\n",&t,&tree);
    printf("EvalUpward::EvalUpward(): {t,tree}.bodies(): %d, %d\n",(int)t.bodies(),(int)tree.bodies());
  };

  template <typename ExpansionContext, typename ChargeContext>
  void execute(ExpansionContext& ec, ChargeContext& cc) const {
    printf("EvalUpward executing...\n");
    // Translate upward(...) to use ec, cc
    printf("EvalUpward::execute(): &tree: %p\n",&tree);
    printf("EvalUpward::execute(): tree.bodies(): %d\n",(int)tree.bodies());
    ec.resize(tree.boxes());
    printf("ec resized\n");

    // EvaluatorBase<Tree,Kernel>::M.resize(10);
    unsigned lowest_level = tree.levels();
    printf("lowest level in tree: %d\n",(int)lowest_level);

    // For the lowest level up to the highest level
    for (unsigned l = tree.levels()-1; l != 0; --l) {
      // For all boxes at this level
      auto b_end = tree.box_end(l);
      for (auto bit = tree.box_begin(l); bit != b_end; ++bit) {
        auto box = *bit;

        // Initialize box data
        unsigned idx = box.index();
        double box_size = box.side_length();
        // K.init_multipole(M[idx], box_size);
        K.init_multipole(ec.multipole_expansion(box), box_size);
        // K.init_local(L[idx], box_size);
        K.init_local(ec.local_expansion(box), box_size);

        if (box.is_leaf()) {
          // If leaf, make P2M calls
          auto body2point = [](typename Octree<point_type>::Body b) { return b.point(); };

          auto p_begin = make_transform_iterator(box.body_begin(), body2point);
          auto p_end   = make_transform_iterator(box.body_end(),   body2point);
          // auto c_begin = charges.begin() + box.body_begin()->index();
          auto c_begin = cc.charge_begin(box);

          printf("P2M: box: %d\n", (int)box.index());
          // K.P2M(p_begin, p_end, c_begin, box.center(), M[idx]);
          K.P2M(p_begin, p_end, c_begin, box.center(), ec.multipole_expansion(box));

        } else {
          // If not leaf, make M2M calls

          // For all the children, M2M
          auto c_end = box.child_end();
          for (auto cit = box.child_begin(); cit != c_end; ++cit) {
            auto cbox = *cit;
            auto translation = box.center() - cbox.center();

            printf("M2M: %d to %d\n", cbox.index(), idx);
            // K.M2M(M[cbox.index()], M[idx], translation);
            K.M2M(ec.multipole_expansion(cbox),ec.multipole_expansion(box), translation);
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
  template <typename BOX, typename ChargeContext>
  void evalP2P(const BOX& b1, 
               const BOX& b2,
               ChargeContext& cc) const {
    // Point iters
    auto body2point = [](const typename Octree<point_type>::Body& b) { return b.point(); };
    auto p1_begin = make_transform_iterator(b1.body_begin(), body2point);
    auto p1_end   = make_transform_iterator(b1.body_end(),   body2point);
    auto p2_begin = make_transform_iterator(b2.body_begin(), body2point);
    auto p2_end   = make_transform_iterator(b2.body_end(),   body2point);

    // Charge iters
    //auto c1_begin = charges_begin + b1.body_begin()->index();
    auto c1_begin = cc.charge_begin(b1);

    // Result iters
    //auto r2_begin = results_begin + b2.body_begin()->index();
    auto r2_begin = cc.result_begin(b2);

    printf("P2P: %d to %d\n",b1.index(),b2.index());

    K.P2P(p1_begin, p1_end, c1_begin,
          p2_begin, p2_end,
          r2_begin);
  }

  template <typename BOX, typename ExpansionContext, typename ChargeContext>
  void evalM2P(const BOX& b1, 
               const BOX& b2,
                     ExpansionContext& ec,
                     ChargeContext& cc) const
  {
    // Target point iters
    auto body2point = [](const typename Octree<point_type>::Body& b) { return b.point(); };
    auto t_begin = make_transform_iterator(b2.body_begin(), body2point);
    auto t_end   = make_transform_iterator(b2.body_end(), body2point);

    // Target result iters
    // auto r_begin = results_begin + b2.body_begin()->index();
    auto r_begin = cc.result_begin(b2);

    printf("M2P: %d to %d\n", b1.index(), b2.index());

    // K.M2P(M[b1.index()], b1.center(),
    K.M2P(ec.multipole_expansion(b1), b1.center(),
          t_begin, t_end,
          r_begin);
  }

  template <typename ExpansionContext>
  void evalM2L(const typename Octree<point_type>::Box& b1, 
               const typename Octree<point_type>::Box& b2,
                              ExpansionContext& ec) const 
  {
    // auto translation = b1.center() - b2.center();
    auto translation = b2.center() - b1.center();

    printf("M2L: %d to %d\n",b2.index(),b1.index());

    K.M2L(ec.multipole_expansion(b1), // M[b1.index()],
          ec.local_expansion(b2), // L[b2.index()],
          translation);
  }

 public:

  EvalInteraction(const Tree& t, const Kernel& k, const Options& options)
        : tree(t), K(k), THETA(options.THETA) {
    // any precomputation here
  }

  template <typename BOX, typename Q, typename ExpansionContext, typename ChargeContext>
  void interact(const BOX& b1, const BOX& b2, Q& pairQ,
                ExpansionContext& ec, ChargeContext& cc) const {
    double r0_norm = norm(b1.center() - b2.center());
    if (r0_norm * THETA > b1.side_length()/2 + b2.side_length()/2) {
      // These boxes satisfy the multipole acceptance criteria
      evalM2L(b1,b2,ec);
    } else if(b1.is_leaf() && b2.is_leaf()) {
      evalP2P(b2,b1,cc);
    } else {
      pairQ.push_back(std::make_pair(b1,b2));
    }
  }

  template <typename ExpansionContext, typename ChargeContext>
  void execute(ExpansionContext& ec, ChargeContext& cc) const {

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
          interact(*cit, b2, pairQ,ec,cc);
      } else {
        // Split the second box into children and interact
        auto c_end = b2.child_end();
        for (auto cit = b2.child_begin(); cit != c_end; ++cit)
          interact(b1, *cit, pairQ,ec,cc);
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
  template <typename ExpansionContext, typename ChargeContext>
  void execute(ExpansionContext& ec, ChargeContext& cc) const {
    // For the highest level down to the lowest level
    for (unsigned l = 1; l < tree.levels(); ++l) {
      // For all boxes at this level
      auto b_end = tree.box_end(l);
      for (auto bit = tree.box_begin(l); bit != b_end; ++bit) {
        auto box = *bit;
        unsigned idx = box.index();

        // Initialize box data
        if (box.is_leaf()) {
          // If leaf, make L2P calls

          // For all the bodies, L2P
          auto body2point = [](typename Octree<point_type>::Body b) { return b.point(); };
          auto t_begin = make_transform_iterator(box.body_begin(), body2point);
          auto t_end   = make_transform_iterator(box.body_end(), body2point);
          auto r_begin = cc.result_begin(box); // results_begin + box.body_begin()->index();

          printf("L2P: %d\n",idx);
          K.L2P(ec.local_expansion(box), box.center(),
                t_begin, t_end,
                r_begin);
        } else {
          // If not leaf, make L2L calls

          // For all the children, L2L
          auto c_end = box.child_end();
          for (auto cit = box.child_begin(); cit != c_end; ++cit) {
            auto cbox = *cit;
            auto translation = cbox.center() - box.center();

            printf("L2L: %d to %d\n",idx,cbox.index());
            // K.L2L(L[idx], L[cbox.index()], translation);
            K.L2L(ec.local_expansion(box), ec.local_expansion(cbox), translation);
          }   
        }   
      }   
    }
  }
};

