#pragma once

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
};



template <typename DerivedType>
struct Evaluator {
  const DerivedType& derived() const {
    return static_cast<const DerivedType&>(*this);
  }
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
  }
  ~Executor() {}

  void execute(ExpansionContext<Tree,Kernel>& ec,
               ChargeContext<Tree,Kernel>& cc) const {
    eval_.execute(ec, cc);
  }

  void execute(const std::vector<typename Kernel::charge_type>& charges,
               std::vector<typename Kernel::result_type>& results) {

  }
};

template <typename Tree, typename Kernel, typename E>
ExecutorBase<Tree,Kernel>* make_executor(const Tree& tree, const Kernel& K,
                                         const Evaluator<E>& eval) {
  return new Executor<Tree, Kernel, decltype(eval.derived())>(tree, K, eval.derived());
}


// Example evaluator
struct NullEval : public Evaluator<NullEval> {
  template <typename ExpansionContext, typename ChargeContext>
  void execute(ExpansionContext&, ChargeContext&) const {}
};



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
    // Do any precomputation
  };

  template <typename ExpansionContext, typename ChargeContext>
  void execute(ExpansionContext& ec, ChargeContext& cc) {
    (void) ec;
    (void) cc;
    // Translate upward(...) to use ec, cc
  }

  //! upward sweep
  void upward(std::vector<multipole_type>& M,
              std::vector<local_type>& L,
              const std::vector<charge_type>& charges)
  {
    // set charges_begin iterator
    typename std::vector<charge_type>::const_iterator charges_begin = charges.begin();

    M.resize(tree.boxes());
    L.resize(tree.boxes());

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
        K.init_multipole(M[idx], box_size);
        K.init_local(L[idx], box_size);

        if (box.is_leaf()) {
          // If leaf, make P2M calls
          auto body2point = [](typename Octree<point_type>::Body b) { return b.point(); };

          auto p_begin = make_transform_iterator(box.body_begin(), body2point);
          auto p_end   = make_transform_iterator(box.body_end(),   body2point);
          auto c_begin = charges.begin() + box.body_begin()->index();

          printf("P2M: box: %d\n", (int)box.index());
          K.P2M(p_begin, p_end, c_begin, box.center(), M[idx]);

        } else {
          // If not leaf, make M2M calls

          // For all the children, M2M
          auto c_end = box.child_end();
          for (auto cit = box.child_begin(); cit != c_end; ++cit) {
            auto cbox = *cit;
            auto translation = box.center() - cbox.center();

            printf("M2M: %d to %d\n", cbox.index(), idx);
            K.M2M(M[cbox.index()], M[idx], translation);
          }
        }
      }
    }
  }

};
