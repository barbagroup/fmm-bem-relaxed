#pragma once

template <typename DerivedType>
struct Evaluator {
  template <typename Kernel, typename STree, typename TTree, typename Options>
  void create(const Kernel&, STree&, TTree&, Options&) {
    // Default Evaluator create does nothing
  }
};

template <typename E1, typename E2>
class EvalPair : public Evaluator<EvalPair<E1,E2>> {
  E1 e1_;
  E2 e2_;
public:
  EvalPair() {
  }
  EvalPair(const E1& e1, const E2& e2)
  : e1_(e1), e2_(e2) {
  }
  template <typename Kernel, typename STree, typename TTree, typename Options>
  void create(const Kernel& K, STree& stree, TTree& ttree, Options& opts) {
    e1_.create(K, stree, ttree, opts);
    e2_.create(K, stree, ttree, opts);
  }
  template <typename BoxContext>
  inline void execute(BoxContext& bc) const {
    e1_.execute(bc);
    e2_.execute(bc);
  }
};

template <typename... Es>
struct evaluator_type;

template <typename E, typename... Es>
struct evaluator_type<E, Es...> {
  typedef EvalPair<E, typename evaluator_type<Es...>::type> type;
};

template <typename E>
struct evaluator_type<E> {
  typedef E type;
};


template <typename E1, typename E2>
EvalPair<E1,E2> make_evaluator(const Evaluator<E1>& e1,
			       const Evaluator<E2>& e2) {
  return new EvalPair<E1,E2>(e1, e2);
}

template <typename E1, typename E2, typename E3>
EvalPair<E1,EvalPair<E2,E3>> make_evaluator(const Evaluator<E1>& e1,
					    const Evaluator<E2>& e2,
					    const Evaluator<E3>& e3) {
  return make_evaluator(e1, make_evaluator(e2, e3));
}
