#pragma once

#include "Evaluator.hpp"

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
