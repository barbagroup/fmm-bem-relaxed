#pragma once

template <typename Context>
struct EvaluatorBase {
  typedef Context                                 context_type;

  typedef typename context_type::kernel_type      kernel_type;
  typedef typename context_type::source_tree_type source_tree_type;
  typedef typename context_type::target_tree_type target_tree_type;

  virtual ~EvaluatorBase() {};
  virtual void execute(context_type&) const = 0;
};
