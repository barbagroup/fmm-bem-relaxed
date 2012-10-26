#pragma once

#include "ExecutorBase.hpp"

#include "MinimalExecutor.hpp"

#include "EvalUpward.hpp"
#include "EvalInteraction.hpp"
#include "EvalDownward.hpp"

//! Type hiding factory for Executor based on FMMOptions
template <typename Tree, typename Kernel, typename Options>
ExecutorBase<Tree,Kernel>* make_executor(const Tree& tree,
                                         const Kernel& K,
                                         const Options& opts) {
  if (opts.evaluator == FMMOptions::FMM) {
    auto upward   = make_upward(tree, K, opts);
    auto inter    = make_fmm_inter(tree, K, opts);
    auto downward = make_downward(tree, K, opts);
    auto eval = make_evaluator(upward, inter, downward);
    return make_minimal_executor(tree, K, eval);
  } else if (opts.evaluator == FMMOptions::TREECODE) {
    auto upward = make_upward(tree, K, opts);
    auto inter  = make_tree_inter(tree, K, opts);
    auto eval = make_evaluator(upward, inter);
    return make_minimal_executor(tree, K, eval);
  }
  return NULL;
}
