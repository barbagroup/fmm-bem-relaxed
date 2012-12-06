#pragma once

#include "ExecutorBase.hpp"

#include "ExecutorSingleTree.hpp"
#include "ExecutorDualTree.hpp"

#include "EvalUpward.hpp"
#include "EvalInteraction.hpp"
#include "EvalDownward.hpp"

#include "EvalInteractionQueue.hpp"
#include "EvalInteractionLazy.hpp"

#include "Octree.hpp"

template <typename Executor, typename Options>
void make_evaluators(Executor& executor, Options& opts)
{
	if (opts.lazy_evaluation) {
		// Custom lazy evaluator
		auto lazy_eval = make_lazy_eval(executor, opts);
		executor.insert_eval(lazy_eval);
	} else {
		// Standard evaluators
		auto upward = make_upward(executor, opts);
		executor.insert_eval(upward);
		auto inter = make_interact(executor, opts);
		executor.insert_eval(inter);
		auto downward = make_downward(executor, opts);
		executor.insert_eval(downward);
	}
}

/** Single tree executor construction
 */
template <typename Kernel,
          typename SourceIter,
          typename Options>
ExecutorBase<Kernel>* make_executor(const Kernel& K,
                                    SourceIter first, SourceIter last,
                                    Options& opts) {
  typedef Octree<typename Kernel::point_type> Tree;

  auto executor = make_executor<Tree>(K,
                                      first, last,
                                      opts);
  make_evaluators(*executor, opts);

  return executor;
}


/** Dual tree executor construction
 */
template <typename Kernel,
          typename SourceIter, typename TargetIter,
          typename Options>
ExecutorBase<Kernel>* make_executor(const Kernel& K,
                                    SourceIter sfirst, SourceIter slast,
                                    TargetIter tfirst, TargetIter tlast,
                                    Options& opts) {
  typedef Octree<typename Kernel::point_type> Tree;

  auto executor = make_executor<Tree>(K,
                                      sfirst, slast,
                                      tfirst, tlast,
                                      opts);
  make_evaluators(*executor, opts);

  return executor;
}
