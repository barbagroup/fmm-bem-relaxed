#pragma once

#include "ExecutorBase.hpp"

#include "ExecutorSingleTree.hpp"
#include "ExecutorDualTree.hpp"

#include "EvalUpward.hpp"
#include "EvalInteraction.hpp"
#include "EvalDownward.hpp"

#include "EvalLocal.hpp"

#include "EvalInteractionQueue.hpp"
#include "EvalInteractionLazy.hpp"

#include "Octree.hpp"

//#define TEMP

template <typename Executor, typename Options>
void make_evaluators(Executor& executor, Options& opts)
{
#ifndef TEMP
	if (opts.lazy_evaluation) {
		// Custom lazy evaluator
		auto lazy_eval = make_lazy_eval(executor, opts);
		executor.insert(lazy_eval);
  } else if (opts.local_evaluation) {
    // only evaluate local field for preconditioner
    auto local_eval = make_local_eval(executor, opts);
    executor.insert(local_eval);
	} else {
#endif
		// Standard evaluators
		auto upward = make_upward(executor, opts);
		executor.insert(upward);
		auto inter = make_interact(executor, opts);
		executor.insert(inter);
		auto downward = make_downward(executor, opts);
		executor.insert(downward);
#ifndef TEMP
	}
#endif
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
