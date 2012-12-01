#pragma once

#include "ExecutorBase.hpp"

#include "SingleTreeExecutor.hpp"

#include "EvalUpward.hpp"
#include "EvalInteraction.hpp"
#include "EvalInteractionQueue.hpp"
#include "EvalInteractionLazy.hpp"
#include "EvalDownward.hpp"

#include "Octree.hpp"

template <typename Kernel, typename SourceIter, typename Options>
ExecutorBase<Kernel>* make_executor(const Kernel& K,
				    SourceIter first, SourceIter last,
				    Options& opts) {
  typedef Octree<typename Kernel::point_type> Tree;

  auto executor = make_executor<Tree>(K, first, last, opts);

  auto upward = make_upward(*executor, opts);
  executor->insert_eval(upward);
  auto inter = make_interact(*executor, opts);
  executor->insert_eval(inter);
  auto downward = make_downward(*executor, opts);
  executor->insert_eval(downward);

  return executor;
}

template <typename Kernel, typename SourceIter, typename TargetIter,
	  typename Options>
ExecutorBase<Kernel>* make_executor(const Kernel& K,
				    SourceIter sfirst, SourceIter slast,
				    TargetIter tfirst, TargetIter tlast,
				    Options& opts) {
  typedef Octree<typename Kernel::point_type> Tree;

  return nullptr;
}

