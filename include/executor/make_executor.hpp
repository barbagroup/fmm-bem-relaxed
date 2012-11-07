#pragma once

#include "ExecutorBase.hpp"

#include "MinimalExecutor.hpp"

#include "EvalUpward.hpp"
#include "EvalInteraction.hpp"
#include "EvalInteractionQueue.hpp"
#include "EvalDownward.hpp"

#include <tree/Octree.hpp>

template <typename Kernel, typename PointIter, typename Options>
ExecutorBase<Kernel>* make_executor(const Kernel& K,
				    PointIter first, PointIter last,
				    Options& opts) {
  typedef Octree<typename Kernel::source_type,
		 typename Kernel::point_type> Tree;

  if (opts.evaluator == FMMOptions::FMM) {
    typedef typename make_evaluator<EvalUpward<Kernel,Tree>,
    EvalInteraction<Kernel,Tree,FMMOptions::FMM>,
    EvalDownward<Kernel,Tree>>::type Evaluator;

    return make_minimal_executor<Tree,Evaluator>(K, first, last, opts);
  } else if (opts.evaluator == FMMOptions::TREECODE) {
    typedef typename make_evaluator<EvalUpward<Kernel,Tree>,
    EvalInteraction<Kernel,Tree,FMMOptions::TREECODE>>::type Evaluator;

    return make_minimal_executor<Tree,Evaluator>(K, first, last, opts);
  }
    return NULL;
}
