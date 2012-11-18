#pragma once

#include "ExecutorBase.hpp"

#include "MinimalExecutor.hpp"

#include "EvalUpward.hpp"
#include "EvalInteraction.hpp"
#include "EvalInteractionQueue.hpp"
#include "EvalInteractionLazy.hpp"
#include "EvalDownward.hpp"

#include <tree/Octree.hpp>

template <typename Kernel, typename PointIter, typename Options>
ExecutorBase<Kernel>* make_executor(const Kernel& K,
				    PointIter first, PointIter last,
				    Options& opts) {
  typedef Octree<typename Kernel::point_type> Tree;

  if (opts.evaluator == FMMOptions::TREECODE && opts.lazy_evaluation) {
    typedef typename make_evaluator<
        EvalInteractionLazy<Tree,FMMOptions::TREECODE>
        >::type Evaluator;

    return make_minimal_executor<Tree,Evaluator>(K, first, last, opts);
  } else  if (opts.evaluator == FMMOptions::FMM && opts.lazy_evaluation) {
    typedef typename make_evaluator<
        EvalInteractionLazy<Tree,FMMOptions::FMM>
        >::type Evaluator;

    return make_minimal_executor<Tree,Evaluator>(K, first, last, opts);
  } else if (opts.evaluator == FMMOptions::FMM) {
    typedef typename make_evaluator<
      EvalUpward,
      EvalInteraction<Tree,FMMOptions::FMM>,
      EvalDownward
      >::type Evaluator;

    return make_minimal_executor<Tree,Evaluator>(K, first, last, opts);
  } else if (opts.evaluator == FMMOptions::TREECODE) {
    typedef typename make_evaluator<
        EvalUpward,
        EvalInteraction<Tree,FMMOptions::TREECODE>
        >::type Evaluator;

    return make_minimal_executor<Tree,Evaluator>(K, first, last, opts);
  } else {
    return NULL;
  }
}
