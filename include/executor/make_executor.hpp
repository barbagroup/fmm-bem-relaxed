#pragma once

#include "ExecutorBase.hpp"

#include "MinimalExecutor.hpp"

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

  if (opts.evaluator == FMMOptions::TREECODE && opts.lazy_evaluation) {
    typedef typename evaluator_type<
        EvalInteractionLazy<Tree,FMMOptions::TREECODE>
        >::type Evaluator;

    return make_minimal_executor<Tree,Evaluator>(K, first, last, opts);
  } else  if (opts.evaluator == FMMOptions::FMM && opts.lazy_evaluation) {
    typedef typename evaluator_type<
        EvalInteractionLazy<Tree,FMMOptions::FMM>
        >::type Evaluator;

    return make_minimal_executor<Tree,Evaluator>(K, first, last, opts);
  } else if (opts.evaluator == FMMOptions::FMM) {
    typedef typename evaluator_type<
      EvalUpward,
      EvalInteraction<Tree,FMMOptions::FMM>,
      EvalDownward
      >::type Evaluator;

    return make_minimal_executor<Tree,Evaluator>(K, first, last, opts);
  } else if (opts.evaluator == FMMOptions::TREECODE) {
    typedef typename evaluator_type<
        EvalUpward,
        EvalInteraction<Tree,FMMOptions::TREECODE>
        >::type Evaluator;

    return make_minimal_executor<Tree,Evaluator>(K, first, last, opts);
  } else {
    return nullptr;
  }
}


template <typename Kernel, typename SourceIter, typename TargetIter,
	  typename Options>
ExecutorBase<Kernel>* make_executor(const Kernel& K,
				    SourceIter sfirst, SourceIter slast,
				    TargetIter tfirst, TargetIter tlast,
				    Options& opts) {
  typedef Octree<typename Kernel::point_type> Tree;

  if (opts.evaluator == FMMOptions::TREECODE && opts.lazy_evaluation) {
    typedef typename evaluator_type<
        EvalInteractionLazy<Tree,FMMOptions::TREECODE>
        >::type Evaluator;

    return make_minimal_executor<Tree,Evaluator>(K,
						 sfirst, slast,
						 tfirst, tlast,
						 opts);
  } else  if (opts.evaluator == FMMOptions::FMM && opts.lazy_evaluation) {
    typedef typename evaluator_type<
        EvalInteractionLazy<Tree,FMMOptions::FMM>
        >::type Evaluator;

    return make_minimal_executor<Tree,Evaluator>(K,
						 sfirst, slast,
						 tfirst, tlast,
						 opts);
  } else if (opts.evaluator == FMMOptions::FMM) {
    typedef typename evaluator_type<
      EvalUpward,
      EvalInteraction<Tree,FMMOptions::FMM>,
      EvalDownward
      >::type Evaluator;

    return make_minimal_executor<Tree,Evaluator>(K,
						 sfirst, slast,
						 tfirst, tlast,
						 opts);
  } else if (opts.evaluator == FMMOptions::TREECODE) {
    typedef typename evaluator_type<
        EvalUpward,
        EvalInteraction<Tree,FMMOptions::TREECODE>
        >::type Evaluator;

    return make_minimal_executor<Tree,Evaluator>(K,
						 sfirst, slast,
						 tfirst, tlast,
						 opts);
  } else {
    return NULL;
  }
}
