#pragma once

#include "ExecutorBase.hpp"

#include "ExecutorSingleTree.hpp"
#include "ExecutorDualTree.hpp"

#include "EvalUpward.hpp"
#include "EvalInteraction.hpp"
#include "EvalDownward.hpp"

#include "EvalLocal.hpp"
#include "EvalLocalSparse.hpp"
#include "EvalDiagonalSparse.hpp"

#include "EvalInteractionQueue.hpp"
#include "EvalInteractionLazy.hpp"
#include "EvalInteractionLazySparse.hpp"

#include "tree/Octree.hpp"



template <typename Executor, typename Options>
void make_evaluators(Executor& executor, Options& opts)
{
	if (opts.lazy_evaluation) {
    if (opts.sparse_local) {
      // sparse local evaluation
      auto lazy_eval = make_lazy_sparse_eval(executor, opts);
      executor.insert(lazy_eval);
    } else {
      // Custom lazy evaluator
      auto lazy_eval = make_lazy_eval(executor, opts);
      executor.insert(lazy_eval);
    }
  } else if (opts.local_evaluation) {
    // only evaluate local field for preconditioner
    if (opts.sparse_local) {
      auto sparse_eval = make_sparse_local_eval(executor, opts);
      executor.insert(sparse_eval);
    }
    else {
      auto local_eval = make_local_eval(executor, opts);
      executor.insert(local_eval);
    }
  } else if (opts.block_diagonal) {
    auto block_diagonal_eval = make_sparse_diagonal_eval(executor, opts);
    executor.insert(block_diagonal_eval);
	} else {
		// Standard evaluators
		auto upward = make_upward(executor, opts);
		executor.insert(upward);
		auto inter = make_interact(executor, opts);
		executor.insert(inter);
		auto downward = make_downward(executor, opts);
		executor.insert(downward);
	}
}

/** Single tree executor construction
 */
template <typename Kernel,
          typename SourceIter,
          typename Options>
// ExecutorBase<Kernel>* make_executor(const Kernel& K,
ExecutorSingleTree<Kernel,Octree<typename Kernel::point_type>>* make_executor(const Kernel& K,
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
