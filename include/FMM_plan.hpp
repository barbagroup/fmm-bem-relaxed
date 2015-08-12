#pragma once

// FMM includes
#include "Vec.hpp"

#include "FMMOptions.hpp"
#include "KernelTraits.hpp"
#include "Logger.hpp"

#include "executor/make_executor.hpp"

//! global logging
Logger Log;

template <class Kernel>
class FMM_plan
{
 public:
	typedef Kernel kernel_type;

  typedef typename kernel_type::point_type point_type;
	typedef typename kernel_type::source_type source_type;
	typedef typename kernel_type::target_type target_type;
	// TODO: Better point support?
	// Want all derived classes to use the following fmmplan::point_type wrapper?
	//typedef typename Vec<kernel_type::dimension, typename kernel_type::point_type> point_type;
	typedef typename kernel_type::charge_type charge_type;
	typedef typename kernel_type::result_type result_type;
  // executor type
  typedef ExecutorSingleTree<kernel_type, Octree<point_type>> executor_type;

	// CONSTRUCTOR

	FMM_plan(const kernel_type& k,
	         const std::vector<source_type>& source,
	         FMMOptions& opts)
      : K(k), opts_(opts) {
		check_kernel();

		executor_ = make_executor(K,
		                          source.begin(), source.end(),
		                          opts_);
	}

	FMM_plan(const Kernel& k,
	         const std::vector<source_type>& source,
	         const std::vector<target_type>& target,
	         FMMOptions& opts)
      : K(k), opts_(opts) {
		check_kernel();

		executor_ = make_executor(K,
		                          source.begin(), source.end(),
		                          target.begin(), target.end(),
		                          opts_);
	}

	// DESTRUCTOR

	~FMM_plan() {
		delete executor_;
	}

  // MODIFIER

  kernel_type& kernel() {
    return K;
  }
  const kernel_type& kernel() const {
    return K;
  }

	// EXECUTE

	std::vector<result_type> execute(const std::vector<charge_type>& charges)
	{
		// Assert that source == target in FMMOptions

		if (!executor_) {
			printf("[E]: Executor not initialised -- returning..\n");
			return std::vector<result_type>(0);
		}

    // XXX: results.size == charges.size()?
		std::vector<result_type> results(charges.size());
		executor_->execute(charges, results);

		// TODO: don't return this, provide accessor
		return results;
	}

  /** Access to the Options this plan is operating with
   */
  FMMOptions& options() {
    return opts_;  // XXX: Need to update the plan with any new settings...
  }

  /** Access to iterators
   */
  typedef typename executor_type::body_source_iterator body_source_iterator;
  body_source_iterator source_begin() {
    return executor_->source_begin(executor_->source_tree().root());
  }

  body_source_iterator source_end() {
    return executor_->source_end(executor_->source_tree().root());
  }

 private:
	// ExecutorBase<kernel_type>* executor_;
  executor_type *executor_;
	kernel_type K;
	FMMOptions opts_;

	void check_kernel() {
		if (opts_.evaluator == FMMOptions::FMM &&
		    !ExpansionTraits<kernel_type>::is_valid_fmm) {
			std::cerr << ExpansionTraits<kernel_type>();
			std::cerr << "[W] Cannot use Kernel for FMM!\n";
		}

		if (opts_.evaluator == FMMOptions::TREECODE &&
		    !ExpansionTraits<kernel_type>::is_valid_treecode) {
			std::cerr << ExpansionTraits<kernel_type>();
			std::cerr << "[W] Cannot use Kernel for treecode!\n";
		}
	}
};
