#pragma once

// FMM includes
#include "FMMOptions.hpp"
#include "Vec.hpp"
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

	typedef typename kernel_type::source_type source_type;
	typedef typename kernel_type::target_type target_type;
	// TODO: Better point support?
	// Want all derived classes to use the following fmmplan::point_type wrapper?
	//typedef typename Vec<kernel_type::dimension, typename kernel_type::point_type> point_type;
	typedef typename kernel_type::charge_type charge_type;
	typedef typename kernel_type::result_type result_type;

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

	// EXECUTE

	/** Execute this FMM plan
	 */
	std::vector<result_type> execute(const std::vector<charge_type>& charges)
	{
		// Assert that source == target in FMMOptions

		if (!executor_) {
			printf("[E]: Executor not initialised -- returning..\n");
			return std::vector<result_type>(0);
		}

		std::vector<result_type> results(charges.size());
		executor_->execute(charges, results);

		// TODO: don't return this, provide accessor
		return results;
	}

	//private:
	ExecutorBase<kernel_type>* executor_;
	kernel_type K;
	FMMOptions opts_;

 private:

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
