#pragma once

// FMM includes
#include <FMMOptions.hpp>
#include <Vec.hpp>

#include <Logger.hpp>

#include <executor/make_executor.hpp>

//! global logging
Logger Log;

template <class Kernel>
class FMM_plan
{
 public:
  typedef typename Kernel::point_type point_type;
  typedef typename Kernel::source_type source_type;
  typedef typename Kernel::target_type target_type;
  // TODO: Better point support?
  // Want all derived classes to use the following fmmplan::point_type wrapper?
  //typedef typename Vec<Kernel::dimension, typename Kernel::point_type> point_type;
  typedef typename Kernel::charge_type charge_type;
  typedef typename Kernel::result_type result_type;

  typedef Kernel kernel_type;

  //private:
  ExecutorBase<kernel_type>* executor_;
  Kernel K;
  FMMOptions opts_;

public:

  // CONSTRUCTOR

  FMM_plan(const Kernel& k,
	   const std::vector<source_type>& source,
           FMMOptions& opts)
    : K(k), opts_(opts) {
    executor_ = make_executor(K, source.begin(), source.end(), opts_);
  }

  // DESTRUCTOR

  ~FMM_plan() {
    delete executor_;
  }

  // EXECUTE

  /** Set the executor strategy of this plan at runtime
   */
  //void set_options(FMMOptions& options) {
  //  executor = make_executor(otree, K, options);
  //}

  /** Execute this FMM plan
   */
  std::vector<result_type> execute(const std::vector<charge_type>& charges,
                                   const std::vector<target_type>& t_points)
  {
    if (!executor_) {
      printf("[E]: Executor not initialised -- returning..\n");
      return std::vector<result_type>(0);
    }

    (void) t_points; // Quiet compiler TODO: Dual tree

    std::vector<result_type> results(charges.size());
    executor_->execute(charges, results);

    // TODO: don't return this, provide accessor
    return results;
  }

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
};
