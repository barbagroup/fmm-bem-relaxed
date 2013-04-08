#pragma once

#include <vector>

/** @class ExecutorBase
 * @brief Abstract Executor class
 */
template <typename Kernel>
struct ExecutorBase {
  //! Kernel type
  typedef Kernel  kernel_type;
  //! Kernel charge type
  typedef typename kernel_type::charge_type charge_type;
  //! Kernel result type
  typedef typename kernel_type::result_type result_type;

  /** Virtual Destructor */
  virtual ~ExecutorBase() {};

  // TODO: improve?
  /** Virtual execute */
  virtual void execute(const std::vector<charge_type>& charges,
                       std::vector<result_type>& results) = 0;
};
