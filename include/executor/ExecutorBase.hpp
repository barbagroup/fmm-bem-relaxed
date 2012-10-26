#pragma once

/** @class ExecutorBase
 * @brief Abstract Executor class
 */
template <typename Tree, typename Kernel>
struct ExecutorBase {
  ExecutorBase() {}
  virtual ~ExecutorBase() {};
  // TODO: improve
  virtual void execute(const std::vector<typename Kernel::charge_type>& charges,
                       std::vector<typename Kernel::result_type>& results) = 0;
};
