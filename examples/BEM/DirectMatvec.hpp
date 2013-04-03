#pragma once

/**
 * Simple wrapper class to use direct evaluation in new GMRES implementation
 */

#include "Direct.hpp"

template <class Kernel>
class DirectMV
{
 public:
  typedef typename Kernel::charge_type charge_type;
  typedef typename Kernel::source_type source_type;
  typedef typename Kernel::source_type target_type;
  typedef typename Kernel::result_type result_type;

 private:
  const Kernel& K;
  // references to sources, targets & results
  const std::vector<source_type>& sources;
  const std::vector<target_type>& targets;

 public:

  DirectMV(const Kernel& k, const std::vector<source_type>& s, const std::vector<target_type>& t)
    : K(k), sources(s), targets(t) {};

  void set_p(int p)
  {
    (void) p;
  }

  std::vector<result_type> execute(std::vector<charge_type>& charges, unsigned dummy)
  {
    (void) dummy;
    std::vector<result_type> results(targets.size(),0);
    Direct::matvec(K, sources.begin(), sources.end(), charges.begin(),
                   targets.begin(), targets.end(), results.begin());
    return results;
  }

};
