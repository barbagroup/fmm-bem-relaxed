#pragma once

/**
 * FMM evaluator based on dual tree traversal
 */

#include <EvaluatorBase.hpp>

//! forward definition
template <typename Tree, typename Kernel>
class EvaluatorBase;

template <typename Tree, typename  Kernel>
class EvaluatorFMM : public EvaluatorBase<Tree,Kernel>
{
 private:
  typedef typename EvaluatorBase<Tree,Kernel>::charge_type charge_type;
  typedef typename EvaluatorBase<Tree,Kernel>::result_type result_type;

 public:
  //! constructor
  EvaluatorFMM(Tree& t, Kernel& k)
        : EvaluatorBase<Tree,Kernel>(t,k) {};
//         : EvaluatorBase<Tree,Kernel>::tree(t), EvaluatorBase<Tree,Kernel>::K(k) {};

  //! upward sweep
  void upward(const std::vector<charge_type>& charges)
  {
    (void) charges;
  }

  //! Box-Box interactions
  void interactions(std::vector<result_type>& results)
  {
    (void) results;
  }

  //! downward sweep
  void downward(std::vector<result_type>& results)
  {
    (void) results;
  }

  //! evaluator name
  std::string name()
  {
    return "FMM";
  }
};

