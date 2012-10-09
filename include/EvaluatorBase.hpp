#pragma once

/**
 * Base class for all evaluators
 */

#include <FMM_plan.hpp>
#include <EvaluatorFMM.hpp>

template <class Tree, class Kernel>
class EvaluatorBase
{
 public:
  //! Kernel type
  typedef Kernel kernel_type;
  //! Point type
  typedef typename Kernel::point_type point_type;
  //! Multipole expansion type
  typedef typename Kernel::multipole_type multipole_type;
  //! Local expansion type
  typedef typename Kernel::local_type local_type;
  //! Kernel source type
  typedef typename Kernel::charge_type charge_type;
  //! Kernel result type
  typedef typename Kernel::result_type result_type;

 protected:
  //! Kernel
  Kernel& K;
  //! Octree
  Tree& tree;

  //! Multipole expansions corresponding to Box indices in Octree
  std::vector<multipole_type> M;
  //! Local expansions corresponding to Box indices in Octree
  std::vector<local_type> L;

 private:
  typename std::vector<result_type>::iterator results_begin;
  typename std::vector<charge_type>::const_iterator charges_begin;

 public:
  EvaluatorBase(Tree& t, Kernel& k) : tree(t), K(k) {};

  //! Upward sweep
  virtual void upward(const std::vector<charge_type>& charges) = 0;
  //! 'Interaction' stage
  virtual void interactions(std::vector<result_type>& results) = 0;
  //! Downward sweep
  virtual void downward(std::vector<result_type>& results) = 0;

  //! Abstract factory
  static EvaluatorBase<Tree,Kernel> *createEvaluator(Tree& t, Kernel& k, FMM_options& options)
  {
    if (options.evaluator == FMM)
    {
      return new EvaluatorFMM<Tree,Kernel>(t,k);
    }
    else
    {
      return NULL;
    }
  }

  //! name of the evaluator
  virtual std::string name() = 0;
};

