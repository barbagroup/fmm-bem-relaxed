#pragma once

// system includes
#include <algorithm>

// FMM includes

// #include <Evaluator.hpp>
#include <FMMOptions.hpp>
// #include <SimpleEvaluator.hpp>
#include <Vec.hpp>
#include "BoundingBox.hpp"
#include "Octree.hpp"
#include <Logger.hpp>

//#include <EvaluatorBase.hpp>
#include <Evaluator.hpp>


//! global logging
Logger Log;

// forward declarations
//template <class Tree, class Kernel>
//class EvaluatorBase;


/** Simple wrapper class for FMM_plan */
/*
struct fmm_wrapper
{
  virtual ~fmm_wrapper() {}

  template <typename charge_type>
  virtual void execute(std::vector<charge_type>& charges, Bodies& jbodies) = 0;
};
*/


template <class Kernel>
class FMM_plan//  : public fmm_wrapper
{
 public:
  typedef typename Kernel::point_type point_type;  // TODO: Better point support
  typedef typename Kernel::charge_type charge_type;
  typedef typename Kernel::result_type result_type;

  typedef Octree<point_type> tree_type;
  typedef Kernel kernel_type;

  //private:
  FMMOptions& options;
  // SimpleEvaluator<Kernel> evaluator;
  ExecutorBase<tree_type,kernel_type>* evaluator;
  Kernel& K;
  Octree<point_type> otree;


  template <typename PointIter>
  BoundingBox<point_type> get_boundingbox(PointIter begin, PointIter end) {
    BoundingBox<point_type> result;
    for ( ; begin != end; ++begin)
      result |= *begin;
    // Make sure the bounding box is square and slightly scaled
    // TODO: improve
    auto dim = result.dimensions();
    auto maxdim = std::max(dim[0], std::max(dim[1], dim[2]));
    result |= result.min() + point_type(maxdim, maxdim, maxdim) * (1 + 1e-6);
    std::cout << "Bounding Box: " << result << "\n";
    return result;
  }

  template <typename T>
  std::vector<T> permute(const std::vector<T>& v,
                         const std::vector<unsigned>& permute) {
    std::vector<T> temp(v.size());
    for (unsigned i=0; i < v.size(); ++i)
      temp[i] = v[permute[i]];
    return temp;
  }

  template <typename T>
  std::vector<T> ipermute(const std::vector<T>& v,
                          const std::vector<unsigned>& permute) {
    std::vector<T> temp(v.size());
    for (unsigned i = 0; i < v.size(); ++i)
      temp[permute[i]] = v[i];
    return temp;
  }

  //! Set the evaluator strategy of this plan at runtime
  void set_evaluator(FMMOptions::EvaluatorType type) {
    (void) type;
    if (type == FMMOptions::FMM) {
      auto eval = make_evaluator(new EvalUpward<tree_type,kernel_type,FMMOptions>(otree,K,options),
                                 new EvalInteraction<tree_type,kernel_type,FMMOptions>(otree,K,options),
                                 new EvalDownward<tree_type,kernel_type,FMMOptions>(otree,K,options));
      evaluator = make_executor(otree, K, eval);
    } else {

    }
    /*
    if (type == FMMOptions::FMM) {
      evaluator = new EvaluatorFMM<tree_type,kernel_type>(otree,K,options.THETA);
    } else if (type == FMMOptions::TREECODE) {
      evaluator = new EvaluatorTreecode<tree_type,kernel_type>(otree,K,options.THETA);
    } else {
      evaluator = NULL;
    }
    */
  }

public:

  // CONSTRUCTOR

  FMM_plan(Kernel& k, const std::vector<point_type>& points,
           FMMOptions& opts)
      : options(opts), K(k), //evaluator(k),
        otree(get_boundingbox(points.begin(), points.end())) {
    // Construct the Octree
    otree.construct_tree(points.begin(),points.end());
  }

  // DESTRUCTOR
  ~FMM_plan() {
    if (evaluator) delete evaluator;
  }

  // EXECUTE

  std::vector<result_type> execute(const std::vector<charge_type>& charges,
                                   const std::vector<point_type>& t_points)
  {
    // setup the evaluator
    set_evaluator(options.evaluator);

    if (!evaluator) {
      printf("[E]: Evaluator not initialised -- returning..\n");
      return std::vector<result_type>(0);
    }

    (void) t_points; // Quiet compiler TODO

    // sort charges to match sorted body array
    auto pcharges = permute(charges, otree.getPermutation());
    std::vector<result_type> results(charges.size());

    evaluator->execute(pcharges, results);


    // run upward sweep based on (permuted) body charges
    //evaluator->upward(pcharges);

    // run evaluator and traverse tree
    //printf("Executing...\n");

    //evaluator->interactions(results);
    //evaluator->downward(results);

    // inverse permute results
    auto ipresults = ipermute(results, otree.getPermutation());
    // TODO, don't return this
    return ipresults;
  }
};

/*
class fmm_plan
{
  fmm_wrapper* plan;
 public:
  // typedef typename Kernel::charge_type charge_type;
  fmm_plan() : plan(NULL) {}

  template <typename Kernel>
  fmm_plan(Kernel K, Bodies& b, FMM_options& opts) {
    plan = new FMM_plan<Kernel>(K, b, opts);
  }

  ~fmm_plan() {
    delete plan;
  }

  template <typename charge_type>
  friend void fmm_execute(fmm_plan& p, std::vector<charge_type>& charges, Bodies& jbodies) {
    if (p.plan)
      p.plan->execute(charges, jbodies);
  }
};
*/
