#pragma once

// system includes
#include <algorithm>

// FMM includes
#include <FMMOptions.hpp>
#include <Vec.hpp>
#include "BoundingBox.hpp"
#include "Octree.hpp"
#include <Logger.hpp>

#include "Executor.hpp"

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
  typedef typename Kernel::point_type point_type;
  // TODO: Better point support?
  // Want all derived classes to use the following fmmplan::point_type wrapper?
  //typedef typename Vec<Kernel::dimension, typename Kernel::point_type> point_type;
  typedef typename Kernel::charge_type charge_type;
  typedef typename Kernel::result_type result_type;

  typedef Octree<point_type> tree_type;
  typedef Kernel kernel_type;

  //private:
  FMMOptions& options;
  ExecutorBase<tree_type,kernel_type>* executor;
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

public:

  // CONSTRUCTOR

  FMM_plan(Kernel& k, const std::vector<point_type>& points,
           FMMOptions& opts)
    : options(opts), K(k), //evaluator(k),
      otree(get_boundingbox(points.begin(), points.end())) {
    // Construct the Octree
    otree.construct_tree(points.begin(),points.end(),opts);
    // setup the executor
    set_options(opts);
  }

  // DESTRUCTOR

  ~FMM_plan() {
    delete executor;
  }

  // EXECUTE

  /** Set the executor strategy of this plan at runtime
   */
  void set_options(FMMOptions& options) {
    executor = make_executor(otree, K, options);
  }

  /** Execute this FMM plan
   */
  std::vector<result_type> execute(const std::vector<charge_type>& charges,
                                   const std::vector<point_type>& t_points)
  {
    if (!executor) {
      printf("[E]: Executor not initialised -- returning..\n");
      return std::vector<result_type>(0);
    }

    (void) t_points; // Quiet compiler TODO

    // sort charges to match sorted body array
    // TODO: not here
    auto pcharges = otree.permute(charges);

    std::vector<result_type> results(pcharges.size());
    executor->execute(pcharges, results);

    // inverse permute results
    // TODO: not here
    auto ipresults = otree.ipermute(results);

    // TODO: don't return this, provide accessor
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
