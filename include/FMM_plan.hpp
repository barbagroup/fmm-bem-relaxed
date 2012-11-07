#pragma once

// FMM includes
#include <FMMOptions.hpp>
#include <Vec.hpp>
#include "BoundingBox.hpp"

#include <Logger.hpp>

#include <tree/Octree.hpp>
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

  typedef Octree<source_type, point_type> tree_type;
  typedef Kernel kernel_type;

  //private:
  ExecutorBase<tree_type,kernel_type>* executor;
  Kernel& K;
  Octree<source_type, point_type> otree;

  template <typename SourceIter>
  BoundingBox<point_type> get_boundingbox(SourceIter begin, SourceIter end) {
    BoundingBox<point_type> result;
    for ( ; begin != end; ++begin)
    result |= static_cast<point_type>(*begin);
    // Make sure the bounding box is square and slightly scaled
    // TODO: improve
    auto dim = result.dimensions();
    auto maxdim = std::max(dim[0], std::max(dim[1], dim[2]));
    result |= result.min() + point_type(maxdim, maxdim, maxdim) * (1 + 1e-6);
    //std::cout << "Bounding Box: " << result << "\n";
    return result;
  }

public:

  // CONSTRUCTOR

  FMM_plan(Kernel& k, const std::vector<source_type>& points,
           FMMOptions& opts)
    : K(k), //evaluator(k),
      otree(get_boundingbox(points.begin(), points.end())) {
    // Construct the Octree
    otree.construct_tree(points.begin(), points.end(), opts.NCRIT);
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
                                   const std::vector<target_type>& t_points)
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
