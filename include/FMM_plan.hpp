#pragma once

// system includes
#include <algorithm>

// FMM includes

// #include <Evaluator.hpp>
#include <SimpleEvaluator.hpp>
#include <Vec.hpp>
#include "BoundingBox.hpp"
#include "Octree.hpp"


#include <Logger.hpp>
Logger Log;

typedef enum {TOPDOWN, BOTTOMUP} TreeType;
typedef enum {FMM} EvaluatorType;

/** Class to define compile-time and run-time FMM options */
class FMM_options
{
public:
  bool symmetric;
  TreeType tree;
  EvaluatorType evaluator;

  FMM_options() : symmetric(false), tree(TOPDOWN), evaluator(FMM) {};
};


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

  //private:
  FMM_options& Opts;
  SimpleEvaluator<Kernel> evaluator;
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
                         const std::vector<unsigned>& permute)
  {
    std::vector<T> temp(v.size());
    for (unsigned i=0; i < v.size(); ++i)
      temp[i] = v[permute[i]];
    return temp;
  }

  template <typename T>
  std::vector<T> ipermute(const std::vector<T>& v,
                          const std::vector<unsigned>& permute)
  {
    std::vector<T> temp(v.size());
    for (unsigned i = 0; i < v.size(); ++i)
      temp[permute[i]] = v[i];
    return temp;
  }

public:

  // CONSTRUCTOR

  FMM_plan(Kernel& k, const std::vector<point_type>& points, FMM_options& opts)
      : Opts(opts), evaluator(k),
        otree(get_boundingbox(points.begin(), points.end()))
  {
    // Construct the Octree
    otree.construct_tree(points.begin(),points.end());
  }

  // EXECUTE

  std::vector<result_type> execute(const std::vector<charge_type>& charges,
                                   const std::vector<point_type>& t_points)
  {
    // sort charges to match sorted body array
    auto pcharges = permute(charges, otree.getPermutation());

    // run upward sweep based on (permuted) body charges
    evaluator.upward(otree, pcharges);

    // run evaluator and traverse tree
    printf("Executing...\n");

    (void) t_points; // Quiet compiler TODO
    std::vector<result_type> results(charges.size());
    evaluator.downward(otree, pcharges, results);

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
