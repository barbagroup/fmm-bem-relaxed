#pragma once

// system includes
#include <algorithm>

// FMM includes
#include <Types.hpp>
#include <Logger.hpp>
#include <Sorter.hpp>
#include <TreeStructure.hpp>
// #include <Evaluator.hpp>
#include <SimpleEvaluator.hpp>
#include <Vec.hpp>
#include "BoundingBox.hpp"
#include "Octree.hpp"

Logger Log;
Sorter sort;

typedef enum {TOPDOWN, BOTTOMUP} treeType;

/** Class to define compile-time and run-time FMM options */
class FMM_options
{
public:
  bool symmetric;
  treeType tree;

  FMM_options() : symmetric(false), tree(TOPDOWN) {};
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
  //typedef typename Kernel::point_type point_type;
  // TODO: Use this as the base vector type?
  typedef typename Kernel::point_type point_type;
  typedef typename Kernel::charge_type charge_type;
  typedef typename Kernel::result_type result_type;


  //private:
  Cells cells, jcells;
  Bodies buffer;
  Kernel &K;
  FMM_options& Opts;
  SimpleEvaluator<Kernel> evaluator;
  TreeStructure<point_type> tree;
  Octree<point_type> otree;
  std::vector<point_type> source_points;
  std::vector<charge_type> charges;
  std::vector<result_type> results;

  BoundingBox<point_type> get_boundingbox(Bodies& bodies) {
    BoundingBox<point_type> result;
    for (B_iter B = bodies.begin(); B != bodies.end(); ++B) {
      result |= point_type(B->X[0], B->X[1], B->X[2]);
    }
    // Make sure the bounding box is square for now TODO: improve
    auto dim = result.dimensions();
    auto maxdim = std::max(dim[0], std::max(dim[1], dim[2]));
    result |= result.min() + point_type(maxdim, maxdim, maxdim);
    std::cout << "Bounding Box: " << result << "\n";
    return result;
  }

  template <typename PointIter>
  BoundingBox<point_type> get_boundingbox(PointIter begin, PointIter end) {
    BoundingBox<point_type> result;
    for ( ; begin != end; ++begin)
      result |= *begin;
    // Make sure the bounding box is square for now TODO: improve
    auto dim = result.dimensions();
    auto maxdim = std::max(dim[0], std::max(dim[1], dim[2]));
    result |= result.min() + point_type(maxdim, maxdim, maxdim);
    std::cout << "Bounding Box: " << result << "\n";
    return result;
  }

  void setDomain(Bodies& bodies, point_type& X0, double& R0) {
    vect xmin,xmax;                                             // Min,Max of domain
    B_iter B = bodies.begin();                                  // Reset body iterator
    xmin = xmax = B->X;                                         // Initialize xmin,xmax
    for( B=bodies.begin(); B!=bodies.end(); ++B ) {             // Loop over bodies
      for( int d=0; d!=3; ++d ) {                               //  Loop over each dimension
        if     (B->X[d] < xmin[d]) xmin[d] = B->X[d];           //   Determine xmin
        else if(B->X[d] > xmax[d]) xmax[d] = B->X[d];           //   Determine xmax
      }                                                         //  End loop over each dimension
    }                                                           // End loop over bodies

    for( int d=0; d!=3; ++d ) {                               // Loop over each dimension
      X0[d] = (xmax[d] + xmin[d]) / 2;                        // Calculate center of domain
      X0[d] = int(X0[d]+.5);                                  //  Shift center to nearest integer
      R0 = std::max(xmax[d] - X0[d], R0);                     //  Calculate max distance from center
      R0 = std::max(X0[d] - xmin[d], R0);                     //  Calculate max distance from center
    }                                                         // End loop over each dimension
    R0 *= 1.000001;                                           // Add some leeway to root radius
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

  void bodies2points(Bodies& bodies, std::vector<point_type>& points)
  {
    points.resize(bodies.size());

    int idx = 0;
    for (auto it=bodies.begin(); it!=bodies.end(); ++it)
    {
      point_type p;
      p[0] = it->X[0];
      p[1] = it->X[1];
      p[2] = it->X[2];
      points[idx] = p;
      idx++;
    }
  }

  FMM_plan(Kernel& k, const std::vector<point_type>& points, FMM_options& opts)
      : K(k), Opts(opts), evaluator(K),
        otree(get_boundingbox(points.begin(), points.end()))
  {
    // do all kernel precomputation
    K.preCalculation();

    // Construct the Octree
    otree.construct_tree(points.begin(),points.end());
  }

  FMM_plan(Kernel& k, Bodies& bodies, FMM_options& opts)
      : K(k), Opts(opts), evaluator(K), otree(get_boundingbox(bodies))
  {
    // do all kernel precomputation
    K.preCalculation();

    // Construct the Octree
    // ... Need point iterator
    // create copy of bodies into array of points
    bodies2points(bodies,source_points);
    otree.construct_tree(source_points.begin(),source_points.end());

#if 0
    // initialise tree & construct
    point_type X0;
    double R0;
    setDomain(bodies, X0, R0);
    tree.init(X0, R0);
    if (opts.tree == TOPDOWN)
      tree.topdown(bodies,cells);
    else
      tree.bottomup(bodies,cells);
    printf("Tree created: %d cells\n",(int)cells.size());
#endif
  }

  ~FMM_plan()
  {
    // all kernel cleanup
    K.postCalculation();
  }

  std::vector<result_type> execute(const std::vector<charge_type>& charges,
                                   const std::vector<point_type>& jbodies)
  {
    // sort charges to match sorted body array
    auto pcharges = permute(charges, otree.getPermutation());

    // run upward sweep based on body charges
    evaluator.upward(otree, charges);

    // run evaluator and traverse tree
    jcells = cells;
    printf("Executing...\n");
    // evaluator.downward(cells,jcells,false);
    std::vector<result_type> results;
    results.resize(charges.size());
    evaluator.downward(otree, charges, results);

    auto ipresults = ipermute(results,otree.getPermutation());
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
