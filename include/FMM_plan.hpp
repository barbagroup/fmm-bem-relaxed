#pragma once

// system includes
#include <algorithm>

// FMM includes
#include <Types.hpp>
#include <Logger.hpp>
#include <Sorter.hpp>
#include <TreeStructure.hpp>
#include <Evaluator.hpp>
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
struct fmm_wrapper
{
  virtual ~fmm_wrapper() {}
  virtual void execute(Bodies& jbodies) = 0;
};



template <class Kernel>
class FMM_plan : public fmm_wrapper
{
 public:
  // TODO: Use this as the base vector type
  typedef Vec<typename Kernel::point_type, Kernel::dimension> point_type;

private:
  Cells cells, jcells;
  Bodies buffer;
  Kernel &K;
  FMM_options& Opts;
  Evaluator<Kernel> evaluator;
  TreeStructure<point_type> tree;
  Octree<point_type> otree;

  BoundingBox<point_type> get_boundingbox(Bodies& bodies) {
    BoundingBox<point_type> result;
    for (B_iter B = bodies.begin(); B != bodies.end(); ++B) {
      result |= point_type(B->X);
    }
    std::cout << result << "\n";
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

public:
  FMM_plan(Kernel& k, Bodies& bodies, FMM_options& opts)
      : K(k), Opts(opts), evaluator(K), otree(get_boundingbox(bodies))
  {
    // do all kernel precomputation
    K.preCalculation();

    // Construct the Octree
    // ... Need point iterator

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

    // initialise evaluator
    evaluator.upward(cells);
  }

  ~FMM_plan()
  {
    // all kernel cleanup
    K.postCalculation();
  }

  void execute(Bodies& jbodies)
  {
    // run evaluator and traverse tree
    jcells = cells;
    printf("Executing...\n");
    evaluator.downward(cells,jcells,false);
  }
};


class fmm_plan
{
  fmm_wrapper* plan;
 public:
  fmm_plan() : plan(NULL) {}

  template <typename Kernel>
  fmm_plan(Kernel K, Bodies& b, FMM_options& opts) {
    plan = new FMM_plan<Kernel>(K, b, opts);
  }

  ~fmm_plan() {
    delete plan;
  }

  friend void fmm_execute(fmm_plan& p, Bodies& jbodies) {
    if (p.plan)
      p.plan->execute(jbodies);
  }
};
