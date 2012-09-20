#pragma once

// system includes
#include <algorithm>

// FMM includes
#include <Types.hpp>
#include <Logger.hpp>
#include <Sorter.hpp>
#include <TreeStructure.hpp>
#include <Evaluator.hpp>

Logger Log;
Sorter sort;

typedef enum {TOPDOWN, BOTTOMUP} treeType;

class FMM_options
{
public:
  bool symmetric;
  treeType tree;

  FMM_options() : symmetric(false), tree(TOPDOWN) {};
};

template <class Kernel>
class FMM_plan
{
private:
  Cells cells, jcells;
  vect X0;
  real R0;
  Bodies buffer;
  TreeStructure tree;
  Evaluator<Kernel> *eval;
  Kernel &K;

  void setDomain(Bodies& bodies, vect x0=0, real r0=M_PI)
  {
    vect xmin,xmax;                                             // Min,Max of domain
    B_iter B = bodies.begin();                                  // Reset body iterator
    xmin = xmax = B->X;                                         // Initialize xmin,xmax
    for( B=bodies.begin(); B!=bodies.end(); ++B ) {             // Loop over bodies
      for( int d=0; d!=3; ++d ) {                               //  Loop over each dimension
        if     (B->X[d] < xmin[d]) xmin[d] = B->X[d];           //   Determine xmin
        else if(B->X[d] > xmax[d]) xmax[d] = B->X[d];           //   Determine xmax
      }                                                         //  End loop over each dimension
    }                                                           // End loop over bodies
    if( IMAGES != 0 ) {                                         // If periodic boundary condition
      if( xmin[0] < x0[0]-r0 || x0[0]+r0 < xmax[0]              //  Check for outliers in x direction
       || xmin[1] < x0[1]-r0 || x0[1]+r0 < xmax[1]              //  Check for outliers in y direction
       || xmin[2] < x0[2]-r0 || x0[2]+r0 < xmax[2] ) {          //  Check for outliers in z direction
         std::cout << "Error: Particles located outside periodic domain : " << std::endl;// Print error message
         std::cout << xmin << std::endl;
         std::cout << xmax << std::endl;
      }                                                         //  End if for outlier checking
      X0 = x0;                                                  //  Center is [0, 0, 0]
      R0 = r0;                                                  //  Radius is r0
    } else {
      for( int d=0; d!=3; ++d ) {                               // Loop over each dimension
        X0[d] = (xmax[d] + xmin[d]) / 2;                        // Calculate center of domain
        X0[d] = int(X0[d]+.5);                                  //  Shift center to nearest integer
        R0 = std::max(xmax[d] - X0[d], R0);                     //  Calculate max distance from center
        R0 = std::max(X0[d] - xmin[d], R0);                     //  Calculate max distance from center
      }                                                         // End loop over each dimension
      R0 *= 1.000001;                                           // Add some leeway to root radius
    }                                                           // Endif for periodic boundary condition
  }

public:
  FMM_plan(Kernel& k, Bodies& bodies, FMM_options& opts) : K(k)
  {
    // set domain of problem (center & radius)
    setDomain(bodies);

    // do all kernel precomputation
    K.preCalculation();

    // initialise tree & construct
    tree.init(X0,R0);
    tree.topdown(bodies,cells,K);
    printf("Tree created: %d cells\n",(int)cells.size());

    // initialise evaluator
    eval = new Evaluator<Kernel>(K,R0);
    eval->upward(cells);
  }

  ~FMM_plan()
  {
    // all kernel cleanup
    K.postCalculation();
    delete eval;
  }

  void execute(Bodies& jbodies)
  {
    // run evaluator and traverse tree
    jcells = cells;
    printf("executing...\n");
    eval->downward(cells,jcells,false);
  }

  void checkError(Bodies& FMM_bodies, Bodies& sources)
  {
    int numTargets = 100;
    FMM_bodies.resize(numTargets);
    Bodies test_bodies = FMM_bodies;

    for (B_iter B=test_bodies.begin(); B!=test_bodies.end(); ++B)
    {
      B->IBODY = B-test_bodies.begin();
      B->TRG = 0;
    }

    eval->evalP2P(test_bodies,sources);

    B_iter B2 = FMM_bodies.begin();

    real diff1=0, diff2=0, norm1=0, norm2=0;
    for (B_iter B=test_bodies.begin(); B!=test_bodies.end(); ++B, ++B2)
    {
      diff1 += (B->TRG[0] - B2->TRG[0]) * (B->TRG[0] - B2->TRG[0]);
      norm1 += B2->TRG[0] * B2->TRG[0];

      diff2 += (B->TRG[1] - B2->TRG[1]) * (B->TRG[1] - B2->TRG[1]);
      diff2 += (B->TRG[2] - B2->TRG[2]) * (B->TRG[2] - B2->TRG[2]);
      diff2 += (B->TRG[3] - B2->TRG[3]) * (B->TRG[3] - B2->TRG[3]);
      norm2 += B2->TRG[1] * B2->TRG[1];
      norm2 += B2->TRG[2] * B2->TRG[2];
      norm2 += B2->TRG[3] * B2->TRG[3];
    }

    printf("Error (pot) : %.4e\n",sqrt(diff1/norm1));
    printf("Error (acc) : %.4e\n",sqrt(diff2/norm2));
  }
};

