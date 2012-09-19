/*
Copyright (C) 2011 by Rio Yokota, Simon Layton, Lorena Barba

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in
all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
THE SOFTWARE.
*/

#pragma once

#include <Types.hpp>

#define splitFirst(Ci,Cj) Cj->NCHILD == 0 || (Ci->NCHILD != 0 && Ci->R > Cj->R)

//! Interface between tree and kernel
template <class Kernel>
class Evaluator
{
private:
  real        timeM2L;                                          //!< M2L execution time
  real        timeM2P;                                          //!< M2P execution time
  real        timeP2P;                                          //!< P2P execution time
  C_iter      Ci0, Cj0;
  real        R0;

protected:
  C_iter      CiB;                                              //!< icells begin per call
  C_iter      CiE;                                              //!< icells end per call
  Lists       listM2L;                                          //!< M2L interaction list
  Lists       listM2P;                                          //!< M2P interaction list
  Lists       listP2P;                                          //!< P2P interaction list

  int         Iperiodic;                                        //!< Periodic image flag (using each bit for images)
  int         Icenter;                                          //!< Periodic image flag at center
  Maps        flagM2L;                                          //!< Existance of periodic image for M2L
  Maps        flagM2P;                                          //!< Existance of periodic image for M2P
  Maps        flagP2P;                                          //!< Existance of periodic image for P2P

  real        NP2P;                                             //!< Number of P2P kernel calls
  real        NM2P;                                             //!< Number of M2P kernel calls
  real        NM2L;                                             //!< Number of M2L kernel calls

  Kernel &K;

private:
  //! Approximate interaction between two cells
  inline void approximate(C_iter Ci, C_iter Cj) {
#if HYBRID
    if( timeP2P*Cj->NDLEAF < timeM2P && timeP2P*Ci->NDLEAF*Cj->NDLEAF < timeM2L) {// If P2P is fastest
      evalP2P(Ci,Cj);                                           //  Evaluate on CPU, queue on GPU
    } else if ( timeM2P < timeP2P*Cj->NDLEAF && timeM2P*Ci->NDLEAF < timeM2L ) {// If M2P is fastest
      evalM2P(Ci,Cj);                                           //  Evaluate on CPU, queue on GPU
    } else {                                                    // If M2L is fastest
      evalM2L(Ci,Cj);                                           //  Evaluate on CPU, queue on GPU
    }                                                           // End if for kernel selection
#elif TREECODE
    evalM2P(Ci,Cj);                                             // Evaluate on CPU, queue on GPU
#else
    evalM2L(Ci,Cj);                                             // Evalaute on CPU, queue on GPU
#endif
  }

  //! Traverse a pair of trees using a queue
  void traverseQueue(Pair pair) {
    PairQueue pairQueue;                                        // Queue of interacting cell pairs
    pairQueue.push_back(pair);                                  // Push pair to queue
    while( !pairQueue.empty() ) {                               // While dual traversal queue is not empty
      pair = pairQueue.front();                                 //  Get interaction pair from front of queue
      pairQueue.pop_front();                                    //  Pop dual traversal queue
      if(splitFirst(pair.first,pair.second)) {                  //  If first cell is larger
        C_iter C = pair.first;                                  //   Split the first cell
        for( C_iter Ci=Ci0+C->CHILD; Ci!=Ci0+C->CHILD+C->NCHILD; ++Ci ) {// Loop over first cell's children
          interact(Ci,pair.second,pairQueue);                   //    Calculate interaction between cells
        }                                                       //   End loop over fist cell's children
      } else {                                                  //  Else if second cell is larger
        C_iter C = pair.second;                                 //   Split the second cell
        for( C_iter Cj=Cj0+C->CHILD; Cj!=Cj0+C->CHILD+C->NCHILD; ++Cj ) {// Loop over second cell's children
          interact(pair.first,Cj,pairQueue);                    //    Calculate interaction betwen cells
        }                                                       //   End loop over second cell's children
      }                                                         //  End if for which cell to split
    }                                                           // End while loop for dual traversal queue
  }

  //! Get range of periodic images
  int getPeriodicRange() {
    int prange = 0;                                             //  Range of periodic images
    for( int i=0; i!=IMAGES; ++i ) {                            //  Loop over periodic image sublevels
      prange += int(pow(3,i));                                  //   Accumulate range of periodic images
    }                                                           //  End loop over perioidc image sublevels
    return prange;                                              // Return range of periodic images
  }

protected:
  //! Get level from cell index
  int getLevel(bigint index) {
    int i = index;                                              // Copy to dummy index
    int level = -1;                                             // Initialize level counter
    while( i >= 0 ) {                                           // While cell index is non-negative
      level++;                                                  //  Increment level
      i -= 1 << 3*level;                                        //  Subtract number of cells in that level
    }                                                           // End while loop for cell index
    return level;                                               // Return the level
  }

  void timeKernels();                                           //!< Time all kernels for auto-tuning

  //! Upward phase for periodic cells
  void upwardPeriodic(Cells &jcells) {
    Cells pccells, pjcells;                                     // Periodic center cell and jcell
    pccells.push_back(jcells.back());                           // Root cell is first periodic cell
    for( int level=0; level<IMAGES-1; ++level ) {               // Loop over sublevels of tree
      Cell cell;                                                //  New periodic cell at next sublevel
      C_iter C = pccells.end() - 1;                             //  Set previous periodic center cell as source
      for( int ix=-1; ix<=1; ++ix ) {                           //  Loop over x periodic direction
        for( int iy=-1; iy<=1; ++iy ) {                         //   Loop over y periodic direction
          for( int iz=-1; iz<=1; ++iz ) {                       //    Loop over z periodic direction
            if( ix != 0 || iy != 0 || iz != 0 ) {               //     If periodic cell is not at center
              for( int cx=-1; cx<=1; ++cx ) {                   //      Loop over x periodic direction (child)
                for( int cy=-1; cy<=1; ++cy ) {                 //       Loop over y periodic direction (child)
                  for( int cz=-1; cz<=1; ++cz ) {               //        Loop over z periodic direction (child)
                    cell.X[0]  = C->X[0] + (ix * 6 + cx * 2) * C->R;//     Set new x coordinate for periodic image
                    cell.X[1]  = C->X[1] + (iy * 6 + cy * 2) * C->R;//     Set new y cooridnate for periodic image
                    cell.X[2]  = C->X[2] + (iz * 6 + cz * 2) * C->R;//     Set new z coordinate for periodic image
                    cell.M     = C->M;                          //         Copy multipoles to new periodic image
                    cell.NCLEAF = cell.NDLEAF = cell.NCHILD = 0;//         Initialize NCLEAF, NDLEAF, & NCHILD
                    jcells.push_back(cell);                     //         Push cell into periodic jcell vector
                  }                                             //        End loop over z periodic direction (child)
                }                                               //       End loop over y periodic direction (child)
              }                                                 //      End loop over x periodic direction (child)
            }                                                   //     Endif for periodic center cell
          }                                                     //    End loop over z periodic direction
        }                                                       //   End loop over y periodic direction
      }                                                         //  End loop over x periodic direction
      for( int ix=-1; ix<=1; ++ix ) {                           //  Loop over x periodic direction
        for( int iy=-1; iy<=1; ++iy ) {                         //   Loop over y periodic direction
          for( int iz=-1; iz<=1; ++iz ) {                       //    Loop over z periodic direction
            cell.X[0] = C->X[0] + ix * 2 * C->R;                //     Set new x coordinate for periodic image
            cell.X[1] = C->X[1] + iy * 2 * C->R;                //     Set new y cooridnate for periodic image
            cell.X[2] = C->X[2] + iz * 2 * C->R;                //     Set new z coordinate for periodic image
            cell.M = C->M;                                      //     Copy multipoles to new periodic image
            pjcells.push_back(cell);                            //     Push cell into periodic jcell vector
          }                                                     //    End loop over z periodic direction
        }                                                       //   End loop over y periodic direction
      }                                                         //  End loop over x periodic direction
      cell.X = C->X;                                            //  This is the center cell
      cell.R = 3 * C->R;                                        //  The cell size increases three times
      pccells.pop_back();                                       //  Pop periodic center cell from vector
      pccells.push_back(cell);                                  //  Push cell into periodic cell vector
      C_iter Ci = pccells.end() - 1;                            //  Set current cell as target for M2M
      Ci->CHILD = 0;                                            //  Set child cells for periodic M2M
      Ci->NCHILD = 27;                                          //  Set number of child cells for periodic M2M
      evalM2M(pccells,pjcells);                                 // Evaluate periodic M2M kernels for this sublevel
      pjcells.clear();                                          // Clear periodic jcell vector
    }                                                           // End loop over sublevels of tree
  }

  //! Traverse tree for periodic cells
  void traversePeriodic(Cells &cells, Cells &jcells) {
    Xperiodic = 0;                                              // Set periodic coordinate offset
    Iperiodic = Icenter;                                        // Set periodic flag to center
    C_iter Cj = jcells.end()-1;                                 // Initialize iterator for periodic source cell
    for( int level=0; level<IMAGES-1; ++level ) {               // Loop over sublevels of tree
      for( int I=0; I!=26*27; ++I, --Cj ) {                     //  Loop over periodic images (exclude center)
#if TREECODE
        for( C_iter Ci=cells.begin(); Ci!=cells.end(); ++Ci ) { //   Loop over cells
          if( Ci->NCHILD == 0 ) {                               //    If cell is twig
            evalM2P(Ci,Cj);                                     //     Perform M2P kernel
          }                                                     //    Endif for twig
        }                                                       //   End loop over cells
#else
        C_iter Ci = cells.end() - 1;                            //   Set root cell as target iterator
        evalM2L(Ci,Cj);                                         //   Perform M2P kernel
#endif
      }                                                         //  End loop over x periodic direction
    }                                                           // End loop over sublevels of tree
  }


public:
  //! Constructor
  Evaluator() : R0(0), Icenter(1 << 13), NP2P(0), NM2P(0), NM2L(0), K(Kernel()) {};
  Evaluator(Kernel &k, real r0) : R0(r0),  Icenter(1 << 13), NP2P(0), NM2P(0), NM2L(0), K(k){};
  //! Destructor
  ~Evaluator() {}

  //! Add single list for kernel unit test
  void addM2L(C_iter Cj) {
    listM2L.resize(1);                                          // Resize vector of M2L interation lists
    flagM2L.resize(1);                                          // Resize vector of M2L periodic image flags
    listM2L[0].push_back(Cj);                                   // Push single cell into list
    flagM2L[0][Cj] |= Icenter;                                  // Flip bit of periodic image flag
  }

  //! Add single list for kernel unit test
  void addM2P(C_iter Cj) {
    listM2P.resize(1);                                          // Resize vector of M2P interation lists
    flagM2P.resize(1);                                          // Resize vector of M2L periodic image flags
    listM2P[0].push_back(Cj);                                   // Push single cell into list
    flagM2P[0][Cj] |= Icenter;                                  // Flip bit of periodic image flag
  }

  //! Create periodic images of bodies
  Bodies periodicBodies(Bodies &bodies) {
    Bodies jbodies;                                             // Vector for periodic images of bodies
    int prange = getPeriodicRange();                            // Get range of periodic images
    for( int ix=-prange; ix<=prange; ++ix ) {                   // Loop over x periodic direction
      for( int iy=-prange; iy<=prange; ++iy ) {                 //  Loop over y periodic direction
        for( int iz=-prange; iz<=prange; ++iz ) {               //   Loop over z periodic direction
          for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {//    Loop over bodies
            Body body = *B;                                     //     Copy current body
            body.X[0] += ix * 2 * R0;                           //     Shift x position
            body.X[1] += iy * 2 * R0;                           //     Shift y position
            body.X[2] += iz * 2 * R0;                           //     Shift z position
            jbodies.push_back(body);                            //     Push shifted body into jbodies
          }                                                     //    End loop over bodies
        }                                                       //   End loop over z periodic direction
      }                                                         //  End loop over y periodic direction
    }                                                           // End loop over x periodic direction
    return jbodies;                                             // Return vector for periodic images of bodies
  }

  //! Use multipole acceptance criteria to determine whether to approximate, do P2P, or subdivide
  void interact(C_iter Ci, C_iter Cj, PairQueue &pairQueue) {
    vect dX = Ci->X - Cj->X - Xperiodic;                        // Distance vector from source to target
    real Rq = std::sqrt(norm(dX));                              // Scalar distance
    if( Rq * THETA > Ci->R + Cj->R ) {                          // If distance if far enough
      approximate(Ci,Cj);                                       //  Use approximate kernels, e.g. M2L, M2P
    } else if(Ci->NCHILD == 0 && Cj->NCHILD == 0) {             // Else if both cells are leafs
      evalP2P(Ci,Cj);                                           //  Use P2P
    } else {                                                    // If cells are close but not leafs
      Pair pair(Ci,Cj);                                         //  Form a pair of cell iterators
      pairQueue.push_back(pair);                                //  Push pair to queue
    }                                                           // End if for multipole acceptance
  }

  //! Dual tree traversal
  void traverse(Cells &cells, Cells &jcells) {
    C_iter root = cells.end() - 1;                              // Iterator for root target cell
    C_iter jroot = jcells.end() - 1;                            // Iterator for root source cell
    if( IMAGES != 0 ) {                                         // If periodic boundary condition
      jroot = jcells.end() - 1 - 26 * 27 * (IMAGES - 1);        //  The root is not at the end
    }                                                           // Endif for periodic boundary condition
    Ci0 = cells.begin();                                        // Set begin iterator for target cells
    Cj0 = jcells.begin();                                       // Set begin iterator for source cells
    K.Ci0 = Ci0;
    K.Cj0 = Cj0;
    if( IMAGES == 0 ) {                                         // If free boundary condition
      Iperiodic = Icenter;                                      //  Set periodic image flag to center
      Xperiodic = 0;                                            //  Set periodic coordinate offset
      Pair pair(root,jroot);                                    //  Form pair of root cells
      traverseQueue(pair);                                      //  Traverse a pair of trees
    } else {                                                    // If periodic boundary condition
      int I = 0;                                                //  Initialize index of periodic image
      for( int ix=-1; ix<=1; ++ix ) {                           //  Loop over x periodic direction
        for( int iy=-1; iy<=1; ++iy ) {                         //   Loop over y periodic direction
          for( int iz=-1; iz<=1; ++iz, ++I ) {                  //    Loop over z periodic direction
            Iperiodic = 1 << I;                                 //     Set periodic image flag
            Xperiodic[0] = ix * 2 * R0;                         //     Coordinate offset for x periodic direction
            Xperiodic[1] = iy * 2 * R0;                         //     Coordinate offset for y periodic direction
            Xperiodic[2] = iz * 2 * R0;                         //     Coordinate offset for z periodic direction
            Pair pair(root,jroot);                              //     Form pair of root cells
            traverseQueue(pair);                                //     Traverse a pair of trees
          }                                                     //    End loop over z periodic direction
        }                                                       //   End loop over y periodic direction
      }                                                         //  End loop over x periodic direction
      for( C_iter Ci=cells.begin(); Ci!=cells.end(); ++Ci ) {   //  Loop over target cells
        listM2L[Ci-Ci0].sort();                                 //  Sort interaction list
        listM2L[Ci-Ci0].unique();                               //  Eliminate duplicate periodic entries
        listM2P[Ci-Ci0].sort();                                 //  Sort interaction list
        listM2P[Ci-Ci0].unique();                               //  Eliminate duplicate periodic entries
        listP2P[Ci-Ci0].sort();                                 //  Sort interaction list
        listP2P[Ci-Ci0].unique();                               //  Eliminate duplicate periodic entries
      }                                                         //  End loop over target cells
    }                                                           // Endif for periodic boundary condition
  }

  //! Downward phase (M2L,M2P,P2P,L2L,L2P evaluation)
  void downward(Cells &cells, Cells &jcells, bool periodic=true) {
#if HYBRID
    timeKernels();                                              // Time all kernels for auto-tuning
#endif
    for( C_iter C=cells.begin(); C!=cells.end(); ++C ) {        // Initialize local coefficients
      for (size_t i=0; i<C->L.size(); i++) C->L[i] = 0;
    }
    if( IMAGES != 0 ) {                                         // If periodic boundary condition
      Log.startTimer("Upward P");                                   //  Start timer
      upwardPeriodic(jcells);                                   //  Upward phase for periodic images
      Log.stopTimer("Upward P");                           //  Stop timer & print
    }                                                           // Endif for periodic boundary condition
    Log.startTimer("Traverse");                                     // Start timer
    traverse(cells,jcells);                                     // Traverse tree to get interaction list
    Log.stopTimer("Traverse");                             // Stop timer & print
    if( IMAGES != 0 && periodic ) {                             // If periodic boundary condition
      Log.startTimer("Traverse P");                                 // Start timer
      traversePeriodic(cells,jcells);                           // Traverse tree for periodic images
      Log.stopTimer("Traverse P");                         // Stop timer & print
    }                                                           // Endif for periodic boundary condition
    evalL2L(cells);                                             // Evaluate all L2L kernels
    evalL2P(cells);                                             // Evaluate all L2P kernels
    if(true) std::cout << "P2P: "  << NP2P
                           << " M2P: " << NM2P
                           << " M2L: " << NM2L << std::endl;
  }

  void setSourceBody();                                         //!< Set source buffer for bodies (for GPU)
  void setSourceCell(bool isM);                                 //!< Set source buffer for cells (for GPU)
  void setTargetBody(Lists lists, Maps flags);                  //!< Set target buffer for bodies (for GPU)
  void setTargetCell(Lists lists, Maps flags);                  //!< Set target buffer for cells (for GPU)
  void getTargetBody(Lists &lists);                             //!< Get body values from target buffer (for GPU)
  void getTargetCell(Lists &lists, bool isM);                   //!< Get cell values from target buffer (for GPU)
  void clearBuffers();                                          //!< Clear GPU buffers

  void evalP2P(Bodies &ibodies, Bodies &jbodies, bool onCPU=false);//!< Evaluate all P2P kernels (all pairs)
  void evalP2M(Cells &cells);                                   //!< Evaluate all P2M kernels
  void evalM2M(Cells &cells, Cells &jcells);                    //!< Evaluate all M2M kernels
  void evalM2L(C_iter Ci, C_iter Cj);                           //!< Evaluate on CPU, queue on GPU
  void evalM2L(Cells &cells);                                   //!< Evaluate queued M2L kernels
  void evalM2P(C_iter Ci, C_iter Cj);                           //!< Evaluate on CPU, queue on GPU
  void evalM2P(Cells &cells);                                   //!< Evaluate queued M2P kernels
  void evalP2P(C_iter Ci, C_iter Cj);                           //!< Evaluate on CPU, queue on GPU
  void evalP2P(Cells &cells);                                   //!< Evaluate queued P2P kernels (near field)
  void evalL2L(Cells &cells);                                   //!< Evaluate all L2L kernels
  void evalL2P(Cells &cells);                                   //!< Evaluate all L2P kernels
};

#undef splitFirst
template <class Kernel>
void Evaluator<Kernel>::evalP2P(Bodies &ibodies, Bodies &jbodies, bool) {// Evaluate all P2P kernels
  Xperiodic = 0;                                                // Set periodic coordinate offset
  Cells cells;                                                  // Cells to put target and source bodies
  cells.resize(2);                                              // Resize cells to put target and source bodies
  cells[0].LEAF = ibodies.begin();                              // Iterator of first target leaf
  cells[0].NDLEAF = ibodies.size();                             // Number of target leafs
  cells[1].LEAF = jbodies.begin();                              // Iterator of first source leaf
  cells[1].NDLEAF = jbodies.size();                             // Number of source leafs
  C_iter Ci = cells.begin(), Cj = cells.begin()+1;              // Iterator of target and source cells
  printf("evaluating %d x %d P2P\n",(int)ibodies.size(),(int)jbodies.size());
  K.P2P(Ci,Cj);                                                   // Perform P2P kernel
}

template <class Kernel>
void Evaluator<Kernel>::evalP2M(Cells &cells) {               // Evaluate all P2M kernels
  Log.startTimer("evalP2M");                                        // Start timer
  for( C_iter C=cells.begin(); C!=cells.end(); ++C ) {          // Loop over cells
    for (size_t i=0; i<C->M.size(); i++) C->M[i] = 0;
    //C->M = 0;                                                   //  Initialize multipole coefficients
    for (size_t i=0; i<C->L.size(); i++) C->L[i] = 0;
    //C->L = 0;                                                   //  Initialize local coefficients
    if( C->NCHILD == 0 ) {                                      //  If cell is a twig
      K.P2M(C);                                                   //   Perform P2M kernel
    }                                                           //  Endif for twig
  }                                                             // End loop over cells
  Log.stopTimer("evalP2M");                                         // Stop timer
}

template <class Kernel>
void Evaluator<Kernel>::evalM2M(Cells &cells, Cells &jcells) {// Evaluate all M2M kernels
  Cj0 = jcells.begin();                                         // Set begin iterator
  K.Cj0 = Cj0;
  for( C_iter Ci=cells.begin(); Ci!=cells.end(); ++Ci ) {       // Loop over target cells bottomup
    int level = getLevel(Ci->ICELL);                            // Get current level
    std::stringstream eventName;                                // Declare event name
    eventName << "evalM2M: " << level << "   ";                 // Set event name with level
    Log.startTimer(eventName.str());                                // Start timer
    K.M2M(Ci);                                                    // Perform M2M kernel
    Log.stopTimer(eventName.str());                                 // Stop timer
  }                                                             // End loop target over cells
}

template <class Kernel>
void Evaluator<Kernel>::evalM2L(C_iter Ci, C_iter Cj) {       // Evaluate single M2L kernel
  K.M2L(Ci,Cj);                                                   // Perform M2L kernel
  NM2L++;                                                       // Count M2L kernel execution
}

template <class Kernel>
void Evaluator<Kernel>::evalM2L(Cells &cells) {               // Evaluate queued M2L kernels
  Ci0 = cells.begin();                                          // Set begin iterator
  K.Ci0 = Ci0;
  for( C_iter Ci=cells.begin(); Ci!=cells.end(); ++Ci ) {       // Loop over cells
    int level = getLevel(Ci->ICELL);                            // Get current level
    std::stringstream eventName;                                // Declare event name
    eventName << "evalM2L: " << level << "   ";                 // Set event name with level
    Log.startTimer(eventName.str());                                // Start timer
    while( !listM2L[Ci-Ci0].empty() ) {                         //  While M2L interaction list is not empty
      C_iter Cj = listM2L[Ci-Ci0].back();                       //   Set source cell iterator
      Iperiodic = flagM2L[Ci-Ci0][Cj];                          //   Set periodic image flag
      int I = 0;                                                //   Initialize index of periodic image
      for( int ix=-1; ix<=1; ++ix ) {                           //   Loop over x periodic direction
        for( int iy=-1; iy<=1; ++iy ) {                         //    Loop over y periodic direction
          for( int iz=-1; iz<=1; ++iz, ++I ) {                  //     Loop over z periodic direction
            if( Iperiodic & (1 << I) ) {                        //      If periodic flag is on
              Xperiodic[0] = ix * 2 * R0;                       //       Coordinate offset for x periodic direction
              Xperiodic[1] = iy * 2 * R0;                       //       Coordinate offset for y periodic direction
              Xperiodic[2] = iz * 2 * R0;                       //       Coordinate offset for z periodic direction
              K.M2L(Ci,Cj);                                       //       Perform M2L kernel
            }                                                   //      Endif for periodic flag
          }                                                     //     End loop over x periodic direction
        }                                                       //    End loop over y periodic direction
      }                                                         //   End loop over z periodic direction
      listM2L[Ci-Ci0].pop_back();                               //   Pop last element from M2L interaction list
    }                                                           //  End while for M2L interaction list
    Log.stopTimer(eventName.str());                                 // Stop timer
  }                                                             // End loop over cells topdown
  listM2L.clear();                                              // Clear interaction lists
  flagM2L.clear();                                              // Clear periodic image flags
}

template <class Kernel>
void Evaluator<Kernel>::evalM2P(C_iter Ci, C_iter Cj) {       // Evaluate single M2P kernel
  K.M2P(Ci,Cj);                                                   // Perform M2P kernel
  NM2P++;                                                       // Count M2P kernel execution
}

template <class Kernel>
void Evaluator<Kernel>::evalM2P(Cells &cells) {               // Evaluate queued M2P kernels
  Ci0 = cells.begin();                                          // Set begin iterator
  K.Ci0 = Ci0;
  for( C_iter Ci=cells.begin(); Ci!=cells.end(); ++Ci ) {       // Loop over cells
    int level = getLevel(Ci->ICELL);                            // Get current level
    std::stringstream eventName;                                // Declare event name
    eventName << "evalM2P: " << level << "   ";                 // Set event name with level
    Log.startTimer(eventName.str());                                // Start timer
    while( !listM2P[Ci-Ci0].empty() ) {                         //  While M2P interaction list is not empty
      C_iter Cj = listM2P[Ci-Ci0].back();                       //   Set source cell iterator
      Iperiodic = flagM2P[Ci-Ci0][Cj];                          //   Set periodic image flag
      int I = 0;                                                //   Initialize index of periodic image
      for( int ix=-1; ix<=1; ++ix ) {                           //   Loop over x periodic direction
        for( int iy=-1; iy<=1; ++iy ) {                         //    Loop over y periodic direction
          for( int iz=-1; iz<=1; ++iz, ++I ) {                  //     Loop over z periodic direction
            if( Iperiodic & (1 << I) ) {                        //      If periodic flag is on
              Xperiodic[0] = ix * 2 * R0;                       //       Coordinate offset for x periodic direction
              Xperiodic[1] = iy * 2 * R0;                       //       Coordinate offset for y periodic direction
              Xperiodic[2] = iz * 2 * R0;                       //       Coordinate offset for z periodic direction
              K.M2P(Ci,Cj);                                       //       Perform M2P kernel
            }                                                   //      Endif for periodic flag
          }                                                     //     End loop over x periodic direction
        }                                                       //    End loop over y periodic direction
      }                                                         //   End loop over z periodic direction
      listM2P[Ci-Ci0].pop_back();                               //   Pop last element from M2P interaction list
    }                                                           //  End while for M2P interaction list
    Log.stopTimer(eventName.str());                                 // Stop timer
  }                                                             // End loop over cells topdown
  listM2P.clear();                                              // Clear interaction lists
  flagM2P.clear();                                              // Clear periodic image flags
}

template <class Kernel>
void Evaluator<Kernel>::evalP2P(C_iter Ci, C_iter Cj) {       // Evaluate single P2P kernel
  K.P2P(Ci,Cj);                                                   // Perform P2P kernel
  NP2P++;                                                       // Count P2P kernel execution
}

template <class Kernel>
void Evaluator<Kernel>::evalP2P(Cells &cells) {               // Evaluate queued P2P kernels
  Log.startTimer("evalP2P");                                        // Start timer
  Ci0 = cells.begin();                                          // Set begin iterator
  K.Ci0 = Ci0;
  for( C_iter Ci=cells.begin(); Ci!=cells.end(); ++Ci ) {       // Loop over cells
    while( !listP2P[Ci-Ci0].empty() ) {                         //  While M2P interaction list is not empty
      C_iter Cj = listP2P[Ci-Ci0].back();                       //   Set source cell iterator
      Iperiodic = flagP2P[Ci-Ci0][Cj];                          //   Set periodic image flag
      int I = 0;                                                //   Initialize index of periodic image
      for( int ix=-1; ix<=1; ++ix ) {                           //   Loop over x periodic direction
        for( int iy=-1; iy<=1; ++iy ) {                         //    Loop over y periodic direction
          for( int iz=-1; iz<=1; ++iz, ++I ) {                  //     Loop over z periodic direction
            if( Iperiodic & (1 << I) ) {                        //      If periodic flag is on
              Xperiodic[0] = ix * 2 * R0;                       //       Coordinate offset for x periodic direction
              Xperiodic[1] = iy * 2 * R0;                       //       Coordinate offset for y periodic direction
              Xperiodic[2] = iz * 2 * R0;                       //       Coordinate offset for z periodic direction
              K.P2P(Ci,Cj);                                       //       Perform P2P kernel
            }                                                   //      Endif for periodic flag
          }                                                     //     End loop over x periodic direction
        }                                                       //    End loop over y periodic direction
      }                                                         //   End loop over z periodic direction
      listP2P[Ci-Ci0].pop_back();                               //   Pop last element from M2P interaction list
    }                                                           //  End while for M2P interaction list
  }                                                             // End loop over cells topdown
  listP2P.clear();                                              // Clear interaction lists
  flagP2P.clear();                                              // Clear periodic image flags
  Log.stopTimer("evalP2P");                                         // Stop timer
}

template <class Kernel>
void Evaluator<Kernel>::evalL2L(Cells &cells) {               // Evaluate all L2L kernels
  Ci0 = cells.begin();                                          // Set begin iterator
  K.Ci0 = Ci0;
  for( C_iter Ci=cells.end()-2; Ci!=cells.begin()-1; --Ci ) {   // Loop over cells topdown (except root cell)
    int level = getLevel(Ci->ICELL);                            // Get current level
    std::stringstream eventName;                                // Declare event name
    eventName << "evalL2L: " << level << "   ";                 // Set event name with level
    Log.startTimer(eventName.str());                                // Start timer
    K.L2L(Ci);                                                    // Perform L2L kernel
    Log.stopTimer(eventName.str());                                 // Stop timer
  }                                                             // End loop over cells topdown
}

template <class Kernel>
void Evaluator<Kernel>::evalL2P(Cells &cells) {               // Evaluate all L2P kernels
  Log.startTimer("evalL2P");                                        // Start timer
  for( C_iter C=cells.begin(); C!=cells.end(); ++C ) {          // Loop over cells
    if( C->NCHILD == 0 ) {                                      //  If cell is a twig
      K.L2P(C);                                                   //   Perform L2P kernel
    }                                                           //  Endif for twig
  }                                                             // End loop over cells topdown
  Log.stopTimer("evalL2P");                                         // Stop timer
}

template <class Kernel>
void Evaluator<Kernel>::timeKernels() {                       // Time all kernels for auto-tuning
  Bodies ibodies(1000), jbodies(1000);                          // Artificial bodies
  for( B_iter Bi=ibodies.begin(),Bj=jbodies.begin(); Bi!=ibodies.end(); ++Bi, ++Bj ) {// Loop over artificial bodies
    Bi->X = 0;                                                  //  Set coordinates of target body
    Bj->X = 1;                                                  //  Set coordinates of source body
  }                                                             // End loop over artificial bodies
  Cells cells;                                                  // Artificial cells
  cells.resize(2);                                              // Two artificial cells
  C_iter Ci = cells.begin(), Cj = cells.begin()+1;              // Artificial target & source cell
  Ci->X = 0;                                                    // Set coordinates of target cell
  Ci->NDLEAF = 10;                                              // Number of leafs in target cell
  Ci->LEAF = ibodies.begin();                                   // Leaf iterator in target cell
  Cj->X = 1;                                                    // Set coordinates of source cell
  Cj->NDLEAF = 1000;                                            // Number of leafs in source cell
  Cj->LEAF = jbodies.begin();                                   // Leaf iterator in source cell
  Log.startTimer("P2P kernel");                                     // Start timer
  for( int i=0; i!=1; ++i ) K.P2P(Ci,Cj);                         // Perform P2P kernel
  timeP2P = Log.stopTimer("P2P kernel") / 10000;                    // Stop timer
  Log.startTimer("M2L kernel");                                     // Start timer
  for( int i=0; i!=1000; ++i ) K.M2L(Ci,Cj);                      // Perform M2L kernel
  timeM2L = Log.stopTimer("M2L kernel") / 1000;                     // Stop timer
  Log.startTimer("M2P kernel");                                     // Start timer
  for( int i=0; i!=100; ++i ) K.M2P(Ci,Cj);                       // Perform M2P kernel
  timeM2P = Log.stopTimer("M2P kernel") / 1000;                     // Stop timer
}
