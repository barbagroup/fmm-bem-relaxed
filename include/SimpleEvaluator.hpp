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
#include <Vec.hpp>
#include <Octree.hpp>

#define splitFirst(Ci,Cj) Cj->NCHILD == 0 || (Ci->NCHILD != 0 && Ci->R > Cj->R)

//! Interface between tree and kernel
template <class Kernel>
class SimpleEvaluator
{
public:
  //! Point type
  typedef typename Kernel::point_type point_type;
  //! Multipole expansion type
  typedef typename Kernel::multipole_type multipole_type;
  //! Local expansion type
  typedef typename Kernel::local_type local_type;
  //! Kernel source type
  typedef typename Kernel::charge_type charge_type;

private:
  // kernel & expansions
  Kernel &K;
  std::vector<multipole_type> M;
  std::vector<local_type> L;

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
    //K.M2P(*Cj,M[Cj->ICELL],*Ci);
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

public:
  //! Constructor
  Evaluator() : R0(0), Icenter(1 << 13), NP2P(0), NM2P(0), NM2L(0), K(Kernel()), M(0), L(0) {};
  Evaluator(Kernel& k) : R0(0),  Icenter(1 << 13), NP2P(0), NM2P(0), NM2L(0), K(k), M(0), L(0) {};
  //! Destructor
  ~Evaluator() {}

  // upward sweep using new tree structure
  void upward(Octree<point_type>& otree, std::vector<charge_type>& charges)
  {
    M.resize(otree.boxes());
    L.resize(otree.boxes());

    unsigned lowest_level = otree.levels();
    printf("lowest level in tree: %d\n",(int)lowest_level);

    // P2M calls
    for (auto it=otree.box_begin(); it!=otree.box_end(); ++it)
    {
      if (it->is_leaf()) {
        K.init_multipole(M[it->index()], 0.);
        K.init_local(L[it->index()],0.);

        K.P2M(it->body_begin(),it->body_end(),charges.begin()+it->body_begin()->index(),it->center(),M[it->index()]);
      }
    }

    // M2M calls
    for (auto it=otree.box_begin(); it!=otree.box_end(); ++it)
    {
      if (!it->is_leaf()) {
        K.init_multipole(M[it->index()],0.);
        K.init_local(L[it->index()],0.);
  
        for (auto child_it=it->child_begin(); child_it!=it->child_end(); ++child_it)
        {
          auto translation = it->center() - child_it->center();
          K.M2M(M[child_it->index()],M[it->index()],translation);
        }
      }
    }
  }

  //! Use multipole acceptance criteria to determine whether to approximate, do P2P, or subdivide
  void interact(C_iter Ci, C_iter Cj, PairQueue& pairQueue) {
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
  void traverse(Cells& cells, Cells& jcells) {
    C_iter root = cells.end() - 1;                              // Iterator for root target cell
    C_iter jroot = jcells.end() - 1;                            // Iterator for root source cell
    if( IMAGES != 0 ) {                                         // If periodic boundary condition
      jroot = jcells.end() - 1 - 26 * 27 * (IMAGES - 1);        //  The root is not at the end
    }                                                           // Endif for periodic boundary condition
    Ci0 = cells.begin();                                        // Set begin iterator for target cells
    Cj0 = jcells.begin();                                       // Set begin iterator for source cells
    // Non-periodic
    Iperiodic = Icenter;                                      //  Set periodic image flag to center
    Xperiodic = 0;                                            //  Set periodic coordinate offset
    Pair pair(root,jroot);                                    //  Form pair of root cells
    traverseQueue(pair);                                      //  Traverse a pair of trees
  }

  //! Downward phase (M2L,M2P,P2P,L2L,L2P evaluation)
  void downward(Cells& cells, Cells& jcells, bool periodic=true) {
    Log.startTimer("Traverse");                                     // Start timer
    traverse(cells,jcells);                                     // Traverse tree to get interaction list
    Log.stopTimer("Traverse");                             // Stop timer & print

    evalL2L();
    evalL2P();
  }

  static void evalP2P(const Kernel& K, Bodies& ibodies, Bodies& jbodies) {
    Xperiodic = 0;                                                // Set periodic coordinate offset
    Cells cells;                                                  // Cells to put target and source bodies
    cells.resize(2);                                              // Resize cells to put target and source bodies
    cells[0].LEAF = ibodies.begin();                              // Iterator of first target leaf
    cells[0].NDLEAF = ibodies.size();                             // Number of target leafs
    cells[1].LEAF = jbodies.begin();                              // Iterator of first source leaf
    cells[1].NDLEAF = jbodies.size();                             // Number of source leafs
    C_iter Ci = cells.begin(), Cj = cells.begin()+1;              // Iterator of target and source cells
    printf("evaluating %d x %d P2P\n",(int)ibodies.size(),(int)jbodies.size());
    K.P2P(Ci,Cj);
  }
};

#undef splitFirst

template <class Kernel>
void Evaluator<Kernel>::evalP2P(Bodies& ibodies, Bodies& jbodies, bool) {// Evaluate all P2P kernels
  evalP2P(K, ibodies, jbodies); // Perform P2P kernel
}

template <class Kernel>
void Evaluator<Kernel>::evalM2L(C_iter Ci, C_iter Cj) {       // Evaluate single M2L kernel
  vect translation = Ci->X - Cj->X;
  K.M2L(M[Cj->ICELL],L[Ci->ICELL],translation);
  NM2L++;                                                       // Count M2L kernel execution
}

template <class Kernel>
void Evaluator<Kernel>::evalM2P(C_iter Ci, C_iter Cj) {       // Evaluate single M2P kernel
  // K.M2P(*Cj,M[Cj->ICELL],*Ci);
  K.M2P(Cj->X,M[Cj->ICELL],*Ci);
  NM2P++;                                                       // Count M2P kernel execution
}

template <class Kernel>
void Evaluator<Kernel>::evalP2P(C_iter Ci, C_iter Cj) {       // Evaluate single P2P kernel
  K.P2P(Ci,Cj);                                                   // Perform P2P kernel
  NP2P++;                                                       // Count P2P kernel execution
}

template <class Kernel>
void Evaluator<Kernel>::evalL2L(Cells &cells) {               // Evaluate all L2L kernels
  // for each cell in range, get parent & perform L2L
  for( C_iter Ci=cells.end()-2; Ci!=cells.begin()-1; --Ci )
  {
    C_iter parent = cells.begin()+Ci->PARENT;
    vect dist = Ci->X-parent->X;
    K.L2L(L[parent->ICELL],L[Ci->ICELL],dist);
  }
}

template <class Kernel>
void Evaluator<Kernel>::evalL2P(Cells &cells) {               // Evaluate all L2P kernels
  Log.startTimer("evalL2P");                                        // Start timer
  for( C_iter C=cells.begin(); C!=cells.end(); ++C ) {          // Loop over cells
    if( C->NCHILD == 0 ) {                                      //  If cell is a twig
      K.L2P(*C,L[C->ICELL]);
    }                                                           //  Endif for twig
  }                                                             // End loop over cells topdown
  Log.stopTimer("evalL2P");                                         // Stop timer
}

