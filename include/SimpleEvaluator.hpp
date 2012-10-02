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

    // For the lowest level up to the highest level
    for (unsigned L = otree.levels()-1; L != 1; --L) {
      // For all boxes at this level
      auto b_end = otree.box_end(L);
      for (auto bit = otree.box_begin(L); bit != b_end; ++bit) {
        auto box = *bit;

        // Initialize box data
        unsigned idx = box.index();
        double box_size = box.side_length();
        K.init_multipole(M[idx], box_size);
        K.init_local(L[idx], box_size);

        if (box.is_leaf()) {
          // If leaf, make P2M calls

          // For all the bodies, P2M
          auto p_begin = box.body_begin();
          K.P2M(p_begin, box.body_end(),
                charges.begin() + p_begin->index(),
                box.center(),
                M[idx]);
        } else {
          // If not leaf, make M2M calls

          // For all the children, M2M
          auto c_end = box.child_end();
          for (auto cit = box.child_begin(); cit != c_end; ++cit) {
            auto cbox = *cit;
            auto translation = box.center() - cbox.center();

            K.M2M(M[cbox.index()], M[idx], translation);
          }
        }
      }
    }
  }

  template <typename BOX, typename Q>
  void interact(const BOX& b1, const BOX& b2, Q& pairQ) {
    point_type r0 = b1.center() - b2.center();
    double r0_norm = std::sqrt(norm(r0));
    if (r0_norm * THETA > b1.side_length() + b2.side_length()) {
      // These boxes satisfy the multipole acceptance criteria
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
    } else if(b1.is_leaf() && b2.is_leaf()) {
      evalP2P(b1,b2);
    } else {
      pairQ.push_back(std::make_pair(b1,b2));
    }
  }


  void downward(Octree<point_type>& octree) {

    typedef Octree<point_type>::Box Box;
    typedef std::pair<Box, Box> box_pair;
    std::deque<box_pair> pairQ;

    // Queue based tree traversal for P2P, M2P, and/or M2L operations
    pairQ.push_back(box_pair(octree.root(), octree.root()));

    while (!pairQ.empty()) {
      box_pair boxes = pairQ.front();
      pairQ.pop_front();
      bool is_leaf1 = boxes.first.is_leaf();
      bool is_leaf2 = boxes.second.is_leaf();

      if (is_leaf2 || (!is_leaf1 && boxes.first > boxes.second)) {
        // Split the first box into children and interact
        auto c_end = boxes.first.child_end();
        for (auto cit = boxes.first.child_begin(); cit != c_end; ++cit)
          interact(*cit, boxes.second, pairQ);
      } else {
        // Split the second box into children and interact
        auto c_end = boxes.second.child_end();
        for (auto cit = boxes.second.child_begin(); cit != c_end; ++cit)
          interact(boxes.first, *cit, pairQ);
      }
    }

    //

    // For the highest level down to the lowest level
    for (unsigned L = 2; L < octree.levels(); ++L) {
      // For all boxes at this level
      auto b_end = otree.box_end(L);
      for (auto bit = otree.box_begin(L); bit != b_end; ++bit) {
        auto box = *bit;
        unsigned idx = box.index();

        // Initialize box data
        if (box.is_leaf()) {
          // If leaf, make L2P calls

          // For all the bodies, L2P
          auto p_begin = box.body_begin();
          K.L2P(p_begin, box.body_end(),
                box.center(),
                L[idx]);
        } else {
          // If not leaf, make L2L calls

          // For all the children, L2L
          auto c_end = box.child_end();
          for (auto cit = box.child_begin(); cit != c_end; ++cit) {
            auto cbox = *cit;
            auto translation = cbox.center() - box.center();

            K.L2L(L[idx], L[cbox.index()], translation);
          }
        }
      }
    }
  }

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
    Ci0 = cells.begin();                                        // Set begin iterator for target cells
    Cj0 = jcells.begin();                                       // Set begin iterator for source cells
    // Non-periodic
    Iperiodic = Icenter;                                      //  Set periodic image flag to center
    Xperiodic = 0;                                            //  Set periodic coordinate offset
    Pair pair(root,jroot);                                    //  Form pair of root cells
    traverseQueue(pair);                                      //  Traverse a pair of trees
  }

  //! Traverse a pair of trees using a queue
  void traverseQueue(Pair pair) {
    PairQueue pairQueue;                                        // Queue of interacting cell pairs
    pairQueue.push_back(pair);                                  // Push pair to queue
    while( !pairQueue.empty() ) {                               // While dual traversal queue is not empty
      pair = pairQueue.front();                                 //  Get interaction pair from front of queue
      pairQueue.pop_front();                                    //  Pop dual traversal queue
      if(splitFirst(pair.first,pair.second)) {                  //  If first cell is larger
        // Cj->NCHILD == 0 || (Ci->NCHILD != 0 && Ci->R > Cj->R)
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

