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
#include <Logger.hpp>
#include <Sorter.hpp>
#include <Evaluator.hpp>
#include <Kernel.hpp>

extern Logger Log;
extern Sorter sort;

//! Base class for tree structure
class TreeStructure
{
private:
  vect X0;
  real R0;

public:
  Bodies buffer;                                                //!< Buffer for MPI communication & sorting

  TreeStructure() : X0(0), R0(0) {};
  TreeStructure(vect &X, real &R) : X0(X), R0(R) {};

  void init(vect &x, real &r)
  {
    X0 = x;
    R0 = r;
  }

  // topdown tree construction
  template <class Kernel>
  void topdown(Bodies &bodies, Cells &cells, Evaluator<Kernel> *eval)
  {
    grow(bodies);
    setIndex();

    buffer.resize(bodies.size());
    sort.sortBodies(bodies,buffer,false);

    Cells twigs;
    bodies2twigs(bodies,twigs);
    eval->evalP2M(twigs);

    // here twigs contains all twig cells for P2M

    Cells sticks;
    twigs2cells(twigs,cells,sticks);
    eval->evalM2M(cells,cells);
  }

private:

  int getLevel(bigint index)
  {
    int i = index;                                              // Copy to dummy index
    int level = -1;                                             // Initialize level counter
    while( i >= 0 ) {                                           // While cell index is non-negative
      level++;                                                  //  Increment level
      i -= 1 << 3*level;                                        //  Subtract number of cells in that level
    }                                                           // End while loop for cell index
    return level;                                               // Return the level
  }

  //! Get parent cell index from current cell index
  bigint getParent(bigint index) {
    int level = getLevel(index);                                // Get level from cell index
    bigint cOff = ((1 << 3 *  level   ) - 1) / 7;               // Cell index offset of current level
    bigint pOff = ((1 << 3 * (level-1)) - 1) / 7;               // Cell index offset of parent level
    bigint i = ((index-cOff) >> 3) + pOff;                      // Cell index of parent cell
    return i;                                                   // Return cell index of parent cell
  }

  //! Merge sticks with cells (levelwise)
  void unique(Cells &cells, Cells &sticks, int begin, int &end) {
    int c_old = begin;                                          // Initialize old cell counter
    for( int c=begin; c!=end; ++c ) {                           // Loop over cells in level
      if( cells[c].ICELL != cells[c_old].ICELL ) {              //  If current cell index is different from previous
        c_old = c;                                              //   Update old cell counter
      } else if( c != c_old ) {                                 //  If cell index is repeated
        if( cells[c].NCHILD != 0 ) {                            //   Stick-cell collision
          cells[c_old].NCHILD = cells[c].NCHILD;                //    Copy number of children
          cells[c_old].NCLEAF = cells[c].NCLEAF;                //    Copy number of leafs
          cells[c_old].NDLEAF = cells[c].NDLEAF;                //    Copy number of leafs
          cells[c_old].PARENT = cells[c].PARENT;                //    Copy parent link
          cells[c_old].CHILD = cells[c].CHILD;                  //    Copy child link
          cells[c_old].LEAF = cells[c].LEAF;                    //    Copy iterator of first leaf
          sticks.push_back(cells[c_old]);                       //    Push stick into vector
        }                                                       //   Endif for collision type
        // cells[c_old].M += cells[c].M;                           //   Accumulate multipole
        for (size_t i=0; i<cells[c_old].M.size(); i++) cells[c_old].M[i] += cells[c].M[i];
        cells.erase(cells.begin()+c);                           //   Erase colliding cell
        c--;                                                    //   Decrement counter to account for erase
        end--;                                                  //   Decrement end to account for erase
      }                                                         //  Endif for repeated cell index
    }                                                           // End loop over cells in level
  }

  //! Form parent-child mutual link
  void linkParent(Cells &cells, int &begin, int &end) {
    Cell parent = Cell();                                                // Parent cell
    //parent.alloc_multipole(Kernel<equation>::multipole_size(getLevel(parent.ICELL)));
    //parent.alloc_local(Kernel<equation>::local_size(getLevel(parent.ICELL)));
    parent.alloc_multipole(NTERM);
    parent.alloc_local(NTERM);
    Cells parents;                                              // Parent cell vector;
    int oldend = end;                                           // Save old end counter
    parent.ICELL = getParent(cells[begin].ICELL);               // Set cell index
    // parent.M = 0;                                               // Initialize multipole coefficients
    for (size_t i=0; i<parent.M.size(); i++) parent.M[i] = 0;
    // parent.L = 0;                                               // Initlalize local coefficients
    for (size_t i=0; i<parent.L.size(); i++) parent.L[i] = 0;
    parent.NCLEAF = parent.NDLEAF = parent.NCHILD = 0;          // Initialize NCLEAF, NDLEAF, & NCHILD
    parent.LEAF = cells[begin].LEAF;                            // Set pointer to first leaf
    parent.CHILD = begin;                                       // Link to child
    getCenter(parent);                                          // Set cell center and radius
    for( int i=begin; i!=oldend; ++i ) {                        // Loop over cells at this level
      if( getParent(cells[i].ICELL) != parent.ICELL ) {         //  If it belongs to a new parent cell
        cells.push_back(parent);                                //   Push cells into vector
        end++;                                                  //   Increment cell counter
        parent.ICELL = getParent(cells[i].ICELL);               //   Set cell index
        // parent.M = 0;                                           //   Initialize multipole coefficients
        for (size_t l=0; l<parent.M.size(); l++) parent.M[l] = 0;
        // parent.L = 0;                                           //   Initialize local coefficients
        for (size_t l=0; l<parent.L.size(); l++) parent.L[l] = 0;
        parent.NCLEAF = parent.NDLEAF = parent.NCHILD = 0;      //   Initialize NCLEAF, NDLEAF, & NCHILD
        parent.LEAF = cells[i].LEAF;                            //   Set pointer to first leaf
        parent.CHILD = i;                                       //   Link to child
        getCenter(parent);                                      //   Set cell center and radius
      }                                                         //  Endif for new parent cell
      for( int c=0; c!=cells[i].NCHILD; ++c ) {                 //  Loop over child cells
        (cells.begin()+cells[i].CHILD+c)->PARENT = i;           //   Link child to current
      }                                                         //  End loop over child cells
      cells[i].PARENT = end;                                    //  Link to current to parent
      parent.NDLEAF += cells[i].NDLEAF;                         //  Add nleaf of child to parent
      parents.push_back(parent);                                //  Push parent cell into vector
      parent = parents.back();                                  //  Copy information from vector
      parents.pop_back();                                       //  Pop parent cell from vector
      parent.NCHILD++;                                          //  Increment child counter
    }                                                           // End loop over cells at this level
    cells.push_back(parent);                                    // Push cells into vector
    end++;                                                      // Increment cell counter
    begin = oldend;                                             // Set new begin index to old end index
  }

protected:
  //! Get cell center and radius from cell index
  void getCenter(Cell &cell) {
    int level = getLevel(cell.ICELL);                           // Get level from cell index
    bigint index = cell.ICELL - ((1 << 3*level) - 1) / 7;       // Subtract cell index offset of current level
    cell.R = R0 / (1 << level);                                 // Cell radius
    int d = level = 0;                                          // Initialize dimension and level
    vec<3,int> nx = 0;                                          // Initialize 3-D cell index
    while( index != 0 ) {                                       // Deinterleave bits while index is nonzero
      nx[d] += (index % 2) * (1 << level);                      //  Add deinterleaved bit to 3-D cell index
      index >>= 1;                                              //  Right shift the bits
      d = (d+1) % 3;                                            //  Increment dimension
      if( d == 0 ) level++;                                     //  If dimension is 0 again, increment level
    }                                                           // End while loop for deinterleaving bits
    for( d=0; d!=3; ++d ) {                                     // Loop over dimensions
      cell.X[d] = (X0[d]-R0) + (2 *nx[d] + 1) * cell.R;         //  Calculate cell center from 3-D cell index
    }                                                           // End loop over dimensions
  }

public:
  //! Group bodies into twig cells
  void bodies2twigs(Bodies &bodies, Cells &twigs) {
    Log.startTimer("Bodies2twigs");                                 // Start timer
    int nleaf = 0;                                              // Initialize number of leafs
    bigint index = bodies[0].ICELL;                             // Initialize cell index
    B_iter firstLeaf = bodies.begin();                          // Initialize body iterator for first leaf
    Cell cell = Cell();                                                  // Cell structure
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over bodies
      if( B->ICELL != index ) {                                 //  If it belongs to a new cell
        cell.NCLEAF = nleaf;                                    //   Set number of child leafs
        cell.NDLEAF = nleaf;                                    //   Set number of decendant leafs
        cell.CHILD  = 0;                                        //   Set pointer offset to first child
        cell.NCHILD = 0;                                        //   Set number of child cells
        cell.ICELL  = index;                                    //   Set cell index
        cell.LEAF   = firstLeaf;                                //   Set pointer to first leaf
        getCenter(cell);                                        //   Set cell center and radius

        //cell.alloc_multipole(Kernel<equation>::multipole_size(getLevel(cell.ICELL)));
        //cell.alloc_local(Kernel<equation>::multipole_size(getLevel(cell.ICELL)));
        cell.alloc_multipole(NTERM);
        cell.alloc_local(NTERM);

        twigs.push_back(cell);                                  //   Push cells into vector
        firstLeaf = B;                                          //   Set new first leaf
        nleaf = 0;                                              //   Reset number of bodies
        index = B->ICELL;                                       //   Set new cell
      }                                                         //  Endif for new cell
      nleaf++;                                                  //  Increment body counter
    }                                                           // End loop over bodies
    cell.NCLEAF = nleaf;                                        // Set number of child leafs
    cell.NDLEAF = nleaf;                                        // Set number of decendant leafs
    cell.CHILD  = 0;                                            //   Set pointer offset to first child
    cell.NCHILD = 0;                                            // Set number of child cells
    cell.ICELL  = index;                                        // Set cell index
    cell.LEAF   = firstLeaf;                                    // Set pointer to first leaf
    getCenter(cell);                                            // Set cell center and radius
    twigs.push_back(cell);                                      // Push cells into vector
    Log.stopTimer("Bodies2twigs");                              // Stop timer & print
    // evalP2M(twigs);                                             // Evaluate all P2M kernels

  }

  //! Link twigs bottomup to create all cells in tree
  void twigs2cells(Cells &twigs, Cells &cells, Cells &sticks) {
    int begin = 0, end = 0;                                     // Initialize range of cell vector
    int level = getLevel(twigs.back().ICELL);                   // Initialize level of tree
    Log.startTimer("Sort resize");                                  // Start timer
    Cells cbuffer;                                              // Sort buffer for cells
    cbuffer.resize(2*twigs.size());                             // Resize sort buffer for cells
    Log.stopTimer("Sort resize");                                   // Stop timer
    while( !twigs.empty() ) {                                   // Keep poppig twigs until the vector is empty
      while( getLevel(twigs.back().ICELL) != level ) {          //  While cell belongs to a higher level
        sort.sortCells(cells,cbuffer,false,begin,end);               //   Sort cells at this level
        Log.startTimer("Twigs2cells");                              //   Start timer
        unique(cells,sticks,begin,end);                         //   Get rid of duplicate cells
        linkParent(cells,begin,end);                            //   Form parent-child mutual link
        level--;                                                //   Go up one level
        Log.stopTimer("Twigs2cells");                               //   Stop timer
      }                                                         //  End while for higher level
      Log.startTimer("Twigs2cells");                                //  Start timer
      cells.push_back(twigs.back());                            //  Push cells into vector
      twigs.pop_back();                                         //  Pop twigs from vector
      end++;                                                    //  Increment cell counter
      Log.stopTimer("Twigs2cells");                                 //  Stop timer
    }                                                           // End while for popping twigs
    for( int l=level; l>0; --l ) {                              // Once all the twigs are done, do the rest
      sort.sortCells(cells,cbuffer,false,begin,end);                 //  Sort cells at this level
      Log.startTimer("Twigs2cells");                                //  Start timer
      unique(cells,sticks,begin,end);                           //  Get rid of duplicate cells
      linkParent(cells,begin,end);                              //  Form parent-child mutual link
      Log.stopTimer("Twigs2cells");                                 //  Stop timer
    }                                                           // End loop over levels
    Log.startTimer("Twigs2cells");                                  // Start timer
    unique(cells,sticks,begin,end);                             // Just in case there is a collision at root
    Log.stopTimer("Twigs2cells");                          // Stop timer & print
    // evalM2M(cells,cells);                                       // Evaluate all M2M kernels
  }

#if 0
  //! Downward phase (M2L,M2P,P2P,L2L,L2P evaluation)
  void downward(Cells &cells, Cells &jcells, bool periodic=true) {
#if HYBRID
    timeKernels();                                              // Time all kernels for auto-tuning
#endif
    for( C_iter C=cells.begin(); C!=cells.end(); ++C ) {        // Initialize local coefficients
      //C->L = 0;
      for (size_t i=0; i<C->L.size(); i++) C->L[i] = 0;
    }
    if( IMAGES != 0 ) {                                         // If periodic boundary condition
      Log.startTimer("Upward P");                                   //  Start timer
      upwardPeriodic(jcells);                                   //  Upward phase for periodic images
      Log.stopTimer("Upward P",printNow);                           //  Stop timer & print
    }                                                           // Endif for periodic boundary condition
    Log.startTimer("Traverse");                                     // Start timer
    traverse(cells,jcells);                                     // Traverse tree to get interaction list
    Log.stopTimer("Traverse",printNow);                             // Stop timer & print
    if( IMAGES != 0 && periodic ) {                             // If periodic boundary condition
      Log.startTimer("Traverse P");                                 // Start timer
      traversePeriodic(cells,jcells);                           // Traverse tree for periodic images
      Log.stopTimer("Traverse P",printNow);                         // Stop timer & print
    }                                                           // Endif for periodic boundary condition
    evalL2L(cells);                                             // Evaluate all L2L kernels
    evalL2P(cells);                                             // Evaluate all L2P kernels
    if(printNow) std::cout << "P2P: "  << NP2P
                           << " M2P: " << NM2P
                           << " M2L: " << NM2L << std::endl;
  }
#endif

  // topdown tree construction stuff
  //! Nodes are primitive cells
  struct Node {
    int LEVEL;                                                  //!< Level of node
    int ICHILD;                                                 //!< Flag of empty child nodes
    int NLEAF;                                                  //!< Number of leafs in node
    bigint I;                                                   //!< Cell index
    bigint CHILD[8];                                            //!< Iterator offset of child nodes
    B_iter LEAF[NCRIT];                                         //!< Iterator for leafs
    vect X;                                                     //!< Node center
    real R;                                                     //!< Node radius
  };  
  std::vector<Node> nodes;                                      //!< Nodes in the tree

//! Calculate octant from position
  int getOctant(const vect X, int i) {
    int octant = 0;                                             // Initialize octant
    for( int d=0; d!=3; ++d ) {                                 // Loop over dimensions
      octant += (X[d] > nodes[i].X[d]) << d;                    //  interleave bits and accumulate octant
    }                                                           // End loop over dimensions
    return octant;                                              // Return octant
  }

//! Add child node and link it
  void addChild(const int octant, int i) {
    bigint pOff = ((1 << 3* nodes[i].LEVEL   ) - 1) / 7;        // Parent cell index offset
    bigint cOff = ((1 << 3*(nodes[i].LEVEL+1)) - 1) / 7;        // Current cell index offset
    vect x = nodes[i].X;                                        // Initialize new center position with old center
    real r = nodes[i].R/2;                                      // Initialize new size
    for( int d=0; d!=3; ++d ) {                                 // Loop over dimensions
      x[d] += r * (((octant & 1 << d) >> d) * 2 - 1);           //  Calculate new center position
    }                                                           // End loop over dimensions
    Node node;                                                  // Node structure
    node.NLEAF = node.ICHILD = 0;                               // Initialize child node counters
    node.X = x;                                                 // Initialize child node center
    node.R = r;                                                 // Initialize child node radius
    node.LEVEL = nodes[i].LEVEL + 1;                            // Level of child node
    node.I = ((nodes[i].I-pOff) << 3) + octant + cOff;          // Cell index of child node
    nodes[i].ICHILD |= (1 << octant);                           // Flip bit of octant
    nodes[i].CHILD[octant] = nodes.end()-nodes.begin();         // Link child node to its parent
    nodes.push_back(node);                                      // Push child node into vector
  }

//! Add leaf to node
  void addLeaf(B_iter b, int i) {
    nodes[i].LEAF[nodes[i].NLEAF] = b;                          // Assign body iterator to leaf
    nodes[i].NLEAF++;                                           // Increment leaf counter
  }

//! Split node and reassign leafs to child nodes
  void splitNode(int i) {
    for( int l=0; l!=NCRIT; ++l ) {                             // Loop over leafs in parent node
      int octant = getOctant(nodes[i].LEAF[l]->X,i);            //  Find the octant where the body belongs
      if( !(nodes[i].ICHILD & (1 << octant)) ) {                //  If child doesn't exist in this octant
        addChild(octant,i);                                     //   Add new child to list
      }                                                         //  Endif for octant
      int c = nodes[i].CHILD[octant];                           //  Set counter for child node
      addLeaf(nodes[i].LEAF[l],c);                              //  Add leaf to child
      if( nodes[c].NLEAF >= NCRIT ) {                           //  If there are still too many leafs
        splitNode(c);                                           //   Split the node into smaller ones
      }                                                         //  Endif for too many leafs
    }                                                           // End loop over leafs
  }

//! Traverse tree
  void traverse(typename std::vector<Node>::iterator N) {
    if( N->NLEAF >= NCRIT ) {                                   // If node has children
      for( int i=0; i!=8; ++i ) {                               // Loop over children
        if( N->ICHILD & (1 << i) ) {                            //  If child exists in this octant
          traverse(nodes.begin()+N->CHILD[i]);                  //   Recursively search child node
        }                                                       //  Endif for octant
      }                                                         // End loop over children
    } else {                                                    //  If child doesn't exist
      for( int i=0; i!=N->NLEAF; ++i ) {                        //   Loop over leafs
        N->LEAF[i]->ICELL = N->I;                               //    Store cell index in bodies
      }                                                         //   End loop over leafs
    }                                                           //  Endif for child existence
  }

private:
//! Grow tree from root
  void grow(Bodies &bodies) {
    Log.startTimer("Grow tree");                                    // Start timer
    int octant;                                                 // In which octant is the body located?
    Node node;                                                  // Node structure
    node.LEVEL = node.NLEAF = node.ICHILD = node.I = 0;         // Initialize root node counters
    node.X = X0;                                                // Initialize root node center
    node.R = R0;                                                // Initialize root node radius
    nodes.push_back(node);                                      // Push child node into vector
    for( B_iter B=bodies.begin(); B!=bodies.end(); ++B ) {      // Loop over bodies
      int i = 0;                                                //  Reset node counter
      while( nodes[i].NLEAF >= NCRIT ) {                        //  While the nodes have children
        nodes[i].NLEAF++;                                       //   Increment the cumulative leaf counter
        octant = getOctant(B->X,i);                             //   Find the octant where the body belongs
        if( !(nodes[i].ICHILD & (1 << octant)) ) {              //   If child doesn't exist in this octant
          addChild(octant,i);                                   //    Add new child to list
        }                                                       //   Endif for child existence
        i = nodes[i].CHILD[octant];                             //    Update node iterator to child
      }                                                         //  End while loop
      addLeaf(B,i);                                             //  Add body to node as leaf
      if( nodes[i].NLEAF >= NCRIT ) {                           //  If there are too many leafs
        splitNode(i);                                           //   Split the node into smaller ones
      }                                                         //  Endif for splitting
    }                                                           // End loop over bodies
    Log.stopTimer("Grow tree");                            // Stop timer
  }

//! Store cell index of all bodies
  void setIndex() {
    Log.startTimer("Set index");                                    // Start timer
    traverse(nodes.begin());                                    // Traverse tree
    Log.stopTimer("Set index");                            // Stop timer 
  }
};

