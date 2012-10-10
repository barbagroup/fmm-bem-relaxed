#pragma once

/**
 * FMM evaluator based on dual tree traversal
 */

#include <TransformIterator.hpp>

//! forward definition
template <typename Tree, typename Kernel>
class EvaluatorBase;

template <typename Tree, typename  Kernel>
class EvaluatorTreecode : public EvaluatorBase<Tree,Kernel>
{
  typedef EvaluatorBase<Tree,Kernel> Base;
 public:
  typedef typename Base::charge_type charge_type;
  typedef typename Base::result_type result_type;
  typedef typename Base::point_type  point_type;

  template <typename BOX, typename Q>
  void interact(const BOX& b1, const BOX& b2, Q& pairQ) {
    double r0_norm = std::sqrt(normSq(b1.center() - b2.center()));
    if (r0_norm * Base::THETA > b1.side_length()/2 + b2.side_length()/2) {
      // These boxes satisfy the multipole acceptance criteria
      Base::evalM2P(b1,b2);
    } else if(b1.is_leaf() && b2.is_leaf()) {
      Base::evalP2P(b2,b1);
    } else {
      pairQ.push_back(std::make_pair(b1,b2));
    }
  }

 public:
  //! constructor
  EvaluatorTreecode(Tree& t, Kernel& k, double theta)
        : EvaluatorBase<Tree,Kernel>(t,k,theta) {};

  //! upward sweep
  void upward(const std::vector<charge_type>& charges)
  {
    (void) charges; // Quiet warning

    // EvaluatorBase<Tree,Kernel>::M.resize(10);
    unsigned lowest_level = Base::tree.levels();
    printf("lowest level in tree: %d\n",(int)lowest_level);

    // For the lowest level up to the highest level
    for (unsigned l = Base::tree.levels()-1; l != 0; --l) {
      // For all boxes at this level
      auto b_end = Base::tree.box_end(l);
      for (auto bit = Base::tree.box_begin(l); bit != b_end; ++bit) {
        auto box = *bit;

        // Initialize box data
        unsigned idx = box.index();
        double box_size = box.side_length();
        Base::K.init_multipole(Base::M[idx], box_size);
        Base::K.init_local(Base::L[idx], box_size);

        if (box.is_leaf()) {
          // If leaf, make P2M calls
          Base::evalP2M(box);
        } else {
          // If not leaf, make M2M calls
          Base::evalM2M(box);
        }
      }
    }
  }

  //! Box-Box interactions
  void interactions(std::vector<result_type>& results)
  {
    (void) results; // Quiet warning

    typedef typename Octree<point_type>::Box Box;
    typedef typename std::pair<Box, Box> box_pair;
    std::deque<box_pair> pairQ;

    // Queue based tree traversal for P2P, M2P, and/or M2L operations
    pairQ.push_back(box_pair(Base::tree.root(), Base::tree.root()));

    while (!pairQ.empty()) {
      auto b1 = pairQ.front().first;
      auto b2 = pairQ.front().second;
      pairQ.pop_front();

      if (b2.is_leaf() || (!b1.is_leaf() && b1.side_length() > b2.side_length())) {
        // Split the first box into children and interact
        auto c_end = b1.child_end();
        for (auto cit = b1.child_begin(); cit != c_end; ++cit)
          interact(*cit, b2, pairQ);
      } else {
        // Split the second box into children and interact
        auto c_end = b2.child_end();
        for (auto cit = b2.child_begin(); cit != c_end; ++cit)
          interact(b1, *cit, pairQ);
      }
    }
  }

  //! downward sweep
  void downward(std::vector<result_type>& results)
  {
    // not needed for treecode
    (void) results; // Quiet warning
  }

  //! evaluator name
  std::string name() {
    return "Treecode";
  }
};

