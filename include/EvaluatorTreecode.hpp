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
 private:
  typedef typename EvaluatorBase<Tree,Kernel>::charge_type charge_type;
  typedef typename EvaluatorBase<Tree,Kernel>::result_type result_type;
  typedef typename EvaluatorBase<Tree,Kernel>::point_type  point_type;

  typedef EvaluatorBase<Tree,Kernel> Base;

  typename std::vector<result_type>::iterator results_begin;
  typename std::vector<charge_type>::const_iterator charges_begin;

  /** One-sided P2P!!
   */
  void evalP2P(const typename Octree<point_type>::Box& b1, 
               const typename Octree<point_type>::Box& b2) {
    // Point iters
    auto body2point = [](const typename Octree<point_type>::Body& b) { return b.point(); };
    auto p1_begin = make_transform_iterator(b1.body_begin(), body2point);
    auto p1_end   = make_transform_iterator(b1.body_end(),   body2point);
    auto p2_begin = make_transform_iterator(b2.body_begin(), body2point);
    auto p2_end   = make_transform_iterator(b2.body_end(),   body2point);

    // Charge iters
    auto c1_begin = charges_begin + b1.body_begin()->index();

    // Result iters
    auto r2_begin = results_begin + b2.body_begin()->index();

    printf("P2P: %d to %d\n",b1.index(),b2.index());

    Base::K.P2P(p1_begin, p1_end, c1_begin,
                p2_begin, p2_end,
                r2_begin);
  }

  void evalM2P(const typename Octree<point_type>::Box& b1, 
               const typename Octree<point_type>::Box& b2) 
  {
    // Target point iters
    auto body2point = [](const typename Octree<point_type>::Body& b) { return b.point(); };
    auto t_begin = make_transform_iterator(b2.body_begin(), body2point);
    auto t_end   = make_transform_iterator(b2.body_end(), body2point);

    // Target result iters
    auto r_begin = results_begin + b2.body_begin()->index();

    printf("M2P: %d to %d\n", b1.index(), b2.index());

    Base::K.M2P(Base::M[b1.index()], b1.center(),
                t_begin, t_end,
                r_begin);
  }

  void evalM2L(const typename Octree<point_type>::Box& b1, 
               const typename Octree<point_type>::Box& b2) 
  {
    // unused -- quiet compiler warnings
    (void) b1;
    (void) b2;
  }

  template <typename BOX, typename Q>
  void interact(const BOX& b1, const BOX& b2, Q& pairQ) {
    double r0_norm = std::sqrt(normSq(b1.center() - b2.center()));
    if (r0_norm * Base::THETA > b1.side_length()/2 + b2.side_length()/2) {
      // These boxes satisfy the multipole acceptance criteria
      evalM2P(b1,b2);
    } else if(b1.is_leaf() && b2.is_leaf()) {
      evalP2P(b2,b1);
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

    // set charges_begin iterator
    charges_begin = charges.begin();

    Base::M.resize(Base::tree.boxes());
    Base::L.resize(Base::tree.boxes());
    
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
          auto body2point = [](typename Octree<point_type>::Body b) { return b.point(); };

          auto p_begin = make_transform_iterator(box.body_begin(), body2point);
          auto p_end   = make_transform_iterator(box.body_end(),   body2point);
          auto c_begin = charges.begin() + box.body_begin()->index();

          printf("P2M: box: %d\n", (int)box.index());
          Base::K.P2M(p_begin, p_end, c_begin, box.center(), Base::M[idx]);

        } else {
          // If not leaf, make M2M calls

          // For all the children, M2M
          auto c_end = box.child_end();
          for (auto cit = box.child_begin(); cit != c_end; ++cit) {
            auto cbox = *cit;
            auto translation = box.center() - cbox.center();

            printf("M2M: %d to %d\n", cbox.index(), idx);
            Base::K.M2M(Base::M[cbox.index()], Base::M[idx], translation);
          }   
        }   
      }   
    }
  }

  //! Box-Box interactions
  void interactions(std::vector<result_type>& results)
  {
    // set reference to beginning of results
    results_begin = results.begin();

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
    // not needed for treecode -- quiet warning
    (void) results;
  }

  //! evaluator name
  std::string name()
  {
    return "Treecode";
  }
};

