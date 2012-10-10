#pragma once

/**
 * Base class for all evaluators
 */

#include <EvaluatorFMM.hpp>
#include <EvaluatorTreecode.hpp>

template <class Tree, class Kernel>
class EvaluatorBase
{
 public:
  //! Kernel type
  typedef Kernel kernel_type;
  //! Tree type
  typedef Tree tree_type;
  //! Point type
  typedef typename Kernel::point_type point_type;
  //! Multipole expansion type
  typedef typename Kernel::multipole_type multipole_type;
  //! Local expansion type
  typedef typename Kernel::local_type local_type;
  //! Kernel source type
  typedef typename Kernel::charge_type charge_type;
  //! Kernel result type
  typedef typename Kernel::result_type result_type;
  //! Box type used in the Tree
  typedef typename tree_type::box_type box_type;
  //! Body type used in the Tree
  typedef typename tree_type::body_type body_type;

 protected:
  //! Handle to Tree
  Tree& tree;
  //! Handle to the Kernel
  Kernel& K;

  //! Multipole expansions corresponding to Box indices in Octree
  std::vector<multipole_type> M;
  //! Local expansions corresponding to Box indices in Octree
  std::vector<local_type> L;

  // THETA for multipole acceptance criteria
  double THETA;
 private:
  // TODO: Still don't like this...
  typename std::vector<charge_type>::const_iterator charges_begin;
  typename std::vector<result_type>::iterator results_begin;

 public:
  EvaluatorBase(Tree& t, Kernel& k, double theta)
      : tree(t), K(k), THETA(theta) {
  };
  virtual ~EvaluatorBase() {};

  void execute(const std::vector<charge_type>& charges,
               std::vector<result_type>& results) {
    M.resize(tree.boxes());
    L.resize(tree.boxes());

    charges_begin = charges.begin();
    results_begin = results.begin();

    upward(charges);
    interactions(results);
    downward(results);
  }

  //! Upward sweep
  virtual void upward(const std::vector<charge_type>& charges) = 0;
  //! 'Interaction' stage
  virtual void interactions(std::vector<result_type>& results) = 0;
  //! Downward sweep
  virtual void downward(std::vector<result_type>& results) = 0;

  //! set the value of THETA
  void setTheta(double th) {
    THETA = th;
  }
  //! get the value of THETA
  double getTheta() {
    return THETA;
  }

  //! Name of the evaluator
  virtual std::string name() = 0;

  /** One-sided P2P!!
   */
  void evalP2P(const box_type& b1, const box_type& b2) {
    // Point iters
    auto body2point = [](const body_type& b) { return b.point(); };
    auto p1_begin = make_transform_iterator(b1.body_begin(), body2point);
    auto p1_end   = make_transform_iterator(b1.body_end(),   body2point);
    auto p2_begin = make_transform_iterator(b2.body_begin(), body2point);
    auto p2_end   = make_transform_iterator(b2.body_end(),   body2point);

    // Charge iters
    auto c1_begin = charges_begin + b1.body_begin()->index();

    // Result iters
    auto r2_begin = results_begin + b2.body_begin()->index();

    printf("P2P: %d to %d\n", b1.index(), b2.index());

    K.P2P(p1_begin, p1_end, c1_begin,
          p2_begin, p2_end,
          r2_begin);
  }

  void evalP2M(const box_type& b) {
    // Source point iters
    auto body2point = [](const body_type& b) { return b.point(); };
    auto p_begin = make_transform_iterator(b.body_begin(), body2point);
    auto p_end   = make_transform_iterator(b.body_end(),   body2point);
    auto c_begin = charges_begin + b.body_begin()->index();

    printf("P2M: box: %d\n", b.index());

    K.P2M(p_begin, p_end, c_begin,
          b.center(),
          M[b.index()]);
  }

  void evalM2M(const box_type& box) {
    unsigned idx = box.index();

    // For all the children, M2M
    auto c_end = box.child_end();
    for (auto cit = box.child_begin(); cit != c_end; ++cit) {
      auto cbox = *cit;
      auto translation = box.center() - cbox.center();

      printf("M2M: %d to %d\n", cbox.index(), idx);
      K.M2M(M[cbox.index()],
            M[idx],
            translation);
    }
  }

  void evalM2P(const box_type& b1, const box_type& b2) {
    // Target point iters
    auto body2point = [](const body_type& b) { return b.point(); };
    auto t_begin = make_transform_iterator(b2.body_begin(), body2point);
    auto t_end   = make_transform_iterator(b2.body_end(), body2point);

    // Target result iters
    auto r_begin = results_begin + b2.body_begin()->index();

    printf("M2P: %d to %d\n", b1.index(), b2.index());

    K.M2P(M[b1.index()], b1.center(),
          t_begin, t_end,
          r_begin);
  }

  void evalM2L(const box_type& b1, const box_type& b2) {
    auto translation = b2.center() - b1.center();

    printf("M2L: %d to %d\n", b2.index(), b1.index());

    K.M2L(M[b1.index()],
          L[b2.index()],
          translation);
  }

  void evalL2L(const box_type& box) {
    // For all the children, L2L
    auto c_end = box.child_end();
    for (auto cit = box.child_begin(); cit != c_end; ++cit) {
      auto cbox = *cit;
      auto translation = cbox.center() - box.center();

      printf("L2L: %d to %d\n", box.index(), cbox.index());
      K.L2L(L[box.index()],
            L[cbox.index()],
            translation);
    }
  }

  void evalL2P(const box_type& box) {
    // For all the bodies, L2P
    auto body2point = [](const body_type& b) { return b.point(); };
    auto t_begin = make_transform_iterator(box.body_begin(), body2point);
    auto t_end   = make_transform_iterator(box.body_end(), body2point);
    auto r_begin = results_begin + box.body_begin()->index();

    printf("L2P: %d\n", box.index());
    K.L2P(L[box.index()], box.center(),
          t_begin, t_end,
          r_begin);
  }
};
