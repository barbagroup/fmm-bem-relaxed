#pragma once

/**
 * Storage class for all FMM options
 */


/** Class to define compile-time and run-time FMM options */
class FMMOptions
{
public:
  typedef enum {TOPDOWN, BOTTOMUP} TreeType;
  typedef enum {FMM, TREECODE} EvaluatorType;

  bool symmetric;
  TreeType tree;
  EvaluatorType evaluator;
  double THETA;
  unsigned NCRIT;

  struct DefaultMAC {
    double theta_;
    DefaultMAC(double theta) : theta_(theta) {}

    template <typename BOX>
    bool operator()(const BOX& b1, const BOX& b2) const {
      double r0_norm = norm(b1.center() - b2.center());
      return r0_norm * theta_ > b1.radius() + b2.radius();
    }
  };

  // TODO: Generalize type
  DefaultMAC MAC;

  FMMOptions()
    : symmetric(false),
      tree(TOPDOWN),
      evaluator(FMM),
      NCRIT(126),
      MAC(DefaultMAC(0.5)) {
  };

  void set_theta(double theta) {
    MAC = DefaultMAC(theta);
  }
};

