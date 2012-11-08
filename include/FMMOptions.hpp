#pragma once

/**
 * Storage class for all FMM options
 */


/** Class to define compile-time and run-time FMM options */
class FMMOptions
{
public:
  bool symmetric; // TODO: Call this Galerkin?  s == t

  //! Evaluation type
  enum EvalType {FMM, TREECODE};
  EvalType evaluator;



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
  DefaultMAC MAC_;
  unsigned NCRIT_;

  FMMOptions()
    : symmetric(false),
      evaluator(FMM),
      MAC_(DefaultMAC(0.5)),
      NCRIT_(126) {
  };

  void set_mac_theta(double theta) {
    MAC_ = DefaultMAC(theta);
  }

  DefaultMAC MAC() {
    return MAC_;
  }

  void set_max_per_box(unsigned ncrit) {
    NCRIT_ = ncrit;
  }

  unsigned max_per_box() const {
    return NCRIT_;
  }
};

