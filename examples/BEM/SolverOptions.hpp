#pragma once

/**
 * Central place to store Solver config options
 */

#include <cstdio>
#include <vector>
#include <cmath>

struct SolverOptions
{
  double residual;
  int max_iters, restart;
  unsigned max_p, p_min;
  bool variable_p;

  enum relaxation_type { SIMONCINI, BOURAS };

  relaxation_type relax_type;

  SolverOptions(double r, int m_iters, unsigned p) : residual(r), max_iters(m_iters), restart(50), max_p(p), variable_p(false), relax_type(BOURAS) {};
  SolverOptions() : residual(1e-5), max_iters(500), restart(500), max_p(16), p_min(5), variable_p(true), relax_type(BOURAS) {};

  unsigned predict_p(double eps) const {
    // if no relaxation, return default p
    if (!this->variable_p) return max_p;
    // relaxation from Bouras & Fraysse
    if (this->relax_type == BOURAS) {
      double alpha = 1. / std::min(eps, 1.);
      double nu = std::min(alpha*this->residual,1.);
      // predict p for Spherical Laplace kernel -- abstract out
      return std::min((unsigned)ceil(-log(nu)),max_p);
    } else if (this->relax_type == SIMONCINI) {
      return std::min((unsigned)ceil(-log(eps)),max_p);
    }
    return max_p;
  }
};

