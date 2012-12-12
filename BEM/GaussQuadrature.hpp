#pragma once

#include <map>
#include <vector>
#include "Vec.hpp"

/** holder for Gaussian quadrature weights / points
 *  underlying implementation uses a map */

template <typename T=double>
struct GaussQuadrature
{
 private:
  // 1D raw data
  std::vector<T>        GQ_1D_points_n5 = { -0.9061798549, -0.5384693101, 0., 0.5384693101, 0.9061798549 };
  std::vector<T>        GQ_1D_weight_n5 = { 0.2369268850, 0.4786286704, 0.5688888888, 0.4786286704, 0.2369268850 };

  // 2D raw data
  std::vector<Vec<3,T>> GQ_2D_points_n1 = { Vec<3,T>(1./3,1./3,1./3) };
  std::vector<T>        GQ_2D_weight_n1 = { 1. };

  std::vector<Vec<3,T>> GQ_2D_points_n3 = { Vec<3,T>(0.5,0.5,0.),
                                         Vec<3,T>(0.,0.5,0.5),
                                         Vec<3,T>(0.5,0.,0.5) };
  std::vector<T>        GQ_2D_weight_n3 = { 1./3, 1./3, 1./3 };

  std::vector<Vec<3,T>> GQ_2D_points_n4 = { Vec<3,T>(1./3,1./3,1./3),
                                         Vec<3,T>(.6,.2,.2),
                                         Vec<3,T>(.2,.6,.2),
                                         Vec<3,T>(.2,.2,.6) };
  std::vector<T>        GQ_2D_weight_n4 = { -27./48, 25./48, 25./48, 25./48 };

  std::vector<Vec<3,T>> GQ_2D_points_n7 = { Vec<3,T>(1./3, 1./3, 1./3),
                                         Vec<3,T>(0.79742699, 0.10128651, 0.10128561),
                                         Vec<3,T>(0.10128651, 0.79742699, 0.10128651),
                                         Vec<3,T>(0.10128651, 0.10128651, 0.79742699),
                                         Vec<3,T>(0.05971587, 0.47014206, 0.47014206),
                                         Vec<3,T>(0.47014206, 0.05971587, 0.47014206),
                                         Vec<3,T>(0.47014206, 0.47014206, 0.13239415) };
  std::vector<T>        GQ_2D_weight_n7 = { 0.225, 0.12593918, 0.12593918, 0.12593918, 0.13239415, 0.13239415, 0.13239415 };

 public:
  // constructor
  GaussQuadrature() {
    // 1D: K = 5
    weights1D_[5] = GQ_1D_weight_n5;
    points1D_[5]  = GQ_1D_points_n5;
    // 2D: K = 1
    weights2D_[1] = GQ_2D_weight_n1;
    points2D_[1]  = GQ_2D_points_n1;
    // K = 3
    weights2D_[3] = GQ_2D_weight_n3;
    points2D_[3]  = GQ_2D_points_n3;
    // K = 4
    weights2D_[4] = GQ_2D_weight_n4;
    points2D_[4]  = GQ_2D_points_n4;
    // K = 7
    weights2D_[7] = GQ_2D_weight_n4;
    points2D_[7]  = GQ_2D_points_n4;
  }

  std::vector<T>& weights(unsigned k) {
    auto it = weights2D_.find(k);
    // error -- invalid # of quadrature points requested value
    if (it == weights2D_.end()) {
      std::cout << "[E] Invalid # of quadrature points requested - valid values are: "; //  << std::endl;
      for (auto mit : weights2D_) std::cout << mit.first << "  ";
      std::cout << std::endl;
      std::exit(0);
    }
    return it->second;
  };
  std::vector<Vec<3,T>>& points(unsigned k) {
    auto it = points2D_.find(k);
    // error -- invalid # of quadrature points requested value
    if (it == points2D_.end()) {
      std::cout << "[E] Invalid # of quadrature points requested - valid values are: "; //  << std::endl;
      for (auto mit : points2D_) std::cout << mit.first << "  ";
      std::cout << std::endl;
      std::exit(0);
    }
    return it->second;
  };

 private:
  std::map<unsigned,std::vector<T>> weights1D_;
  std::map<unsigned,std::vector<T>> points1D_;
  std::map<unsigned,std::vector<T>> weights2D_;
  std::map<unsigned,std::vector<Vec<3,T>>> points2D_;
};
