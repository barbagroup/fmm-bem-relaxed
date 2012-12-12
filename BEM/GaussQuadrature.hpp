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
  // raw data
  std::vector<Vec<3,T>> GQ_points_n1 = { Vec<3,T>(1./3,1./3,1./3) };
  std::vector<T>        GQ_weight_n1 = { 1. };

  std::vector<Vec<3,T>> GQ_points_n3 = { Vec<3,T>(0.5,0.5,0.),
                                             Vec<3,T>(0.,0.5,0.5),
                                             Vec<3,T>(0.5,0.,0.5) };
  std::vector<T>        GQ_weight_n3 = { 1./3, 1./3, 1./3 };

  std::vector<Vec<3,T>> GQ_points_n4 = { Vec<3,T>(1./3,1./3,1./3),
                                             Vec<3,T>(.6,.2,.2),
                                             Vec<3,T>(.2,.6,.2),
                                             Vec<3,T>(.2,.2,.6) };
  std::vector<T>        GQ_weight_n4 = { -27./48, 25./48, 25./48, 25./48 };

 public:
  // constructor
  GaussQuadrature() {
    weightsMap_[1] = GQ_weight_n1;
    pointsMap_[1]  = GQ_points_n1;
    weightsMap_[3] = GQ_weight_n3;
    pointsMap_[3]  = GQ_points_n3;
    pointsMap_[4]  = GQ_points_n4;
    pointsMap_[4]  = GQ_points_n4;
  }

  std::vector<T>& weights(unsigned k) {
    auto it = weightsMap_.find(k);
    // error -- invalid # of quadrature points requested value
    if (it == weightsMap_.end()) {
      std::cout << "[E] Invalid # of quadrature points requested - valid values are:" << std::endl;
      for (auto mit : weightsMap_) std::cout << mit.first << std::endl;
      std::exit(0);
    }
    return it->second;
  };
  std::vector<Vec<3,T>>& points(unsigned k) {
    auto it = pointsMap_.find(k);
    // error -- invalid # of quadrature points requested value
    if (it == pointsMap_.end()) {
      std::cout << "[E] Invalid # of quadrature points requested - valid values are:" << std::endl;
      for (auto mit : pointsMap_) std::cout << mit.first << std::endl;
      std::exit(0);
    }
    return it->second;
  };

 private:
  std::map<unsigned,std::vector<T>> weightsMap_;
  std::map<unsigned,std::vector<Vec<3,T>>> pointsMap_;
};
