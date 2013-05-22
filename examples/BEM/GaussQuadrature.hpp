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
    // K = 13
    {
      double a, b1, b2, c1, c2, d1, d2, d3;
      a = 1/3.;
      b1 = 0.479308067841920; b2 = 0.260345966079040;
      c1 = 0.869739794195568; c2 = 0.065130102902216;
      d1 = 0.048690315425316; d2 = 0.312865496004874; d3 = 0.638444188569810;

      double wa, wb, wc, wd;
      wa = -0.149570044467682;
      wb = 0.175615257433208;
      wc = 0.053347235608838;
      wd = 0.077113760890257;

      std::vector<Vec<3,T>> GQ_2D_points_n13 = { Vec<3,T>(a,a,a),
                                                 Vec<3,T>(b1,b2,b2), Vec<3,T>(b2,b1,b2), Vec<3,T>(b2,b2,b1),
                                                 Vec<3,T>(c1,c2,c2), Vec<3,T>(c2,c1,c2), Vec<3,T>(c2,c2,c1),
                                                 Vec<3,T>(d1,d2,d3), Vec<3,T>(d1,d3,d2), Vec<3,T>(d2,d1,d3), Vec<3,T>(d2,d3,d1), Vec<3,T>(d3,d1,d2), Vec<3,T>(d3,d2,d1) };
      std::vector<T> GQ_2D_weights_n13 = { wa,
                                           wb, wb, wb,
                                           wc, wc, wc,
                                           wd, wd, wd, wd, wd, wd };

      points2D_[13] = GQ_2D_points_n13;
      weights2D_[13] = GQ_2D_weights_n13;
    }
    // K = 17
    {
      double a, b1, b2, c1, c2, d1, d2, e1, e2, e3;
      a = 1/3.;
      b1 = 0.081414823414554; b2 = 0.459292588292723;
      c1 = 0.658861384496480; c2 = 0.170569307751760;
      d1 = 0.898905543365938; d2 = 0.050547228317031;
      e1 = 0.008394777409958; e2 = 0.263112829634638; e3 = 0.728492392955404;

      double wa, wb, wc, wd, we;
      wa = 0.144315607677787;
      wb = 0.095091634267285;
      wc = 0.103217370534718;
      wd = 0.032458497623198;
      we = 0.027230314174435;

      std::vector<Vec<3,T>> GQ_2D_points_n17 = { Vec<3,T>(a,a,a),
                                                 Vec<3,T>(b1,b2,b2), Vec<3,T>(b2,b1,b2), Vec<3,T>(b2,b2,b1),
                                                 Vec<3,T>(c1,c2,c2), Vec<3,T>(c2,c1,c2), Vec<3,T>(c2,c2,c1),
                                                 Vec<3,T>(d1,d2,d2), Vec<3,T>(d2,d1,d2), Vec<3,T>(d2,d2,d1),
                                                 Vec<3,T>(e1,e2,e3), Vec<3,T>(e1,e3,e2), Vec<3,T>(e2,e1,e3), Vec<3,T>(e2,e3,e1), Vec<3,T>(e3,e1,e2), Vec<3,T>(e3,e2,e1) };

      std::vector<T> GQ_2D_weights_n17 = { wa,
                                           wb, wb, wb,
                                           wc, wc, wc,
                                           wd, wd, wd,
                                           we, we, we, we, we, we };

      points2D_[17] = GQ_2D_points_n17;
      weights2D_[17] = GQ_2D_weights_n17;
    }
    // K = 19
    {
      double a, b1, b2, c1, c2, d1, d2, e1, e2, f1, f2, f3;
      a = 1/3.;
      b1 = 0.020634961602525; b2 = 0.489682519198738;
      c1 = 0.125820817014127; c2 = 0.437089591492937;
      d1 = 0.623592928761935; d2 = 0.188203535619033;
      e1 = 0.910540973211095; e2 = 0.044729513394453;
      f1 = 0.036838412054736; f2 = 0.221962989160766; f3 = 0.741198598784498;

      double wa, wb, wc, wd, we, wf;
      wa = 0.097135796282799;
      wb = 0.031334700227139;
      wc = 0.077827541004774;
      wd = 0.079647738927210;
      we = 0.025577675658698;
      wf = 0.043283539377289;

      std::vector<Vec<3,T>> GQ_2D_points_n19 = { Vec<3,T>(a,a,a),
                                                 Vec<3,T>(b1,b2,b2), Vec<3,T>(b2,b1,b2), Vec<3,T>(b2,b2,b1),
                                                 Vec<3,T>(c1,c2,c2), Vec<3,T>(c2,c1,c2), Vec<3,T>(c2,c2,c1),
                                                 Vec<3,T>(d1,d2,d2), Vec<3,T>(d2,d1,d2), Vec<3,T>(d2,d2,d1),
                                                 Vec<3,T>(e1,e2,e2), Vec<3,T>(e2,e1,e2), Vec<3,T>(e2,e2,e1),
                                                 Vec<3,T>(f1,f2,f3), Vec<3,T>(f1,f3,f2), Vec<3,T>(f2,f1,f3), Vec<3,T>(f2,f3,f1), Vec<3,T>(f3,f1,f2), Vec<3,T>(f3,f2,f1)};
      std::vector<T> GQ_2D_weights_n19 = { wa,
                                           wb, wb, wb,
                                           wc, wc, wc,
                                           wd, wd, wd,
                                           we, we, we,
                                           wf, wf, wf, wf, wf, wf };

      points2D_[19] = GQ_2D_points_n19;
      weights2D_[19] = GQ_2D_weights_n19;
    }
    // K = 25
    {
      double a, b1, b2, c1, c2, d1, d2, d3, e1, e2, e3, f1, f2, f3;
      a  = 1/3.;
      b1 = 0.028844733232685; b2 = 0.485577633383657;
      c1 = 0.781036849029926; c2 = 0.109481575485037;
      d1 = 0.141707219414880; d2 = 0.307939838764121; d3 = 0.550352941820999;
      e1 = 0.025003534762686; e2 = 0.246672560639903; e3 = 0.728323904597411;
      f1 = 0.009540815400299; f2 = 0.066803251012200; f3 = 0.923655933587500;

      double wa, wb, wc, wd, we, wf;
      wa = 0.090817990382754;
      wb = 0.036725957756467;
      wc = 0.045321059435528;
      wd = 0.072757916845420;
      we = 0.028327242531057;
      wf = 0.009421666963733;

      std::vector<Vec<3,T>> GQ_2D_points_n25 = { Vec<3,T>(a,a,a),
                                                 Vec<3,T>(b1,b2,b2), Vec<3,T>(b2,b1,b2), Vec<3,T>(b2,b2,b1),
                                                 Vec<3,T>(c1,c2,c2), Vec<3,T>(c2,c1,c2), Vec<3,T>(c2,c2,c1),
                                                 Vec<3,T>(d1,d2,d3), Vec<3,T>(d1,d3,d2), Vec<3,T>(d2,d1,d3), Vec<3,T>(d2,d3,d1), Vec<3,T>(d3,d1,d2), Vec<3,T>(d3,d2,d1),
                                                 Vec<3,T>(e1,e2,e3), Vec<3,T>(e1,e3,e2), Vec<3,T>(e2,e1,e3), Vec<3,T>(e2,e3,e1), Vec<3,T>(e3,e1,e2), Vec<3,T>(e3,e2,e1),
                                                 Vec<3,T>(f1,f2,f3), Vec<3,T>(f1,f3,f2), Vec<3,T>(f2,f1,f3), Vec<3,T>(f2,f3,f1), Vec<3,T>(f3,f1,f2), Vec<3,T>(f3,f2,f1) };

      std::vector<T> GQ_2D_weights_n25 = { wa,
                                           wb, wb, wb,
                                           wc, wc, wc,
                                           wd, wd, wd, wd, wd, wd,
                                           we, we, we, we, we, we,
                                           wf, wf, wf, wf, wf, wf };

      points2D_[25] = GQ_2D_points_n25;
      weights2D_[25] = GQ_2D_weights_n25;
    }
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
