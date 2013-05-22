#pragma once

#include "GaussQuadrature.hpp"

/** Singleton pattern to hold BEM config options
 * see: http://stackoverflow.com/questions/2496918/singleton-pattern-in-c */
class BEMConfig {

 public:
  static void Init();
  static BEMConfig *Instance();
  static void Destroy();

  void setK(unsigned k) {
    K = k;
  }

  unsigned getK() {
    return K;
  }

  std::vector<Vec<3,double>>& GaussPoints() {
    return GQ.points(K);
  }
  std::vector<double>& GaussWeights() {
    return GQ.weights(K);
  }

  std::vector<Vec<3,double>>& GaussPoints(int k) {
    return GQ.points(k);
  }
  std::vector<double>& GaussWeights(int k) {
    return GQ.weights(k);
  }

 private:
  static BEMConfig* MyInstance(BEMConfig* pConfig);
  BEMConfig() { };   // static variables
  ~BEMConfig() {};
  BEMConfig(const BEMConfig& config);
  BEMConfig operator=(const BEMConfig& config);

  // actual data
  unsigned K;
  GaussQuadrature<double> GQ;
};

inline void BEMConfig::Init()
{
  BEMConfig* ptr = new BEMConfig();
  MyInstance(ptr);
  atexit(&Destroy);
}

inline BEMConfig* BEMConfig::Instance()
{
  return BEMConfig::MyInstance(NULL);
}

inline BEMConfig* BEMConfig::MyInstance(BEMConfig* ptr)
{
  static BEMConfig* myInstance = NULL;
  if (ptr) myInstance = ptr;
  return myInstance;
}

inline void BEMConfig::Destroy()
{
  BEMConfig* pConfig = MyInstance(NULL);
  if (pConfig) {
    delete pConfig;
    pConfig = 0;
  }
}

