#include "include/Vec.hpp"
#include <iostream>
#include <vector>

struct MyVec {
  std::vector<float> a;

  MyVec()
      : a(4) {
  }

  float& operator[](int k) {
    return a[k];
  }
  const float& operator[](int k) const {
    return a[k];
  }
};



int main()
{
  std::cout << Vec<3,double>(0,0,0) << "\n";

  std::cout << Vec<3,double[3]>(0,0,0) << "\n";

  std::cout << Vec<4,MyVec>(0, 0, 0, 0) << "\n";

  return 0;
}
