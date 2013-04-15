#include "Vec.hpp"

int main() {
  Vec<3,double> v1 = Vec<3,double>(0,1,2);
  Vec<3,double> v2 = Vec<3,double>(3,4,5);

  Vec<3,double> v3 = v1 * (v2 + v2);

  std::cout << v3 << std::endl;
  std::cout << norm(v3) << std::endl;

  return 0;
}
