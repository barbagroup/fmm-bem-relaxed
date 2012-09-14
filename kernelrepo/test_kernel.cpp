#include "Kernel.hpp"
#include "HelmholtzKernel.hpp"
#include "Direct.hpp"

#include <iostream>

double getRandom() {
  return drand48();
}


template <typename T>
std::ostream& operator<<(std::ostream& s, const std::vector<T>& a) {
  for (T ai: a)
    s << ai << std::endl;
  return s;
}


template <typename Kernel>
void simple_test(const Kernel& K) {
  std::cout << K.name() << std::endl;

  int N = 10;

  std::vector<typename Kernel::point_type> p(N);
  std::vector<typename Kernel::charge_type> psi(N);

  for( int k = 0; k < N; ++k ) {
    p[k] = typename Kernel::point_type(getRandom(), getRandom(), getRandom());
    psi[k] = 1;
  }

  std::vector<typename Kernel::range_type> omega(N);

  Direct_Kernel_MatVec(K, p, psi, omega);

  std::cout << omega;
}


int main(int argc, char** argv) {

  simple_test(Kernel());

  std::cout << std::endl;

  simple_test(HelmholtzKernel(1));

  return 0;
}
