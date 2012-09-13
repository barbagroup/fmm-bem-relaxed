#include "Kernel.hpp"
#include "Direct.hpp"

#include <iostream>

template <typename T>
std::ostream& operator<<(std::ostream& s, const std::vector<T> a) {
  for (T ai: a)
    s << ai << std::endl;
  return s;
}

double getRandom() { return drand48(); }

int main(int argc, char** argv) {
  Kernel K;

  int N = 10;

  std::vector<Kernel::point_type> p(N);
  std::vector<Kernel::charge_type> psi(N);

  for( int k = 0; k < N; ++k ) {
    p[k] = Kernel::point_type(getRandom(), getRandom(), getRandom());
    psi[k] = 1;
  }

  std::vector<Kernel::range_type> omega(N);

  Direct(K, p, psi, omega);

  std::cout << omega;

  return 0;
}
