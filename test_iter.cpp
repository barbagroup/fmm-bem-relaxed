

#include <iostream>
#include <vector>

#include "include/TransformIterator.hpp"

struct MyOp {
  int operator()(const int& a) const {
    return a*a;
  }
};

template <typename IT>
transform_iterator<IT,OP>


int main(int argc, char** argv)
{
  std::vector<int> vec(10);

  for (int k = 0; k < vec.size(); ++k) {
    vec[k] = k;
  }

  for (auto v : vec)
    std::cout << v << "\n";

  std::cout << "\n";

  auto begin = make_transform_iterator(vec.begin(), MyOp());
  auto end = make_transform_iterator(vec.end(), MyOp());

  for ( ; begin != end; ++begin)
    std::cout << *begin << "\n";

  std::cout << "\n";

  auto lambda = [](int& a) { return a*a; };
  auto begin2 = make_transform_iterator(vec.begin(), lambda);
  auto end2 = make_transform_iterator(vec.end(), lambda);

  for ( ; begin2 != end2; ++begin2)
    std::cout << *begin2 << "\n";

  return 0;
}
