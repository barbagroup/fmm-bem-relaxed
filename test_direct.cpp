#include "include/Direct.hpp"
#include "include/Vec.hpp"

#include <vector>
#include <iostream>

struct TempKernel {
  //! Multipole expansion type
  typedef unsigned multipole_type;
  //! Local expansion type
  typedef unsigned local_type;

  //! The dimension of the Kernel
  static constexpr unsigned dimension = 3;
  //! Point type
  // typedef vec<dimension,real> point_type;
  typedef Vec<dimension,double> point_type;
  //! Charge type
  typedef double charge_type;
  //! The return type of a kernel evaluation
  typedef unsigned kernel_value_type;
  //! The product of the kernel_value_type and the charge_type
  typedef double result_type;

  kernel_value_type operator()(const point_type& t,
                               const point_type& s) const {
    (void) t;
    (void) s;

    std::cout << "In op()\n";

    return kernel_value_type(1);
  }

  void P2P(int a) const { std::cout << "In P2P(int)\n"; }

  kernel_value_type transpose(const kernel_value_type& kst) const {
    std::cout << "In Transpose\n";
    return kernel_value_type(1);
  }

#if 1
  template <typename PointIter, typename ChargeIter, typename ResultIter>
  void P2P(PointIter s_begin, PointIter s_end, ChargeIter c_begin,
           PointIter t_begin, PointIter t_end, ResultIter r_begin) const {
    (void) s_begin;
    (void) s_end;
    (void) c_begin;
    (void) t_begin;
    (void) t_end;
    (void) r_begin;

    std::cout << "In P2P\n";
  }
#endif
};




int main() {
  typedef TempKernel kernel_type;
  typedef kernel_type::point_type point_type;
  typedef kernel_type::charge_type charge_type;
  typedef kernel_type::result_type result_type;
  typedef kernel_type::multipole_type multipole_type;
  typedef kernel_type::local_type local_type;

  kernel_type K;

  // init source
  std::vector<point_type> points(1);
  points[0] = point_type(0,0,0);

  // init charge
  std::vector<charge_type> charges(1);
  charges[0] = 1.;

  // init target
  std::vector<point_type> target(1);
  target[0] = point_type(1,1,1);

  // init results vectors for exact
  std::vector<result_type> exact(target.size());

  // test direct
  Direct::matvec(K, points, charges, exact);

  std::cout << exact[0] << "\n";
}







