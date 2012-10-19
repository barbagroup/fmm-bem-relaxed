#pragma once

/*
 * Direct calculations for error checking purposes
 * Using K(x,y) interface
 */

class Direct
{
 public:
  // One-sided P2P
  /** Non-symmetric P2P
   * r_j += sum_i K(s_j, t_i) * c_i
   *
   */
  template <typename Kernel, typename PointIter, typename ChargeIter, typename ResultIter>
  static void matvec(Kernel& K,
                     PointIter s_begin, PointIter s_end, ChargeIter c_begin,
                     PointIter t_begin, PointIter t_end, ResultIter r_begin)
  {
    typedef typename Kernel::result_type result_type;
    for ( ; t_begin!=t_end; ++t_begin, ++r_begin) {
      result_type r(0);

      for (auto s = s_begin; s!=s_end; ++s, ++c_begin) {
        // r += K(*s,*t_begin)*(*c_begin);
        r += (*c_begin) * K(*t_begin, *s);
      }
      *r_begin = r;
    }
  }

  // Two-sided P2P
  // TODO: Need to know if the kernel is symmetric or anti-symmetric
  template <typename Kernel, typename PointIter, typename ChargeIter, typename ResultIter>
  static void matvec(Kernel& K,
                     PointIter p1_begin, PointIter p1_end, ChargeIter c1_begin,
                     PointIter p2_begin, PointIter p2_end, ChargeIter c2_begin,
                     ResultIter r1_begin, ResultIter r2_begin)
  {
    typedef typename Kernel::result_type result_type;
    for ( ; p1_begin != p1_end; ++p1_begin, ++c1_begin, ++r1_begin) {
      result_type r(0);

      auto p2 = p2_begin;
      auto c2 = c2_begin;
      auto r2 = r2_begin;
      for ( ; p2 != p2_end; ++p2, ++c2, ++r2) {
        auto kts = K(*p1_begin, *p2);

        if (Kernel::is_symmetric) {

        } else {

        }
      }
    }
  }
};
