#pragma once

/*
 * Direct calculations for error checking purposes
 * Using K(x,y) interface
 */

class Direct
{
 public:
  template <typename Kernel, typename PointIter, typename ChargeIter, typename ResultIter>
  static void matvec(Kernel& K,
                     PointIter s_begin, PointIter s_end, ChargeIter c_begin,
                     PointIter t_begin, PointIter t_end, ResultIter r_begin)
  {
    typedef typename Kernel::result_type result_type;
    for ( ; t_begin!=t_end; ++t_begin, ++r_begin)
    {
      result_type r(0);

      auto s = s_begin;
      for ( ; s!=s_end; ++s, ++c_begin)
      {
        // r += K(*s,*t_begin)*(*c_begin);
        r += (*c_begin)*K(*s,*t_begin);
      }
      *r_begin = r;
    }
  }
};
