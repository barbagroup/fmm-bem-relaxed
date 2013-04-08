#pragma once
/** @file P2P.hpp
 * @brief Dispatch methods for the P2P stage
 *
 */

#include "Direct.hpp"

struct P2P
{
  //////////////////////////////////////
  /////// Static Dispatchers ///////////
  //////////////////////////////////////

  struct ONE_SIDED {};
  struct TWO_SIDED {};

  /** One sided P2P
   */
  template <typename Kernel, typename Context>
  inline static void eval(const Kernel& K,
                          Context& bc,
                          const typename Context::box_type& source,
                          const typename Context::box_type& target,
                          const ONE_SIDED&)
  {
#ifdef DEBUG
    printf("P2P: %d to %d\n", source.index(), target.index());
#endif

    Direct::matvec(K,
                   bc.source_begin(source), bc.source_end(source),
                   bc.charge_begin(source),
                   bc.target_begin(target), bc.target_end(target),
                   bc.result_begin(target));
  }

  /** Two sided P2P
   */
  template <typename Kernel, typename Context>
  inline static void eval(const Kernel& K,
                          Context& bc,
                          const typename Context::box_type& source,
                          const typename Context::box_type& target,
                          const TWO_SIDED&)
  {
#ifdef DEBUG
    printf("P2P: %d to %d\n", source.index(), target.index());
    printf("P2P: %d to %d\n", target.index(), source.index());
#endif

    Direct::matvec(K,
                   bc.source_begin(source), bc.source_end(source),
                   bc.charge_begin(source), bc.result_begin(source),
                   bc.target_begin(target), bc.target_end(target),
                   bc.charge_begin(target), bc.result_begin(target));
  }

  /** Two sided, self P2P
   */
  template <typename Kernel, typename Context>
  inline static void eval(const Kernel& K,
                          Context& bc,
                          const typename Context::box_type& source)
  {
#ifdef DEBUG
    printf("P2P: %d to %d\n", source.index(), source.index());
#endif

    Direct::matvec(K,
                   bc.source_begin(source), bc.source_end(source),
                   bc.charge_begin(source), bc.result_begin(source),
                   bc.target_begin(source), bc.target_end(source),
                   bc.charge_begin(source), bc.result_begin(source));
  }
};



struct P2P_Batch {
  /** Construct a P2P context on the bodies within a single box */
  template <typename Context>
  P2P_Batch(Context& bc,
      const typename Context::box_type& b) {
    // Do nothing
  }

  /** Construct a P2P context on the sources and targets within two boxes */
  template <typename Context>
  P2P_Batch(Context& bc,
      const typename Context::box_type& source,
      const typename Context::box_type& target) {
    // Do nothing
  }

  inline void compute() {
    // Do nothing
  }
};
