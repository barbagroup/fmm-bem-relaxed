#pragma once

/** Base and derived precondtioners for BEM problems */

namespace Preconditioners {

/** Identity preconditioner -- no change */
class Identity
{
 public:
  template <typename VecType>
  void operator()(const VecType& x, VecType& y) const {
    auto yit = y.begin();
    for (auto it=x.begin(); it!=x.end(); ++it, ++yit) *yit = *it;
  };
};

/** Diagonal PC -- divide by Panel self-interactions */
template <typename ValueType>
class Diagonal
{
 public:
  template <typename Kernel, typename SourceIter>
  Diagonal(Kernel& K, SourceIter SourceBegin, SourceIter SourceEnd) {
    Reciprocals.resize(SourceEnd-SourceBegin);
    auto Rit = Reciprocals.begin();
    for ( ; SourceBegin!=SourceEnd; ++SourceBegin, ++Rit) {
      // generate the reciprocal of the self-interaction
      *Rit = 1./K(*SourceBegin, *SourceBegin);
    }
  }

  template <typename VecType>
  void operator()(const VecType& x, VecType& y) const {
    auto yit = y.begin();
    auto Rit = Reciprocals.begin();
    for (auto it=x.begin(); it!=x.end(); ++it, ++yit, ++Rit) *yit = (*Rit) * (*it);
  }
 private:
  std::vector<ValueType> Reciprocals;
};

}; // end namespace Preconditioners
