#pragma once

#include <iostream>
#include <iterator>
#include <algorithm>

#include <vector>
#include <string>

/** Helper class for type debugging in a static_asset
 * Usage:
 * static_assert(AssertValue<
 *                     std::is_same<Complicated1, Complicated2>
 *               >::value,
 *               "Something horrible happened...");
 * Note the ::value is left off the template parameter to be used in AssertValue
 * AssertValue will fail instantiation and print out the full type of
 * the template parameter.
 */
template <typename Assertion>
struct AssertValue {
  static_assert(Assertion::value,
                "Assertion failed <see below for more information>");
  static bool const value = Assertion::value;
};



/** Bucket sort using pigeonhole sorting
 *
 * @param[in] first,last Iterator pair to the sequence to be bucketed
 * @param[in] num_buckets The number of buckets to be used
 * @param[in] map Functor that maps elements in [first,last) to [0,num_buckets)
 * @returns vector of iterators:
 *          bucket_off[i],bucket_off[i+1] are the [first,last) of the ith bucket.
 * @pre For all elements i in [first,last), 0 <= map(i) < num_buckets.
 * @post For all elements i and j such that first <= i < j < last,
 *       then 0 <= map(i) <= map(j) < num_buckets.
 */
template <typename Iterator, typename BucketMap>
std::vector<Iterator> bucket_sort(Iterator first, Iterator last,
                                  int num_buckets, BucketMap map) {
  typedef typename std::iterator_traits<Iterator>::value_type value_type;
  std::vector<std::vector<value_type>> buckets(num_buckets);

  // Push each element into a bucket
  std::for_each(first, last, [&buckets, &map] (const value_type& v) {
      buckets[map(v)].push_back(v);
    });

  // Copy the buckets back to the range and keep each offset iterator
  std::vector<Iterator> bucket_off(num_buckets+1);
  auto offset = bucket_off.begin();
  (*offset++) = first;
  std::accumulate(buckets.begin(), buckets.end(), first,
                  [&offset](Iterator out, const std::vector<value_type>& bucket)
                  { return (*offset++) =
                        std::copy(bucket.begin(), bucket.end(), out); });

  return bucket_off;
}


