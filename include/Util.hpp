#pragma once

#include <iostream>
#include <iterator>
#include <algorithm>

#include <vector>
#include <string>

#include <sys/time.h>

/** Clock class, useful when timing code.
 */
class Clock {
 public:
  /** Construct a Clock and start timing. */
  Clock() {
    start();
  }
  /** Start the clock. */
  inline void start() {
    time_ = now();
  }
  /** Return the seconds elapsed since the last start. */
  inline double elapsed() const {
    timeval tv = now();
    timersub(&tv, &time_, &tv);
    return seconds(tv);
  }
  /** Return the seconds difference between the Clocks */
  inline double operator-(const Clock& clock) const {
    timeval tv;
    timersub(&time_, &clock.time_, &tv);
    return seconds(tv);
  }
 private:
  timeval time_;
  inline static timeval now() {
    timeval tv;
    gettimeofday(&tv, nullptr);
    return tv;
  }
  inline static double seconds(const timeval& tv) {
    return tv.tv_sec + 1e-6 * tv.tv_usec;
  }
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
