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
 * @param[in] begin,end Iterator pair to the sequence to be bucketed
 * @param[in] num_buckets The number of buckets to be used
 * @param[in] map Functor that maps elements in [begin,end) to [0,num_buckets)
 * @returns Vector of iterators:
 *          bucket_off[i],bucket_off[i+1] are the [begin,end) of ith bucket
 * @post For all elements i and j such that begin <= i < j < end,
 *       then 0 <= map(i) <= map(j) < num_buckets
 */
template <typename Iterator, typename BucketMap>
std::vector<typename std::vector<typename Iterator::value_type>::iterator>
bucket_sort(Iterator begin, Iterator end, int num_buckets, BucketMap map) {
  typedef typename std::iterator_traits<Iterator>::value_type value_type;
  std::vector<std::vector<value_type>> buckets(num_buckets);

  // Push each element into a bucket
  std::for_each(begin, end, [&buckets, &map] (value_type& v) {
      buckets[map(v)].push_back(v);
    });

  // Copy the buckets back to the range and keep each iterator
  std::vector<typename std::vector<value_type>::iterator> bucket_off(num_buckets+1);
  auto off_iter = bucket_off.begin();
  (*off_iter++) = begin;
  std::accumulate(buckets.begin(), buckets.end(), begin,
                  [&off_iter](typename std::vector<value_type>::iterator out,
                              std::vector<value_type>& bucket)
                  { return (*off_iter++) =
                        std::copy(bucket.begin(), bucket.end(), out); });

  return bucket_off;
}
