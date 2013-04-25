#pragma once

#include <iostream>
#include <iomanip>

#include <map>


#include <sys/time.h>

/** Clock class, useful when timing code.
 */
struct Clock {
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



//! Timer and Trace logger
class Logger {

  struct EventData {
    Clock start_time;
    double total_time;
    int hit;
    EventData()
        : start_time(), total_time(0), hit(0) {
    }
    void start() {
      start_time.start();
    }
    double log(const Clock& end_time) {
      double elapsed = (end_time - start_time);
      total_time += elapsed;
      hit += 1;
      return elapsed;
    }
    friend std::ostream& operator<<(std::ostream& s, const EventData& e) {
      return s << e.hit << " (calls) * " << e.total_time/e.hit << " (sec/call) = "
               << e.total_time << " (secs)";
    }
  };

  std::map<std::string, EventData> data_;

 public:

  //! Start a clock for an event
  inline void start(const std::string& event) {
    data_[event].start();
  }

  //! Return the elasped time for given event
  double stop(const std::string& event, bool print_event = false) {
    Clock end_time;      // Stop the clock
    double elapsed = data_[event].log(end_time);

    if (print_event)
      std::cout << data_[event] << std::endl;

    return elapsed;
  }

  //! Erase entry in timer
  inline void clear(const std::string& event) {
    data_.erase(event);
  }

  //! Erase all events in timer
  inline void clear() {
    data_.clear();
  }

  // Get an event's data?
  //PublicEventData operator[](const std::string& event) { ... }

  //! Print all events and timing to an ostream
  friend std::ostream& operator<<(std::ostream& s, const Logger& log) {
    for (auto it : log.data_)
      s << std::setw(20) << std::left << it.first << " : "
        << it.second << std::endl;
    return s;
  }
};
