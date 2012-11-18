#pragma once

#include "Util.hpp"

#include <iostream>
#include <iomanip>

#include <map>

//! Timer and Trace logger
class Logger {
  std::map<std::string, Clock> event_start_;
  std::map<std::string, double> event_time_;

  //! Print event name and time to output stream
  inline static void print(std::ostream& s,
                           const std::string& event, double time) {
    s << std::setw(20) << std::left
      << event << " : " << time << std::endl;
  }

 public:

  //! Print event name and time to output stream
  inline void print(std::ostream& s,
                    const std::string& event) {
    print(s, event, event_time_[event]);
  }

  //! Start a clock for an event
  inline void start(const std::string& event) {
    event_start_[event].start();
  }

  //! Return the elasped time for given event
  double stop(const std::string& event, bool print_event = false) {
    Clock event_end;
    double event_elapsed = event_end - event_start_[event];

    // Accumulate event times in the log
    event_time_[event] += event_elapsed;

    if (print_event) print(std::cout, event);
    return event_elapsed;
  }

  //! Erase entry in timer
  inline void clear(const std::string& event) {
    event_time_.erase(event);
  }

  //! Erase all events in timer
  inline void clear() {
    event_time_.clear();
  }

  //! Print all events and timing to an ostream
  friend std::ostream& operator<<(std::ostream& s, const Logger& log) {
    for (auto it : log.event_time_)
      print(s, it.first, it.second);
    return s;
  }
};

