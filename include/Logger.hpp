#pragma once

#include "Util.hpp"

#include <iostream>
#include <iomanip>

#include <map>

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
    EventData& event_data = data_[event];

    // Accumulate event times in the log
    double elapsed = end_time - event_data.start_time;
    event_data.total_time += elapsed;
    event_data.hit += 1;

    if (print_event)
      std::cout << event_data << std::endl;

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
