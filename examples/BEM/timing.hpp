#pragma once

#include <sys/time.h>

double get_time() {
  struct timeval tv;
  gettimeofday(&tv,NULL);
  return (double)(tv.tv_sec + 1e-6*tv.tv_usec);
}

