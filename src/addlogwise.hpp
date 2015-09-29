#pragma once

#include <cmath>
#include <limits>

template <class T = double>
struct log_zero{
  constexpr T operator()() const {return std::numeric_limits<T>::lowest();}
};

template <class T>
inline double addlogwise(const T a, const T b) {
  T log_max;
  T log_abs;
  if (a > b) {
    log_max = a;
    log_abs = b - a;
  } else {
    log_max = b;
    log_abs = a - b;
  }
  return log_max + std::log1p(std::exp(log_abs));
}
