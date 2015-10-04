#pragma once

#include <cmath>
#include <limits>

template <class T = double>
constexpr T log_zero() {
  return std::numeric_limits<T>::lowest();
}

// TODO: maybe use __builtin_expect since a > b in almost all cases
// TODO: can be optimized if std::exp(x) = 0, return value don't need log1p(exp())
template <class T>
inline double addlogwise(T a, T b) noexcept {
  return (a > b) ? (a + std::log1p(std::exp(b - a)))
                 : (b + std::log1p(std::exp(a - b)));
}

//template <class T>
//inline T addlogwise(T a, T b) noexcept {
//  if (a < b){
//    std::swap(a, b);
//  }
//  if (b-a <= -1000.0){
//    return a;
//  }else{
//    return a + std::log1p(std::exp(b-a));
//  }
//}
