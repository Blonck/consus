#pragma once

#include <vector>
#include <tuple>
#include <algorithm>
#include <chrono>
#include "reweight.hpp"
#include "../DiscreteAxis2D.hpp"
#include "../vector/helper.hpp"
#include "../addlogwise.hpp"

namespace consus
{

namespace WHAM2D
{

template <class T = double>
T reweight(const DiscreteAxis2D<T>& DOS, const DiscreteAxis2D<T>& MicroMean,
           const std::pair<T, T>& Parameter) {
  T min = std::numeric_limits<T>::max();
  for (size_t i = 0; i < MicroMean.size(); ++i) {
    if (std::isfinite(MicroMean[i])) {
      min = std::min(min, MicroMean[i]);
    }
  }
  T shift = min - 0.01;

  T lnZ = log_zero<T>();
  T tmp = log_zero<T>();
  for (int i = 0; i < DOS.get_num_bins_first(); ++i) {
    const auto E1 = DOS.get_value_first(i);
    for (int j = 0; j < DOS.get_num_bins_second(); ++j) {
      const auto E2 = DOS.get_value_second(j);
      const int index = DOS.get_index(i, j);
      if (std::isfinite(DOS[index])) {
        lnZ = addlogwise(
            lnZ, DOS[index] - Parameter.first * (E1 + Parameter.second * E2));
      }
    }
  }

  for (int i = 0; i < DOS.get_num_bins_first(); ++i) {
    const auto E1 = DOS.get_value_first(i);
    for (int j = 0; j < DOS.get_num_bins_second(); ++j) {
      const auto E2 = DOS.get_value_second(j);
      const int index = DOS.get_index(i, j);
      if (std::isfinite(DOS[index]) and std::isfinite(MicroMean[index])) {
        tmp = addlogwise(
            tmp, DOS[index] + std::log(MicroMean[index] - shift) -
                     Parameter.first * (E1 + Parameter.second * E2) - lnZ);
      }
    }
  }
  return std::exp(tmp) + shift;
}

template <class T = double>
vec1<T> reweight(const DiscreteAxis2D<T>& DOS, const DiscreteAxis2D<T>& MicroMean,
               const vec1<std::pair<T, T>>& Parameters) {
  T min = std::numeric_limits<T>::max();
  for (size_t i = 0; i < MicroMean.size(); ++i) {
    if (std::isfinite(MicroMean[i])) {
      min = std::min(min, MicroMean[i]);
    }
  }
  const T shift = min - 0.01;

  vec1<T> tmp(Parameters.size(), log_zero<T>());
  for (size_t k = 0; k < Parameters.size(); ++k) {
    auto lnZ = log_zero<T>();
    for (int i = 0; i < DOS.get_num_bins_first(); ++i) {
      const auto E1 = DOS.get_value_first(i);
      for (int j = 0; j < DOS.get_num_bins_second(); ++j) {
        const auto E2 = DOS.get_value_second(j);
        const int index = DOS.get_index(i, j);
        if (std::isfinite(DOS[index])) {
          lnZ = addlogwise(
              lnZ, DOS[index] -
                       Parameters[k].first * (E1 + Parameters[k].second * E2));
        }
      }
    }
    for (int i = 0; i < DOS.get_num_bins_first(); ++i) {
      const auto E1 = DOS.get_value_first(i);
      for (int j = 0; j < DOS.get_num_bins_second(); ++j) {
        const auto E2 = DOS.get_value_second(j);
        const int index = DOS.get_index(i, j);
        if (std::isfinite(DOS[index]) and std::isfinite(MicroMean[index])) {
          tmp[k] = addlogwise(
              tmp[k],
              DOS[index] + std::log(MicroMean[index] - shift) -
                  Parameters[k].first * (E1 + Parameters[k].second * E2) - lnZ);
        }
      }
    }
    tmp[k] = std::exp(tmp[k]) + shift;
  }
  return tmp;
}


  
} /* end of namespace WHAM2D */ 
  
} /* end of namespace consus */ 
