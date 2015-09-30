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

template <class TEnsemble>
double reweight(const DiscreteAxis2D& DOS, const DiscreteAxis2D& MicroMean,
                const std::pair<double, double>& Parameter) {
  double min = std::numeric_limits<double>::max();
  for (int i = 0; i < MicroMean.size(); ++i) {
    if (std::isfinite(MicroMean[i])) {
      min = std::min(min, MicroMean[i]);
    }
  }
  double shift = min - 0.01;

  double lnZ = log_zero<double>();
  double tmp = log_zero<double>();
  for (int i = 0; i < DOS.get_num_bins_first(); ++i) {
    const auto E1 = DOS.get_value_first(i);
    for (int j = 0; j < DOS.get_num_bins_second(); ++j) {
      const auto E2 = DOS.get_value_second(j);
      const int index = DOS.get_index(i, j);
      if (std::isfinite(DOS[index])) {
        lnZ = addlogwise(
            lnZ, DOS[index] + TEnsemble::log_weight(Parameter.first, E1,
                                                    Parameter.second, E2));
      }
    }
  }

  for (int i = 0; i < DOS.get_num_bins_first(); ++i) {
    const auto E1 = DOS.get_value_first(i);
    for (int j = 0; j < DOS.get_num_bins_second(); ++j) {
      const auto E2 = DOS.get_value_second(j);
      const int index = DOS.get_index(i, j);
      if (std::isfinite(DOS[index]) and std::isfinite(MicroMean[index])) {
        tmp = addlogwise(tmp, DOS[index] + std::log(MicroMean[index] - shift) +
                                  TEnsemble::log_weight(Parameter.first, E1,
                                                        Parameter.second, E2) -
                                  lnZ);
      }
    }
  }
  return std::exp(tmp) + shift;
}

template <class TEnsemble>
vec1<double> reweight(const DiscreteAxis2D& DOS,
                      const DiscreteAxis2D& MicroMean,
                      const vec1<std::pair<double, double>>& Parameters) {
  double min = std::numeric_limits<double>::max();
  for (int i = 0; i < MicroMean.size(); ++i) {
    if (std::isfinite(MicroMean[i])) {
      min = std::min(min, MicroMean[i]);
    }
  }
  const double shift = min - 0.01;

  vec1<double> tmp(Parameters.size(), log_zero<double>());
  for (size_t k = 0; k < Parameters.size(); ++k) {
    auto lnZ = log_zero<double>();
    for (int i = 0; i < DOS.get_num_bins_first(); ++i) {
      const auto E1 = DOS.get_value_first(i);
      for (int j = 0; j < DOS.get_num_bins_second(); ++j) {
        const auto E2 = DOS.get_value_second(j);
        const int index = DOS.get_index(i, j);
        if (std::isfinite(DOS[index])) {
          lnZ = addlogwise(lnZ, DOS[index] + TEnsemble::log_weight(
                                                 Parameters[k].first, E1,
                                                 Parameters[k].second, E2));
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
              tmp[k], DOS[index] + std::log(MicroMean[index] - shift) +
                          TEnsemble::log_weight(Parameters[k].first, E1,
                                                Parameters[k].second, E2) -
                          lnZ);
        }
      }
    }
    tmp[k] = std::exp(tmp[k]) + shift;
  }
  return tmp;
}

} /* end of namespace WHAM2D */ 
  
} /* end of namespace consus */ 
