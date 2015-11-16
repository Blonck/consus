#pragma once

#include <vector>
#include <tuple>
#include <algorithm>
#include <chrono>
#include "reweight.hpp"
#include "../DiscreteAxis.hpp"
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
  #pragma omp parallel for
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

template <class TEnsemble>
double reweight_dT(const DiscreteAxis2D& DOS, const DiscreteAxis2D& MicroMean,
                   const std::pair<double, double>& Parameter) {
  double min = std::numeric_limits<double>::max();
  double minE = std::numeric_limits<double>::max();
  double minOE = std::numeric_limits<double>::max();
  for (int i = 0; i < DOS.get_num_bins_first(); ++i) {
    const double E1 = DOS.get_value_first(i);
    for (int j = 0; j < DOS.get_num_bins_second(); ++j) {
      const int index = DOS.get_index(i, j);
      if (std::isfinite(MicroMean[index])) {
        const double E2 = DOS.get_value_second(j);
        const double E = TEnsemble::ham(E1, Parameter.second, E2);
        min = std::min(min, MicroMean[index]);
        minE = std::min(minE, E);
        minOE = std::min(minOE, MicroMean[index] * E);
      }
    }
  }
  double shiftO = min - 0.1;
  double shiftE = minE - 0.1;
  double shiftOE = minOE - 0.1;

  double lnZ = log_zero<double>();
  for (int i = 0; i < DOS.get_num_bins_first(); ++i){
    const double E1 = DOS.get_value_first(i);
    for (int j = 0; j < DOS.get_num_bins_second(); ++j){
      const double E2 = DOS.get_value_second(j);
      const int index = DOS.get_index(i, j);
      if (std::isfinite(DOS[index])){
        lnZ = addlogwise(
            lnZ, DOS[index] + TEnsemble::log_weight(Parameter.first, E1,
                                                    Parameter.second, E2));
      }
    }
  }

  double O = log_zero<double>();
  double OE =  log_zero<double>();
  double E =  log_zero<double>();
  for (int i = 0; i < DOS.get_num_bins_first(); ++i){
    const double e1 = DOS.get_value_first(i);
    for (int j = 0; j < DOS.get_num_bins_second(); ++j){
      const double e2 = DOS.get_value_second(j);
      const int index = DOS.get_index(i, j);
      if (std::isfinite(DOS[index]) and std::isfinite(MicroMean[index])) {
        double tmp = DOS[index] + TEnsemble::log_weight(Parameter.first, e1,
                                                        Parameter.second, e2) -
                     lnZ;
        double e = TEnsemble::ham(e1, Parameter.second, e2);
        O = addlogwise(O, tmp + std::log(MicroMean[index] - shiftO));
        OE = addlogwise(OE, tmp + std::log((e)*MicroMean[index] - shiftOE));
        E = addlogwise(E, tmp + std::log((e) - shiftE));
      }
    }
  }
  O = std::exp(O) + shiftO;
  OE = std::exp(OE) + shiftOE;
  E = std::exp(E) + shiftE;
  return Parameter.first * Parameter.first * ( OE - O * E );
}

template <class TEnsemble>
vec1<double> reweight_dT(const DiscreteAxis2D& DOS,
                         const DiscreteAxis2D& MicroMean,
                         const vec1<std::pair<double, double>>& Parameters) {
  auto red_params = Parameters;
  typedef std::pair<double, double> TP;
  std::sort(red_params.begin(), red_params.end(),
            [](const TP& a, const TP& b) { return a.second < b.second; });
  red_params.erase(std::unique(red_params.begin(), red_params.end()),
                   red_params.end());
  double min = std::numeric_limits<double>::max();
  double minE = std::numeric_limits<double>::max();
  double minOE = std::numeric_limits<double>::max();
  // TODO: for min no reduction is necessary, could be calculated at k == 0
  #pragma omp parallel for reduction(min : min, minE, minOE)
  for (int i = 0; i < DOS.get_num_bins_first(); ++i) {
    const double E1 = DOS.get_value_first(i);
    for (int j = 0; j < DOS.get_num_bins_second(); ++j) {
      const double E2 = DOS.get_value_second(j);
      const int index = DOS.get_index(i, j);
      if (std::isfinite(MicroMean[index])) {
        min = std::min(min, MicroMean[index]);
        for (size_t k = 0; k < red_params.size(); ++k) {
          const double E = TEnsemble::ham(E1, red_params[k].second, E2);
          minE = std::min(minE, E);
          minOE = std::min(minOE, MicroMean[index] * E);
        }
      }
    }
  }
  double shiftO = min - 0.1;
  double shiftE = minE - 0.1;
  double shiftOE = minOE - 0.1;

  vec1<double> tmp(Parameters.size(), log_zero<double>());
#pragma omp parallel for
  for (size_t k = 0; k < Parameters.size(); ++k) {
    double lnZ = log_zero<double>();
    for (int i = 0; i < DOS.get_num_bins_first(); ++i) {
      const double E1 = DOS.get_value_first(i);
      for (int j = 0; j < DOS.get_num_bins_second(); ++j) {
        const double E2 = DOS.get_value_second(j);
        const int index = DOS.get_index(i, j);
        if (std::isfinite(DOS[index])) {
          lnZ = addlogwise(lnZ, DOS[index] + TEnsemble::log_weight(
                                                 Parameters[k].first, E1,
                                                 Parameters[k].second, E2));
        }
      }
    }

    double O = log_zero<double>();
    double OE = log_zero<double>();
    double E = log_zero<double>();
    for (int i = 0; i < DOS.get_num_bins_first(); ++i) {
      const double e1 = DOS.get_value_first(i);
      for (int j = 0; j < DOS.get_num_bins_second(); ++j) {
        const double e2 = DOS.get_value_second(j);
        const int index = DOS.get_index(i, j);
        if (std::isfinite(DOS[index]) and std::isfinite(MicroMean[index])) {
          double tmp = DOS[index] +
                       TEnsemble::log_weight(Parameters[k].first, e1,
                                             Parameters[k].second, e2) -
                       lnZ;
          double e = TEnsemble::ham(e1, Parameters[k].second, e2);
          O = addlogwise(O, tmp + std::log(MicroMean[index] - shiftO));
          OE = addlogwise(OE, tmp + std::log((e)*MicroMean[index] - shiftOE));
          E = addlogwise(E, tmp + std::log((e)-shiftE));
        }
      }
    }
    O = std::exp(O) + shiftO;
    OE = std::exp(OE) + shiftOE;
    E = std::exp(E) + shiftE;
    tmp[k] = Parameters[k].first * Parameters[k].first * (OE - O * E);
  }
  return tmp;
}

template <class TEnsemble>
double reweight_dT2(const DiscreteAxis2D& DOS, const DiscreteAxis2D& MicroMean,
                    const std::pair<double, double>& Parameter) {
  double min = std::numeric_limits<double>::max();
  double minE = std::numeric_limits<double>::max();
  double minEE = std::numeric_limits<double>::max();
  double minOE = std::numeric_limits<double>::max();
  double minOEE = std::numeric_limits<double>::max();
  for (int i = 0; i < DOS.get_num_bins_first(); ++i) {
    const double E1 = DOS.get_value_first(i);
    for (int j = 0; j < DOS.get_num_bins_second(); ++j) {
      const double E2 = DOS.get_value_second(j);
      const int index = DOS.get_index(i, j);
      if (std::isfinite(MicroMean[index])) {
        double tmpE = TEnsemble::ham(E1, Parameter.second, E2);
        min = std::min(min, MicroMean[index]);
        minE = std::min(minE, tmpE);
        minEE = std::min(minEE, tmpE * tmpE);
        minOE = std::min(minOE, MicroMean[index] * tmpE);
        minOEE = std::min(minOEE, MicroMean[index] * tmpE * tmpE);
      }
    }
  }
  double shiftO = min - 0.1;
  double shiftE = minE - 0.1;
  double shiftEE = minEE - 0.1;
  double shiftOE = minOE - 0.1;
  double shiftOEE = minOEE - 0.1;

  double lnZ = log_zero<double>();
  for (int i = 0; i < DOS.get_num_bins_first(); ++i){
    const double E1 = DOS.get_value_first(i);
    for (int j = 0; j < DOS.get_num_bins_second(); ++j){
      const double E2 = DOS.get_value_second(j);
      const int index = DOS.get_index(i, j);
      if (std::isfinite(DOS[index])){
        lnZ = addlogwise(
            lnZ, DOS[index] +
                     TEnsemble::log_weight(Parameter.first, E1, Parameter.second, E2));
      }
    }
  }

  double O = log_zero<double>();
  double OE = log_zero<double>();
  double OEE = log_zero<double>();
  double E = log_zero<double>();
  double EE = log_zero<double>();
  for (int i = 0; i < DOS.get_num_bins_first(); ++i){
    const double e1 = DOS.get_value_first(i);
    for (int j = 0; j < DOS.get_num_bins_second(); ++j){
      const double e2 = DOS.get_value_second(j);
      const int index = DOS.get_index(i, j);
      if (std::isfinite(DOS[index]) and std::isfinite(MicroMean[index])) {
        double tmpE = TEnsemble::ham(e1, Parameter.second, e2);
        double tmp = DOS[index] - Parameter.first * (tmpE) - lnZ;
        O = addlogwise(O, tmp + std::log(MicroMean[index] - shiftO));
        OE = addlogwise(OE, tmp + std::log(tmpE * MicroMean[index] - shiftOE));
        OEE = addlogwise(
            OEE, tmp + std::log(tmpE * tmpE * MicroMean[index] - shiftOEE));
        E = addlogwise(E, tmp + std::log(tmpE - shiftE));
        EE = addlogwise(EE, tmp + std::log(tmpE * tmpE - shiftEE));
      }
    }
  }
  O = std::exp(O) + shiftO;
  OE = std::exp(OE) + shiftOE;
  OEE = std::exp(OEE) + shiftOEE;
  E = std::exp(E) + shiftE;
  EE = std::exp(EE) + shiftEE;
  return 2 * Parameter.first * Parameter.first * Parameter.first *
             (O * E - OE) +
         Parameter.first * Parameter.first * Parameter.first * Parameter.first *
             (OEE - O * EE - 2 * OE * E + 2 * O * E * E);
}

template <class TEnsemble>
vec1<double> reweight_dT2(const DiscreteAxis2D& DOS,
                          const DiscreteAxis2D& MicroMean,
                          const vec1<std::pair<double, double>>& Parameters) {
  auto red_params = Parameters;
  typedef std::pair<double, double> TP;
  std::sort(red_params.begin(), red_params.end(),
            [](const TP& a, const TP& b) { return a.second < b.second; });
  red_params.erase(std::unique(red_params.begin(), red_params.end()),
                   red_params.end());
  double min = std::numeric_limits<double>::max();
  double minE = std::numeric_limits<double>::max();
  double minEE = std::numeric_limits<double>::max();
  double minOE = std::numeric_limits<double>::max();
  double minOEE = std::numeric_limits<double>::max();
  // TODO: for min no reduction is necessary, could be calculated at k == 0
  #pragma omp parallel for reduction(min : min, minE, minEE, minOE, minOEE)
  for (int i = 0; i < DOS.get_num_bins_first(); ++i) {
    const double E1 = DOS.get_value_first(i);
    for (int j = 0; j < DOS.get_num_bins_second(); ++j) {
      const double E2 = DOS.get_value_second(j);
      const int index = DOS.get_index(i, j);
      if (std::isfinite(MicroMean[index])) {
        min = std::min(min, MicroMean[index]);
        for (size_t k = 0; k < red_params.size(); ++k) {
          const double E = TEnsemble::ham(E1, red_params[k].second, E2);
          minE = std::min(minE, E);
          minEE = std::min(minEE, E * E);
          minOE = std::min(minOE, MicroMean[index] * E);
          minOEE = std::min(minOEE, MicroMean[index] * E * E);
        }
      }
    }
  }
  double shiftO = min - 0.1;
  double shiftE = minE - 0.1;
  double shiftEE = minEE - 0.1;
  double shiftOE = minOE - 0.1;
  double shiftOEE = minOEE - 0.1;

  vec1<double> tmp(Parameters.size(), log_zero<double>());
  #pragma omp parallel for
  for (size_t k = 0; k < Parameters.size(); ++k) {
    double lnZ = log_zero<double>();
    for (int i = 0; i < DOS.get_num_bins_first(); ++i) {
      const double E1 = DOS.get_value_first(i);
      for (int j = 0; j < DOS.get_num_bins_second(); ++j) {
        const double E2 = DOS.get_value_second(j);
        const int index = DOS.get_index(i, j);
        if (std::isfinite(DOS[index])) {
          lnZ = addlogwise(lnZ, DOS[index] + TEnsemble::log_weight(
                                                 Parameters[k].first, E1,
                                                 Parameters[k].second, E2));
        }
      }
    }

    double O = log_zero<double>();
    double OE = log_zero<double>();
    double OEE = log_zero<double>();
    double E = log_zero<double>();
    double EE = log_zero<double>();
    for (int i = 0; i < DOS.get_num_bins_first(); ++i) {
      const double e1 = DOS.get_value_first(i);
      for (int j = 0; j < DOS.get_num_bins_second(); ++j) {
        const double e2 = DOS.get_value_second(j);
        const int index = DOS.get_index(i, j);
        if (std::isfinite(DOS[index]) and std::isfinite(MicroMean[index])) {
          double tmpE = TEnsemble::ham(e1, Parameters[k].second, e2);
          double tmp = DOS[index] - Parameters[k].first * (tmpE)-lnZ;
          O = addlogwise(O, tmp + std::log(MicroMean[index] - shiftO));
          OE =
              addlogwise(OE, tmp + std::log(tmpE * MicroMean[index] - shiftOE));
          OEE = addlogwise(
              OEE, tmp + std::log(tmpE * tmpE * MicroMean[index] - shiftOEE));
          E = addlogwise(E, tmp + std::log(tmpE - shiftE));
          EE = addlogwise(EE, tmp + std::log(tmpE * tmpE - shiftEE));
        }
      }
    }
    O = std::exp(O) + shiftO;
    OE = std::exp(OE) + shiftOE;
    OEE = std::exp(OEE) + shiftOEE;
    E = std::exp(E) + shiftE;
    EE = std::exp(EE) + shiftEE;
    tmp[k] = 2 * Parameters[k].first * Parameters[k].first *
                 Parameters[k].first * (O * E - OE) +
             Parameters[k].first * Parameters[k].first * Parameters[k].first *
                 Parameters[k].first *
                 (OEE - O * EE - 2 * OE * E + 2 * O * E * E);
  }
  return tmp;
}

template <class TEnsemble>
vec1<double> reweight_prob(const DiscreteAxis2D& DOS,
                           const DiscreteAxis2D& MicroMean,
                           const std::pair<double, double>& Parameter,
                           vec1<double>& bounds) {
  std::sort(bounds.begin(), bounds.end());
  vec1<double> result(bounds.size()+1, log_zero<double>());
  double lnZ = log_zero<double>();
  for (int i = 0; i < DOS.get_num_bins_first(); ++i){
    const double E1 = DOS.get_value_first(i);
    for (int j = 0; j < DOS.get_num_bins_second(); ++j){
      const double E2 = DOS.get_value_second(j);
      const int index = DOS.get_index(i, j);
      if (std::isfinite(DOS[index])){
        auto logweight =
            TEnsemble::log_weight(Parameter.first, E1, Parameter.second, E2);
        lnZ = addlogwise(lnZ, DOS[index] + logweight);
        if (std::isfinite(MicroMean[index])) {
          const double value = MicroMean[index];
          auto upper = std::upper_bound(bounds.begin(), bounds.end(), value);
          size_t num = std::distance(bounds.begin(), upper);
          if (upper != bounds.end()) {
            result[num] = addlogwise(result[num], DOS[index] + logweight);
          } else {
            result.back() = addlogwise(result.back(), DOS[index] + logweight);
          }
        }
      }
    }
  }
  for (size_t i = 0; i < result.size(); ++i){
    result[i] -= lnZ;
    result[i] = std::exp(result[i]);
  }
  return result;
}

template <class TEnsemble>
DiscreteAxis2D reweight_hist2d(const DiscreteAxis2D& DOS,
                               const std::pair<double, double>& Parameter) {
  DiscreteAxis2D result(DOS, 0.0);
  double lnZ = log_zero<double>();
  for (int i = 0; i < DOS.get_num_bins_first(); ++i){
    const double E1 = DOS.get_value_first(i);
    for (int j = 0; j < DOS.get_num_bins_second(); ++j){
      const double E2 = DOS.get_value_second(j);
      const int index = DOS.get_index(i, j);
      auto logweight =
          TEnsemble::log_weight(Parameter.first, E1, Parameter.second, E2);
      if (std::isfinite(DOS[index])){
        result[index] = DOS[index] + logweight;
        lnZ = addlogwise(lnZ, result[index]);
      }
    }
  }
  for (int i = 0; i < result.size(); ++i){
    result[i] -= lnZ;
    result[i] = std::exp(result[i]);
  }
  return result;
}

template <class TEnsemble>
DiscreteAxis reweight_hist(const DiscreteAxis2D& DOS, const double step,
                           const std::pair<double, double>& Parameter) {
  double minE = std::numeric_limits<double>::max();
  double maxE = std::numeric_limits<double>::lowest();
  for (int i = 0; i < DOS.get_num_bins_first(); ++i){
    const double E1 = DOS.get_value_first(i);
    for (int j = 0; j < DOS.get_num_bins_second(); ++j){
      const double E2 = DOS.get_value_second(j);
      const double E = TEnsemble::ham(E1, Parameter.second, E2);
      minE = std::min(minE, E);
      maxE = std::max(maxE, E);
    }
  }

  DiscreteAxis result(minE, maxE, step, log_zero<double>());
  double lnZ = log_zero<double>();
  for (int i = 0; i < DOS.get_num_bins_first(); ++i){
    const double E1 = DOS.get_value_first(i);
    for (int j = 0; j < DOS.get_num_bins_second(); ++j){
      const double E2 = DOS.get_value_second(j);
      auto logweight =
          TEnsemble::log_weight(Parameter.first, E1, Parameter.second, E2);
      const double E = TEnsemble::ham(E1, Parameter.second, E2);
      const int index = DOS.get_index(i, j);
      const int rindex = result.get_bin(E);
      if (std::isfinite(DOS[index])){
        result[rindex] = addlogwise(result[rindex], DOS[index] + logweight);
        lnZ = addlogwise(lnZ, DOS[index] + logweight);
      }
    }
  }
  for (int i = 0; i < result.size(); ++i){
      result[i] -= lnZ;
      result[i] = std::exp(result[i]);
  }
  return result;
}

} /* end of namespace WHAM2D */ 
  
} /* end of namespace consus */ 
