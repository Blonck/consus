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
  double min = std::numeric_limits<double>::max();
  double minE = std::numeric_limits<double>::max();
  double minOE = std::numeric_limits<double>::max();
  // TODO: for min no reduction is necessary, could be calculated at k == 0
  #pragma omp parallel for reduction(min : min, minE, minOE)
  for (size_t k = 0; k < Parameters.size(); ++k) {
    for (int i = 0; i < DOS.get_num_bins_first(); ++i) {
      const double E1 = DOS.get_value_first(i);
      for (int j = 0; j < DOS.get_num_bins_second(); ++j) {
        const int index = DOS.get_index(i, j);
        if (std::isfinite(MicroMean[index])) {
          const double E2 = DOS.get_value_second(j);
          const double E = TEnsemble::ham(E1, Parameters[k].second, E2);
          min = std::min(min, MicroMean[index]);
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

} /* end of namespace WHAM2D */ 
  
} /* end of namespace consus */ 
