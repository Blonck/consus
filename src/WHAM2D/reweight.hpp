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

//vec1d reweight_dT(const DiscreteAxis2D& DOS, const DiscreteAxis2D& MicroMean,
//                  const vec1d& Betas, const double Parameter) {
//  double min = std::numeric_limits<double>::max();
//  double minE = std::numeric_limits<double>::max();
//  double minOE = std::numeric_limits<double>::max();
//  for (int i = 0; i < DOS.get_NumBins_first(); ++i) {
//    const double E1 = DOS.get_Value_first(i);
//    for (int j = 0; j < DOS.get_NumBins_second(); ++j) {
//      const int index = DOS.get_index(i, j);
//      if (std::isfinite(MicroMean[index])) {
//        const double E2 = DOS.get_Value_second(j);
//        min = std::min(min, MicroMean[index]);
//        minE = std::min(minE, E1 + Parameter * E2);
//        minOE = std::min(minOE, MicroMean[index] * (E1 + Parameter *E2));
//      }
//    }
//  }
//  double shiftO = min - 0.1;
//  double shiftE = minE - 0.1;
//  double shiftOE = minOE - 0.1;
//
//  vec1d result(Betas.size());
//
//  for (size_t k = 0; k < Betas.size(); ++k){
//    double lnZ = LOGZERO;
//    double O = LOGZERO;
//    double OE = LOGZERO;
//    double E = LOGZERO;
//    for (int i = 0; i < DOS.get_NumBins_first(); ++i){
//      const double E1 = DOS.get_Value_first(i);
//      for (int j = 0; j < DOS.get_NumBins_second(); ++j){
//        const double E2 = DOS.get_Value_second(j);
//        const int index = DOS.get_index(i, j);
//        if (std::isfinite(DOS[index])){
//          lnZ = addlogwise(lnZ, DOS[index] - Betas[k] * (E1 + Parameter * E2));
//        }
//      }
//    }
//
//    for (int i = 0; i < DOS.get_NumBins_first(); ++i){
//      const double e1 = DOS.get_Value_first(i);
//      for (int j = 0; j < DOS.get_NumBins_second(); ++j){
//        const double e2 = DOS.get_Value_second(j);
//        const int index = DOS.get_index(i, j);
//        if (std::isfinite(DOS[index]) and std::isfinite(MicroMean[index])) {
//          double tmp = DOS[index] - Betas[k] * (e1 + Parameter * e2) - lnZ;
//          double tmpE = e1 + Parameter * e2;
//          O = addlogwise(O, tmp + std::log(MicroMean[index] - shiftO));
//          OE =
//              addlogwise(OE, tmp + std::log(tmpE * MicroMean[index] - shiftOE));
//          E = addlogwise(E, tmp + std::log(tmpE - shiftE));
//        }
//      }
//    }
//    O = std::exp(O) + shiftO;
//    OE = std::exp(OE) + shiftOE;
//    E = std::exp(E) + shiftE;
//    result[k] = Betas[k] * Betas[k] * ( OE - O * E );
//  }
//  return result;
//}

} /* end of namespace WHAM2D */ 
  
} /* end of namespace consus */ 
