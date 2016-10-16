#pragma once

#include <vector>
#include <tuple>
#include <algorithm>
#include <chrono>
#include "reweight.hpp"
#include "../DiscreteAxis.hpp"
#include "../vector/helper.hpp"
#include "../addlogwise.hpp"

namespace consus
{

namespace WHAM
{

template <class TEnsemble>
double reweight(const DiscreteAxis& DOS, const DiscreteAxis& MicroMean,
                double Parameter) {
  double min = std::numeric_limits<double>::max();
  for (int i = 0; i < MicroMean.size(); ++i) {
    if (std::isfinite(MicroMean[i])) {
      min = std::min(min, MicroMean[i]);
    }
  }
  double shift = min - 0.01;

  double lnZ = log_zero<double>();
  double tmp = log_zero<double>();
  for (int i = 0; i < DOS.get_num_bins(); ++i) {
    const auto energy = DOS.get_value(i);
    if (std::isfinite(DOS[i])) {
      lnZ = addlogwise(lnZ, DOS[i] + TEnsemble::log_weight(Parameter, energy));
    }
  }

  for (int i = 0; i < DOS.get_num_bins(); ++i) {
    const auto energy = DOS.get_value(i);
    if (std::isfinite(DOS[i]) and std::isfinite(MicroMean[i])) {
      tmp = addlogwise(tmp, DOS[i] + std::log(MicroMean[i] - shift) +
                                TEnsemble::log_weight(Parameter, energy) - lnZ);
    }
  }
  return std::exp(tmp) + shift;
}

template <class TEnsemble>
vec1<double> reweight(const DiscreteAxis& DOS,
                      const DiscreteAxis& MicroMean,
                      const vec1<double>& Parameters) {
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
    for (int i = 0; i < DOS.get_num_bins(); ++i) {
      const auto energy = DOS.get_value(i);
      if (std::isfinite(DOS[i])) {
        lnZ = addlogwise(lnZ,
                         DOS[i] + TEnsemble::log_weight(Parameters[k], energy));
      }
    }
    for (int i = 0; i < DOS.get_num_bins(); ++i) {
      const auto energy = DOS.get_value(i);
      if (std::isfinite(DOS[i]) and std::isfinite(MicroMean[i])) {
        tmp[k] = addlogwise(
            tmp[k], DOS[i] + std::log(MicroMean[i] - shift) +
                        TEnsemble::log_weight(Parameters[k], energy) - lnZ);
      }
    }
    tmp[k] = std::exp(tmp[k]) + shift;
  }
  return tmp;
}

template <class TEnsemble>
double reweight_dT(const DiscreteAxis& DOS, const DiscreteAxis& MicroMean,
                   double Parameter) {
  double min = std::numeric_limits<double>::max();
  double minE = std::numeric_limits<double>::max();
  double minOE = std::numeric_limits<double>::max();
  for (int i = 0; i < DOS.get_num_bins(); ++i) {
    const double energy = DOS.get_value(i);
    if (std::isfinite(MicroMean[i])) {
      min = std::min(min, MicroMean[i]);
      minE = std::min(minE, energy);
      minOE = std::min(minOE, MicroMean[i] * energy);
    }
  }
  double shiftO = min - 0.1;
  double shiftE = minE - 0.1;
  double shiftOE = minOE - 0.1;

  double lnZ = log_zero<double>();
  for (int i = 0; i < DOS.get_num_bins(); ++i){
    const double energy = DOS.get_value(i);
      if (std::isfinite(DOS[i])){
        lnZ =
            addlogwise(lnZ, DOS[i] + TEnsemble::log_weight(Parameter, energy));
      }
  }

  double O = log_zero<double>();
  double OE =  log_zero<double>();
  double E =  log_zero<double>();
  for (int i = 0; i < DOS.get_num_bins(); ++i) {
    const double energy = DOS.get_value(i);
    if (std::isfinite(DOS[i]) and std::isfinite(MicroMean[i])) {
      double tmp = DOS[i] + TEnsemble::log_weight(Parameter, energy) - lnZ;
      O = addlogwise(O, tmp + std::log(MicroMean[i] - shiftO));
      OE = addlogwise(OE, tmp + std::log((energy)*MicroMean[i] - shiftOE));
      E = addlogwise(E, tmp + std::log((energy)-shiftE));
    }
  }
  O = std::exp(O) + shiftO;
  OE = std::exp(OE) + shiftOE;
  E = std::exp(E) + shiftE;
  return Parameter * Parameter * (OE - O * E);
}

template <class TEnsemble>
vec1<double> reweight_dT(const DiscreteAxis& DOS,
                         const DiscreteAxis& MicroMean,
                         const vec1<double>& Parameters) {
  double min = std::numeric_limits<double>::max();
  double minE = std::numeric_limits<double>::max();
  double minOE = std::numeric_limits<double>::max();
  // TODO: for min no reduction is necessary, could be calculated at k == 0
  #pragma omp parallel for reduction(min : min, minE, minOE)
  for (int i = 0; i < DOS.get_num_bins(); ++i) {
    const double energy = DOS.get_value(i);
    if (std::isfinite(MicroMean[i])) {
      min = std::min(min, MicroMean[i]);
      minE = std::min(minE, energy);
      minOE = std::min(minOE, MicroMean[i] * energy);
    }
  }
  double shiftO = min - 0.1;
  double shiftE = minE - 0.1;
  double shiftOE = minOE - 0.1;

  vec1<double> tmp(Parameters.size(), log_zero<double>());
#pragma omp parallel for
  for (size_t k = 0; k < Parameters.size(); ++k) {
    double lnZ = log_zero<double>();
    for (int i = 0; i < DOS.get_num_bins(); ++i) {
      const double energy = DOS.get_value(i);
      if (std::isfinite(DOS[i])) {
        lnZ = addlogwise(lnZ,
                         DOS[i] + TEnsemble::log_weight(Parameters[k], energy));
      }
    }

    double O = log_zero<double>();
    double OE = log_zero<double>();
    double E = log_zero<double>();
    for (int i = 0; i < DOS.get_num_bins(); ++i) {
      const double energy = DOS.get_value(i);
      if (std::isfinite(DOS[i]) and std::isfinite(MicroMean[i])) {
        double tmp =
            DOS[i] + TEnsemble::log_weight(Parameters[k], energy) - lnZ;
        O = addlogwise(O, tmp + std::log(MicroMean[i] - shiftO));
        OE = addlogwise(OE, tmp + std::log((energy)*MicroMean[i] - shiftOE));
        E = addlogwise(E, tmp + std::log((energy)-shiftE));
      }
    }
    O = std::exp(O) + shiftO;
    OE = std::exp(OE) + shiftOE;
    E = std::exp(E) + shiftE;
    tmp[k] = Parameters[k] * Parameters[k] * (OE - O * E);
  }
  return tmp;
}

template <class TEnsemble>
double reweight_dT2(const DiscreteAxis& DOS, const DiscreteAxis& MicroMean,
                    double Parameter) {
  double min = std::numeric_limits<double>::max();
  double minE = std::numeric_limits<double>::max();
  double minEE = std::numeric_limits<double>::max();
  double minOE = std::numeric_limits<double>::max();
  double minOEE = std::numeric_limits<double>::max();
  for (int i = 0; i < DOS.get_num_bins(); ++i) {
    const double energy = DOS.get_value(i);
    if (std::isfinite(MicroMean[i])) {
      min = std::min(min, MicroMean[i]);
      minE = std::min(minE, energy);
      minEE = std::min(minEE, energy * energy);
      minOE = std::min(minOE, MicroMean[i] * energy);
      minOEE = std::min(minOEE, MicroMean[i] * energy * energy);
    }
  }
  double shiftO = min - 0.1;
  double shiftE = minE - 0.1;
  double shiftEE = minEE - 0.1;
  double shiftOE = minOE - 0.1;
  double shiftOEE = minOEE - 0.1;

  double lnZ = log_zero<double>();
  for (int i = 0; i < DOS.get_num_bins(); ++i) {
    const double energy = DOS.get_value(i);
    if (std::isfinite(DOS[i])) {
      lnZ = addlogwise(lnZ, DOS[i] + TEnsemble::log_weight(Parameter, energy));
    }
  }

  double O = log_zero<double>();
  double OE = log_zero<double>();
  double OEE = log_zero<double>();
  double E = log_zero<double>();
  double EE = log_zero<double>();
  for (int i = 0; i < DOS.get_num_bins(); ++i){
    const double energy = DOS.get_value(i);
    if (std::isfinite(DOS[i]) and std::isfinite(MicroMean[i])) {
      double tmp = DOS[i] + TEnsemble::log_weight(Parameter, energy) - lnZ;
      O = addlogwise(O, tmp + std::log(MicroMean[i] - shiftO));
      OE = addlogwise(OE, tmp + std::log(energy * MicroMean[i] - shiftOE));
      OEE = addlogwise(
          OEE, tmp + std::log(energy * energy * MicroMean[i] - shiftOEE));
      E = addlogwise(E, tmp + std::log(energy - shiftE));
      EE = addlogwise(EE, tmp + std::log(energy * energy - shiftEE));
    }
  }
  O = std::exp(O) + shiftO;
  OE = std::exp(OE) + shiftOE;
  OEE = std::exp(OEE) + shiftOEE;
  E = std::exp(E) + shiftE;
  EE = std::exp(EE) + shiftEE;
  return 2 * Parameter * Parameter * Parameter * (O * E - OE) +
         Parameter * Parameter * Parameter * Parameter *
             (OEE - O * EE - 2 * OE * E + 2 * O * E * E);
}

template <class TEnsemble>
vec1<double> reweight_dT2(const DiscreteAxis& DOS,
                          const DiscreteAxis& MicroMean,
                          const vec1<double>& Parameters) {
  double min = std::numeric_limits<double>::max();
  double minE = std::numeric_limits<double>::max();
  double minEE = std::numeric_limits<double>::max();
  double minOE = std::numeric_limits<double>::max();
  double minOEE = std::numeric_limits<double>::max();

  #pragma omp parallel for reduction(min : min, minE, minEE, minOE, minOEE)
  for (int i = 0; i < DOS.get_num_bins(); ++i) {
    const double energy = DOS.get_value(i);
    if (std::isfinite(MicroMean[i])) {
      min = std::min(min, MicroMean[i]);
      for (size_t k = 0; k < Parameters.size(); ++k) {
        minE = std::min(minE, energy);
        minEE = std::min(minEE, energy * energy);
        minOE = std::min(minOE, MicroMean[i] * energy);
        minOEE = std::min(minOEE, MicroMean[i] * energy * energy);
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
    for (int i = 0; i < DOS.get_num_bins(); ++i) {
      const double energy = DOS.get_value(i);
        if (std::isfinite(DOS[i])) {
          lnZ = addlogwise(
              lnZ, DOS[i] + TEnsemble::log_weight(Parameters[k], energy));
        }
    }

    double O = log_zero<double>();
    double OE = log_zero<double>();
    double OEE = log_zero<double>();
    double E = log_zero<double>();
    double EE = log_zero<double>();
    for (int i = 0; i < DOS.get_num_bins(); ++i) {
      const double energy = DOS.get_value(i);
      if (std::isfinite(DOS[i]) and std::isfinite(MicroMean[i])) {
        double tmp =
            DOS[i] + TEnsemble::log_weight(Parameters[k], energy) - lnZ;
        O = addlogwise(O, tmp + std::log(MicroMean[i] - shiftO));
        OE = addlogwise(OE, tmp + std::log(energy * MicroMean[i] - shiftOE));
        OEE = addlogwise(
            OEE, tmp + std::log(energy * energy * MicroMean[i] - shiftOEE));
        E = addlogwise(E, tmp + std::log(energy - shiftE));
        EE = addlogwise(EE, tmp + std::log(energy * energy - shiftEE));
      }
    }
    O = std::exp(O) + shiftO;
    OE = std::exp(OE) + shiftOE;
    OEE = std::exp(OEE) + shiftOEE;
    E = std::exp(E) + shiftE;
    EE = std::exp(EE) + shiftEE;
    tmp[k] = 2 * Parameters[k] * Parameters[k] * Parameters[k] * (O * E - OE) +
             Parameters[k] * Parameters[k] * Parameters[k] * Parameters[k] *
                 (OEE - O * EE - 2 * OE * E + 2 * O * E * E);
  }
  return tmp;
}

template <class TEnsemble>
vec1<double> reweight_prob(const DiscreteAxis& DOS,
                           const DiscreteAxis& MicroMean,
                           const std::pair<double, double>& Parameter,
                           vec1<double>& bounds) {
  std::sort(bounds.begin(), bounds.end());
  vec1<double> result(bounds.size()+1, log_zero<double>());
  double lnZ = log_zero<double>();
  for (int i = 0; i < DOS.get_num_bins(); ++i) {
    const double energy = DOS.get_value(i);
    if (std::isfinite(DOS[i])) {
      auto logweight = TEnsemble::log_weight(Parameter.first, energy);
      lnZ = addlogwise(lnZ, DOS[i] + logweight);
      if (std::isfinite(MicroMean[i])) {
        const double value = MicroMean[i];
        auto upper = std::upper_bound(bounds.begin(), bounds.end(), value);
        size_t num = std::distance(bounds.begin(), upper);
        if (upper != bounds.end()) {
          result[num] = addlogwise(result[num], DOS[i] + logweight);
        } else {
          result.back() = addlogwise(result.back(), DOS[i] + logweight);
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
DiscreteAxis reweight_hist(const DiscreteAxis& DOS, double Parameter) {
  DiscreteAxis result(DOS.get_min(), DOS.get_max(), DOS.get_step(), 0.0);
  double lnZ = log_zero<double>();
  for (int i = 0; i < DOS.get_num_bins(); ++i) {
    const double energy = DOS.get_value(i);
    auto logweight = TEnsemble::log_weight(Parameter, energy);
    if (std::isfinite(DOS[i])) {
      result[i] = DOS[i] + logweight;
      lnZ = addlogwise(lnZ, result[i]);
    }
  }
  for (int i = 0; i < result.size(); ++i) {
    result[i] -= lnZ;
    result[i] = std::exp(result[i]);
  }
  return result;
}

} /* end of namespace WHAM2D */ 
  
} /* end of namespace consus */ 
