#pragma once

#include <vector>
#include <tuple>
#include <algorithm>
#include <chrono>
#include "../DiscreteAxis2D.hpp"
#include "../vector/helper.hpp"
#include "../addlogwise.hpp"

namespace consus
{
  
namespace WHAM2D
{

namespace detail {

// checks if two rectangles (in the filled entries in histogram) overlaps
inline bool overlaps(int ffi_first_1, int ffi_first_2, int lfi_first_1,
                     int lfi_first_2, int ffi_second_1, int ffi_second_2,
                     int lfi_second_1, int lfi_second_2) {
  if ((ffi_first_1 >= ffi_first_2 and ffi_first_1 <= lfi_first_2) or
      (lfi_first_1 >= ffi_first_2 and lfi_first_1 <= lfi_first_2)) {
    if ((ffi_second_1 >= ffi_second_2 and ffi_second_1 <= lfi_second_2) or
        (lfi_second_1 >= ffi_second_2 and lfi_second_1 <= lfi_second_2)) {
      return true;
    }
  }
  return false;
}

} /* end of namespace consus */

template <class T>
struct HistInfo2D{
  /// logarithm of sum of all entries
  T log_length;
  /// first filled index for first dimension
  int ffi_first = std::numeric_limits<int>::max();
  /// last filled index for first dimension
  int lfi_first = std::numeric_limits<int>::lowest();
  /// first filled index for second dimension
  int ffi_second = std::numeric_limits<int>::max();
  /// last filled index for second dimension
  int lfi_second = std::numeric_limits<int>::lowest();
};

template <class T = double>
std::tuple<DiscreteAxis2D<T>, vec1<HistInfo2D<T>>, vec1<DiscreteAxis2D<T>>>
make_histogram2d(const vec2<T>& Timeseries, int col1, int col2,
                 const Range<T>& range1, const Range<T>& range2) {
  DiscreteAxis2D<T> Histogram(range1, range2, 0.0);
  vec1<DiscreteAxis2D<T>> MicroMeans(Timeseries.size(),
                                     DiscreteAxis2D<T>(range1, range2));
  vec1<HistInfo2D<T>> HistInfos;
  int ffi_first = std::numeric_limits<int>::max();
  int lfi_first = std::numeric_limits<int>::lowest();
  int ffi_second = std::numeric_limits<int>::max();
  int lfi_second = std::numeric_limits<int>::lowest();
  for (size_t j = 0; j < Timeseries[col1].size(); ++j) {
    const double E1 = Timeseries[col1][j];
    const double E2 = Timeseries[col2][j];
    ffi_first = std::min(ffi_first, Histogram.get_bin_first(E1));
    lfi_first = std::max(lfi_first, Histogram.get_bin_first(E1));
    ffi_second = std::min(ffi_second, Histogram.get_bin_second(E2));
    lfi_second = std::max(lfi_second, Histogram.get_bin_second(E2));
    int index = Histogram.get_Bin(E1, E2);
    Histogram[index] += 1.0;
    for (size_t k = 0; k < Timeseries.size(); ++k) {
      MicroMeans[k][index] += Timeseries[k][j];
    }
  }
  HistInfos.push_back({std::log(Timeseries[col1].size()), ffi_first, lfi_first,
                       ffi_second, lfi_second});
  return std::make_tuple(Histogram, HistInfos, MicroMeans);
}

template <class T = double>
void add_histogram2d(const vec2<T>& Timeseries, int col1, int col2,
                     DiscreteAxis2D<T>& Histogram, HistInfo2D<T>& HistInfo,
                     vec1<HistInfo2D<T>>& HistInfos,
                     vec1<DiscreteAxis2D<T>>& MicroMeans,
                     vec2<bool>& Overlaps) {
  int ffi_first = std::numeric_limits<int>::max();
  int lfi_first = std::numeric_limits<int>::lowest();
  int ffi_second = std::numeric_limits<int>::max();
  int lfi_second = std::numeric_limits<int>::lowest();
  for (size_t j = 0; j < Timeseries[col1].size(); ++j) {
    const double E1 = Timeseries[col1][j];
    const double E2 = Timeseries[col2][j];
    ffi_first = std::min(ffi_first, Histogram.get_bin_first(E1));
    lfi_first = std::max(lfi_first, Histogram.get_bin_first(E1));
    ffi_second = std::min(ffi_second, Histogram.get_bin_second(E2));
    lfi_second = std::max(lfi_second, Histogram.get_bin_second(E2));
    int index = Histogram.get_Bin(E1, E2);
    Histogram[index] += 1.0;
    for (size_t k = 0; k < Timeseries.size(); ++k) {
      MicroMeans[k][index] += Timeseries[k][j];
    }
  }
  vec1<bool> ThisOverlap;
  for (auto& hist: HistInfos){
    if (detail::overlaps(ffi_first, hist.ffi_first, lfi_first, hist.lfi_first,
                         ffi_second, hist.ffi_second, lfi_second,
                         hist.lfi_second)) {
      ThisOverlap.push_back(true);
    }else{
      ThisOverlap.push_back(false);
    }
  }
  Overlaps.push_back(ThisOverlap);
  HistInfo.ffi_first = std::min(ffi_first, HistInfo.ffi_first);
  HistInfo.lfi_first = std::max(lfi_first, HistInfo.lfi_first);
  HistInfo.ffi_second = std::min(ffi_second, HistInfo.ffi_second);
  HistInfo.lfi_second = std::max(lfi_second, HistInfo.lfi_second);
  HistInfo.log_length =
      addlogwise(HistInfo.log_length, std::log(Timeseries[col1].size()));
  HistInfos.push_back({std::log(Timeseries[col1].size()), ffi_first, lfi_first,
                       ffi_second, lfi_second});

}

// normalizes micronanoical mean values
template <class T = double>
void normalize_MicroMeans(DiscreteAxis2D<T>& Histogram,
                          vec1<DiscreteAxis2D<T>>& MicroMeans) {
  for (size_t i = 0; i < MicroMeans.size(); ++i){
    for (int j = 0; j < MicroMeans[i].get_num_bins_first(); ++j){
      for (int k = 0; k < MicroMeans[i].get_num_bins_second(); ++k){
        int index = MicroMeans[i].get_index(j, k);
        if (Histogram[index] != 0){
          MicroMeans[i][index] /= Histogram[index];
        }
      }
    }
  }
}

/// preliminary lnZ for Parameter_Target from histogram
/// measured at Parameter_Target
///
/// Hist - 2d histogram measured at Parameter
/// HistInfo - info object for Hist
/// Parameter - parameters for Hist
/// Parameter_Target - lnZ is estimated for this parameters
template <class T = double>
T estimate_lnZ_from_hist(const DiscreteAxis2D<T>& Hist,
                       const HistInfo2D<T>& HistInfo,
                       const std::pair<T, T>& Parameter,
                       const std::pair<T, T>& Parameter_Target) {
  auto DeltaParam_1 = Parameter_Target.first - Parameter.first;
  auto DeltaParam_2 = Parameter_Target.second - Parameter.second;
  auto log_Z_ratio = log_zero<T>();
  for (int i = HistInfo.ffi_first; i < HistInfo.lfi_first; ++i) {
    auto E1 = Hist.get_value_first(i);
    for (int j = HistInfo.ffi_second; j < Hist.lfi_second; ++j) {
      auto E2 = Hist.get_value_second(j);
      auto index = Hist.get_index(i, j);
      if (Hist[index] != T(0.0)) {
        log_Z_ratio = addlogwise(
            log_Z_ratio,
            std::log(Hist[index]) - (DeltaParam_1) * (E1 + (DeltaParam_2)*E2));
      }
    }
  }
  return log_Z_ratio;
}

template <class T = double>
void normalize_lnZ(vec1<T>& lnZ){
  for (size_t i = 1; i < lnZ.size(); ++i) {
    lnZ[i] -= lnZ[0];
  }
  lnZ[0] = 0;
}

template <class T = double>
void normalize_logDOS(DiscreteAxis2D<T>& logDOS, const vec1<T>& lnZ) {
  for (size_t i = 0; i < logDOS.size(); ++i) {
    logDOS[i] -= lnZ[0];
  }
}

template <class T = double>
T calc_lnZ(const DiscreteAxis2D<T>& DOS, const HistInfo2D<T>& HistInfo,
                 const std::pair<T, T>& Parameter) {
  T lnZ = log_zero<T>();
    for (int i = HistInfo.ffi_first; i < HistInfo.lfi_first; ++i) {
      auto E1 = DOS.get_value_first(i);
      for (int j = HistInfo.ffi_second; j < HistInfo.lfi_second; ++j) {
        auto E2 = DOS.get_value_second(j);
        int index = DOS.get_index(i, j);
        lnZ = addlogwise(
            lnZ, DOS[index] - Parameter.first * (E1 + Parameter.second * E2));
      }
    }
  return lnZ;
}

/// TODO: use HistInfos for each parameter instead of whole area from HistInfo
template <class T = double>
vec1<T> calc_lnZ(const DiscreteAxis2D<T>& DOS, const HistInfo2D<T>& HistInfo,
                 const vec1<std::pair<T, T>>& Parameters) {
  vec1<T> lnZ(Parameters.size(), log_zero<T>());
#pragma omp parallel for
  for (size_t k = 0; k < Parameters.size(); ++k) {
    for (int i = HistInfo.ffi_first; i < HistInfo.lfi_first; ++i) {
      auto E1 = DOS.get_value_first(i);
      for (int j = HistInfo.ffi_second; j < HistInfo.lfi_second; ++j) {
        auto E2 = DOS.get_value_second(j);
        int index = DOS.get_index(i, j);
        lnZ[k] = addlogwise(
            lnZ[k], DOS[index] -
                        Parameters[k].first * (E1 + Parameters[k].second * E2));
      }
    }
  }
  normalize_lnZ(lnZ);
  return lnZ;
}

template <class T = double>
vec1<T> calc_lnZ(const DiscreteAxis2D<T>& DOS,
                 const vec1<HistInfo2D<T>>& HistInfos,
                 const vec1<std::pair<T, T>>& Parameters) {
  vec1<T> lnZ(Parameters.size(), log_zero<T>());
#pragma omp parallel for
  for (size_t k = 0; k < Parameters.size(); ++k) {
    for (int i = HistInfos[k].ffi_first; i < HistInfos[k].lfi_first; ++i) {
      auto E1 = DOS.get_value_first(i);
      for (int j = HistInfos[k].ffi_second; j < HistInfos[k].lfi_second; ++j) {
        auto E2 = DOS.get_value_second(j);
        int index = DOS.get_index(i, j);
        lnZ[k] = addlogwise(
            lnZ[k], DOS[index] -
                        Parameters[k].first * (E1 + Parameters[k].second * E2));
      }
    }
  }
  normalize_lnZ(lnZ);
  return lnZ;
}

template <class T = double>
T deviation(const vec1<T>& lnZ_old, const vec1<T>& lnZ_new) {
  auto dev = T(0.0);
  for (size_t i = 0; i < lnZ_old.size(); ++i) {
    T tmp = (lnZ_new[i] - lnZ_old[i]) / lnZ_new[i];
    dev += tmp * tmp;
  }
  dev = std::sqrt(dev) / lnZ_old.size();
  return dev;
}

template <class T = double>
void iterate_logDOS(const DiscreteAxis2D<T>& Histogram,
                 const HistInfo2D<T>& HistInfo,
                 const vec1<HistInfo2D<T>>& HistInfos,
                 const vec1<std::pair<T, T>>& Parameters, const vec1<T>& lnZ,
                 DiscreteAxis2D<T>& logDOS) {
#pragma omp parallel for
  for (auto i = HistInfo.ffi_first; i < HistInfo.lfi_first; ++i) {
    auto E1 = logDOS.get_value_first(i);
    for (auto j = HistInfo.ffi_second; j < HistInfo.lfi_second; ++j) {
      auto E2 = logDOS.get_value_second(j);
      auto index = logDOS.get_index(i, j);
      auto tmp = log_zero<T>();
      for (size_t k = 0; k < lnZ.size(); ++k) {
        tmp = addlogwise(
            tmp, HistInfos[k].log_length -
                     Parameters[k].first * (E1 + Parameters[k].second * E2) -
                     lnZ[k]);
      }
      if (Histogram[index] != 0) {
        logDOS[index] = std::log(Histogram[index]) - tmp;
      } else {
        logDOS[index] = log_zero<T>() - tmp;
      }
    }
  }
}

template <class T = double>
void calc_logDOS_full(const DiscreteAxis2D<T>& Histogram,
                      const HistInfo2D<T>& HistInfo,
                      const vec1<HistInfo2D<T>>& HistInfos,
                      const vec1<std::pair<T, T>>& Parameters,
                      const double devmax, vec1<T>& lnZ,
                      DiscreteAxis2D<T>& logDOS) {
  T dev;
  do{
    iterate_logDOS(Histogram, HistInfo, HistInfos,  Parameters, devmax, lnZ, logDOS);
    auto new_lnZ = calc_lnZ(logDOS, HistInfo, Parameters);
    dev = deviation(lnZ, new_lnZ);
    std::swap(lnZ, new_lnZ);
    std::cout << dev << " " << devmax << "\n";
  } while(dev > devmax);
  normalize_logDOS(logDOS, lnZ);
  return logDOS;
}

template <class T = double>
void calc_logDOS_reduced(const DiscreteAxis2D<T>& Histogram,
                         const HistInfo2D<T>& HistInfo,
                         const vec1<HistInfo2D<T>>& HistInfos,
                         const vec1<std::pair<T, T>>& Parameters,
                         const double devmax, vec1<T>& lnZ,
                         DiscreteAxis2D<T>& logDOS) {
  T dev;
  do{
    iterate_logDOS(Histogram, HistInfo, HistInfos,  Parameters, devmax, lnZ, logDOS);
    auto new_lnZ = calc_lnZ(logDOS, HistInfos, Parameters);
    dev = deviation(lnZ, new_lnZ);
    std::swap(lnZ, new_lnZ);
    std::cout << dev << " " << devmax << "\n";
  } while(dev > devmax);
  normalize_logDOS(logDOS, lnZ);
  return logDOS;
}

/// call with first timeseries
///
/// returns in this order:
template <class T = double>
void logDOS_iteration_start(const vec2<T>& Timeseries, const int col1,
                            const int col2, const Range<T>& range1,
                            const Range<T>& range2, DiscreteAxis2D<T>& Hist,
                            vec1<DiscreteAxis2D<T>>& MicroMeans,
                            HistInfo2D<T>& HistInfo,
                            vec1<HistInfo2D<T>>& HistInfos, vec1<T>& lnZ) {
  assert(lnZ.size() == 1);
  std::tie(Hist, HistInfos, MicroMeans) =
      make_histogram2d(Timeseries, col1, col2, range1, range2);
  assert(HistInfos.size() == 1);
  HistInfo = HistInfos[0];
}

template <class T = double>
void logDOS_iteration_second(const vec2<T>& Timeseries, const int col1,
                             const int col2,
                             const vec1<std::pair<T, T>>& Parameters,
                             DiscreteAxis2D<T>& Hist,
                             vec1<DiscreteAxis2D<T>>& MicroMeans,
                             HistInfo2D<T>& HistInfo,
                             vec1<HistInfo2D<T>>& HistInfos, vec1<T>& lnZ,
                             DiscreteAxis2D<T>& logDOS, const T devmax) {
  assert(Parameters.size() == 2);
  assert(lnZ.size() == 1);
  assert(HistInfos.size() == 1);
  vec2<bool> Overlaps;
  lnZ.push_back(lnZ[0] + estimate_Z_from_Hist(Hist, HistInfo, Parameters[0],
                                              Parameters[1]));
  add_histogram2d(Timeseries, col1, col2, Hist, HistInfo, HistInfos, MicroMeans,
                  Overlaps);
  assert(lnZ.size() == 2);
  assert(HistInfos.size() == 2);
  assert(Overlaps.size() == 1);
  calc_logDOS_reduced(Hist, HistInfo, HistInfos, Parameters, devmax, lnZ, logDOS);
}

template <class T = double>
void logDOS_iteration_next(
    const vec2<T>& Timeseries, const int col1, const int col2,
    const vec1<std::pair<T, T>>& Parameters, DiscreteAxis2D<T>& Hist,
    vec1<DiscreteAxis2D<T>>& MicroMeans, HistInfo2D<T>& HistInfo,
    vec1<HistInfo2D<T>>& HistInfos, vec1<T>& lnZ, DiscreteAxis2D<T>& logDOS,
    vec2<bool>& Overlaps, const T devmax) {
  assert(Parameters.size() = lnZ.size() + 1);
  assert(lnZ.size() > 1);
  assert(HistInfos.size() > 1);
  lnZ.push_back(calc_lnZ(logDOS, HistInfo, Parameters.back()));
  add_histogram2d(Timeseries, col1, col2, Hist, HistInfo, HistInfos, MicroMeans,
                  Overlaps);
  assert(lnZ.size() > 2);
  assert(HistInfos.size() > 2);
  assert(Overlaps.size() > 1);
  calc_logDOS_reduced(Hist, HistInfo, HistInfos, Parameters, devmax, lnZ, logDOS);
  return std::make_tuple(Hist, MicroMeans, HistInfo, HistInfos, lnZ, logDOS, Overlaps);
}

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

//double reweight_dT(const DiscreteAxis2D& DOS, const DiscreteAxis2D& MicroMean,
//                   const double Beta, const double Parameter) {
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
//  double lnZ = LOGZERO;
//  for (int i = 0; i < DOS.get_NumBins_first(); ++i){
//    const double E1 = DOS.get_Value_first(i);
//    for (int j = 0; j < DOS.get_NumBins_second(); ++j){
//      const double E2 = DOS.get_Value_second(j);
//      const int index = DOS.get_index(i, j);
//      if (std::isfinite(DOS[index])){
//        lnZ = addlogwise(lnZ, DOS[index] - Beta * (E1 + Parameter * E2));
//      }
//    }
//  }
//
//  double O = LOGZERO;
//  double OE = LOGZERO;
//  double E = LOGZERO;
//  for (int i = 0; i < DOS.get_NumBins_first(); ++i){
//    const double e1 = DOS.get_Value_first(i);
//    for (int j = 0; j < DOS.get_NumBins_second(); ++j){
//      const double e2 = DOS.get_Value_second(j);
//      const int index = DOS.get_index(i, j);
//      if (std::isfinite(DOS[index]) and std::isfinite(MicroMean[index])) {
//        double tmp = DOS[index] - Beta * (e1 + Parameter * e2) - lnZ;
//        O = addlogwise(O, tmp + std::log(MicroMean[index] - shiftO));
//        OE = addlogwise(
//            OE,
//            tmp + std::log((e1 + Parameter * e2) * MicroMean[index] - shiftOE));
//        E = addlogwise(E, tmp + std::log((e1 + Parameter * e2) - shiftE));
//      }
//    }
//  }
//  O = std::exp(O) + shiftO;
//  OE = std::exp(OE) + shiftOE;
//  E = std::exp(E) + shiftE;
//  return Beta * Beta * ( OE - O * E );
//}
//
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
//
//double reweight_dT2(const DiscreteAxis2D& DOS, const DiscreteAxis2D& MicroMean,
//                    const double Beta, const double Parameter) {
//  double min = std::numeric_limits<double>::max();
//  double minE = std::numeric_limits<double>::max();
//  double minEE = std::numeric_limits<double>::max();
//  double minOE = std::numeric_limits<double>::max();
//  double minOEE = std::numeric_limits<double>::max();
//  for (int i = 0; i < DOS.get_NumBins_first(); ++i) {
//    const double E1 = DOS.get_Value_first(i);
//    for (int j = 0; j < DOS.get_NumBins_second(); ++j) {
//      const double E2 = DOS.get_Value_second(j);
//      const int index = DOS.get_index(i, j);
//      if (std::isfinite(MicroMean[index])) {
//        double tmpE = (E1 + Parameter * E2);
//        min = std::min(min, MicroMean[index]);
//        minE = std::min(minE, tmpE);
//        minEE = std::min(minEE, tmpE * tmpE);
//        minOE = std::min(minOE, MicroMean[index] * tmpE);
//        minOEE = std::min(minOEE, MicroMean[index] * tmpE * tmpE);
//      }
//    }
//  }
//  double shiftO = min - 0.1;
//  double shiftE = minE - 0.1;
//  double shiftEE = minEE - 0.1;
//  double shiftOE = minOE - 0.1;
//  double shiftOEE = minOEE - 0.1;
//
//  double lnZ = LOGZERO;
//  for (int i = 0; i < DOS.get_NumBins_first(); ++i){
//    const double E1 = DOS.get_Value_first(i);
//    for (int j = 0; j < DOS.get_NumBins_second(); ++j){
//      const double E2 = DOS.get_Value_second(j);
//      const int index = DOS.get_index(i, j);
//      if (std::isfinite(DOS[index])){
//        lnZ = addlogwise(lnZ, DOS[index] - Beta * (E1 + Parameter * E2));
//      }
//    }
//  }
//
//  double O = LOGZERO;
//  double OE = LOGZERO;
//  double OEE = LOGZERO;
//  double E = LOGZERO;
//  double EE = LOGZERO;
//  for (int i = 0; i < DOS.get_NumBins_first(); ++i){
//    const double e1 = DOS.get_Value_first(i);
//    for (int j = 0; j < DOS.get_NumBins_second(); ++j){
//      const double e2 = DOS.get_Value_second(j);
//      const int index = DOS.get_index(i, j);
//      if (std::isfinite(DOS[index]) and std::isfinite(MicroMean[index])) {
//        double tmpE = (e1 + Parameter * e2);
//        double tmp = DOS[index] - Beta * (tmpE) - lnZ;
//        O = addlogwise(O, tmp + std::log(MicroMean[index] - shiftO));
//        OE = addlogwise(OE, tmp + std::log(tmpE * MicroMean[index] - shiftOE));
//        OEE = addlogwise(
//            OEE, tmp + std::log(tmpE * tmpE * MicroMean[index] - shiftOEE));
//        E = addlogwise(E, tmp + std::log(tmpE - shiftE));
//        EE = addlogwise(EE, tmp + std::log(tmpE * tmpE - shiftEE));
//      }
//    }
//  }
//  O = std::exp(O) + shiftO;
//  OE = std::exp(OE) + shiftOE;
//  OEE = std::exp(OEE) + shiftOEE;
//  E = std::exp(E) + shiftE;
//  EE = std::exp(EE) + shiftEE;
//  return 2 * Beta * Beta * Beta * (O * E - OE) +
//         Beta * Beta * Beta * Beta * (OEE - O * EE - 2 * OE * E + 2 * O * E * E);
//}
//
//vec1d reweight_dT2(const DiscreteAxis2D& DOS, const DiscreteAxis2D& MicroMean,
//                   const vec1d& Betas, const double Parameter) {
//  double min = std::numeric_limits<double>::max();
//  double minE = std::numeric_limits<double>::max();
//  double minEE = std::numeric_limits<double>::max();
//  double minOE = std::numeric_limits<double>::max();
//  double minOEE = std::numeric_limits<double>::max();
//  for (int i = 0; i < DOS.get_NumBins_first(); ++i) {
//    const double E1 = DOS.get_Value_first(i);
//    for (int j = 0; j < DOS.get_NumBins_second(); ++j) {
//      const double E2 = DOS.get_Value_second(j);
//      const int index = DOS.get_index(i, j);
//      if (std::isfinite(MicroMean[index])) {
//        min = std::min(min, MicroMean[index]);
//        double tmpE = (E1 + Parameter * E2);
//        minE = std::min(minE, tmpE);
//        minEE = std::min(minEE, tmpE * tmpE);
//        minOE = std::min(minOE, MicroMean[index] * tmpE);
//        minOEE = std::min(minOEE, MicroMean[index] * tmpE * tmpE);
//      }
//    }
//  }
//  double shiftO = min - 0.1;
//  double shiftE = minE - 0.1;
//  double shiftEE = minEE - 0.1;
//  double shiftOE = minOE - 0.1;
//  double shiftOEE = minOEE - 0.1;
//
//  vec1d result(Betas.size());
//
//  for (size_t k = 0; k < Betas.size(); ++k){
//    double lnZ = LOGZERO;
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
//    double O = LOGZERO;
//    double OE = LOGZERO;
//    double OEE = LOGZERO;
//    double E = LOGZERO;
//    double EE = LOGZERO;
//    for (int i = 0; i < DOS.get_NumBins_first(); ++i){
//      const double e1 = DOS.get_Value_first(i);
//      for (int j = 0; j < DOS.get_NumBins_second(); ++j){
//        const double e2 = DOS.get_Value_second(j);
//        const int index = DOS.get_index(i, j);
//        if (std::isfinite(DOS[index]) and std::isfinite(MicroMean[index])) {
//          double tmpE = (e1 + Parameter * e2);
//          double tmp = DOS[index] - Betas[k] * tmpE - lnZ;
//          O = addlogwise(O, tmp + std::log(MicroMean[index] - shiftO));
//          OE =
//              addlogwise(OE, tmp + std::log(tmpE * MicroMean[index] - shiftOE));
//          OEE = addlogwise(
//              OEE, tmp + std::log(tmpE * tmpE * MicroMean[index] - shiftOEE));
//          E = addlogwise(E, tmp + std::log(tmpE - shiftE));
//          EE = addlogwise(EE, tmp + std::log(tmpE * tmpE - shiftEE));
//        }
//      }
//    }
//    O = std::exp(O) + shiftO;
//    OE = std::exp(OE) + shiftOE;
//    OEE = std::exp(OEE) + shiftOEE;
//    E = std::exp(E) + shiftE;
//    EE = std::exp(EE) + shiftEE;
//    double b = Betas[k];
//    result[k] = 2 * b * b * b * (O * E - OE) +
//                b * b * b * b * (OEE - O * EE - 2 * OE * E + 2 * O * E * E);
//  }
//  return result;
//}
//
//DiscreteAxis2D reweight_hist2d(const DiscreteAxis2D& DOS, const double Beta,
//                               const double Parameter) {
//  DiscreteAxis2D result(DOS.get_Range_first(), DOS.get_Range_second());
//  std::fill(result.begin(), result.end(), 0.0);
//  double lnZ = LOGZERO;
//  for (int i = 0; i < DOS.get_NumBins_first(); ++i){
//    const double E1 = DOS.get_Value_first(i);
//    for (int j = 0; j < DOS.get_NumBins_second(); ++j){
//      const double E2 = DOS.get_Value_second(j);
//      const int index = DOS.get_index(i, j);
//      if (std::isfinite(DOS[index])){
//        result[index] = DOS[index] - Beta * ( E1 + Parameter * E2 );
//        lnZ = addlogwise(lnZ, DOS[index] - Beta * ( E1 + Parameter * E2));
//      }
//    }
//  }
//  for (size_t i = 0; i < result.size(); ++i){
//    result[i] -= lnZ;
//    result[i] = std::exp(result[i]);
//  }
//  return result;
//}
//
//DiscreteAxis2D reweight_probdens2d(const DiscreteAxis2D& DOS, const double Beta,
//                                   const double Parameter) {
//  DiscreteAxis2D result(DOS.get_Range_first(), DOS.get_Range_second());
//  std::fill(result.begin(), result.end(), 0.0);
//  double lnZ = LOGZERO;
//  for (int i = 0; i < DOS.get_NumBins_first(); ++i){
//    const double E1 = DOS.get_Value_first(i);
//    for (int j = 0; j < DOS.get_NumBins_second(); ++j){
//      const double E2 = DOS.get_Value_second(j);
//      const int index = DOS.get_index(i, j);
//      if (std::isfinite(DOS[index])){
//        result[index] = DOS[index] - Beta * ( E1 + Parameter * E2 );
//        lnZ = addlogwise(lnZ, DOS[index] - Beta * ( E1 + Parameter * E2));
//      }
//    }
//  }
//  const double area = DOS.get_step_first() * DOS.get_step_second();
//  for (size_t i = 0; i < result.size(); ++i){
//    result[i] -= lnZ;
//    result[i] = std::exp(result[i])/area;
//  }
//  return result;
//}
//
//DiscreteAxis reweight_hist_obs(const DiscreteAxis2D& DOS,
//                               const DiscreteAxis2D& MicroMean, const double NumBins,
//                               const double Beta, const double Parameter) {
//  double min = std::numeric_limits<double>::max();
//  double max = std::numeric_limits<double>::min();
//  for (size_t i = 0; i < MicroMean.size(); ++i) {
//    if (std::isfinite(MicroMean[i])) {
//      min = std::min(min, MicroMean[i]);
//      max = std::max(max, MicroMean[i]);
//    }
//  }
//  //double shift = min - 0.01;
//
//  double step = (max - min)/NumBins;
//
//  DiscreteAxis result(min, step, max, 0.0);
//  double lnZ = LOGZERO;
//  for (int i = 0; i < DOS.get_NumBins_first(); ++i){
//    const double E1 = DOS.get_Value_first(i);
//    for (int j = 0; j < DOS.get_NumBins_second(); ++j){
//      const double E2 = DOS.get_Value_second(j);
//      const int index = DOS.get_index(i, j);
//      const int rindex = result.get_Bin(MicroMean[index]);
//      if (std::isfinite(DOS[index])){
//        result[rindex] += std::exp(DOS[index]) *
//                          std::exp(-1.0 * Beta * (E1 + Parameter * E2));
//        lnZ = addlogwise(lnZ, DOS[index] - Beta * ( E1 + Parameter * E2));
//      }
//    }
//  }
//  const double Z = std::exp(lnZ);
//  for (size_t i = 0; i < result.size(); ++i){
//    result[i] /= Z;
//  }
//  return result;
//}
//
//double reweight_prob(const DiscreteAxis2D& DOS, const DiscreteAxis2D& MicroMean,
//                     const double Beta, const double Parameter,
//                     const double min, const double max) {
//  double result = 0;
//  double lnZ = LOGZERO;
//  for (int i = 0; i < DOS.get_NumBins_first(); ++i){
//    const double E1 = DOS.get_Value_first(i);
//    for (int j = 0; j < DOS.get_NumBins_second(); ++j){
//      const double E2 = DOS.get_Value_second(j);
//      const int index = DOS.get_index(i, j);
//      if (std::isfinite(DOS[index])){
//        if (MicroMean[index] >= min and MicroMean[index] < max){
//        result += std::exp(DOS[index]) *
//                          std::exp(-1.0 * Beta * (E1 + Parameter * E2));
//        }
//        lnZ = addlogwise(lnZ, DOS[index] - Beta * ( E1 + Parameter * E2));
//      }
//    }
//  }
//  result /= std::exp(lnZ);
//  return result;
//}
//
//DiscreteAxis reweight_hist(const DiscreteAxis2D& DOS, const double step,
//                           const double Beta, const double Parameter) {
//  double minE = DOS.get_min_first() + Parameter * DOS.get_min_second();
//  double maxE = DOS.get_max_first() + Parameter * DOS.get_max_second();
//  DiscreteAxis result(minE, step, maxE, LOGZERO);
//  double lnZ = LOGZERO;
//  for (int i = 0; i < DOS.get_NumBins_first(); ++i){
//    const double E1 = DOS.get_Value_first(i);
//    for (int j = 0; j < DOS.get_NumBins_second(); ++j){
//      const double E2 = DOS.get_Value_second(j);
//      const double E = E1 + Parameter * E2;
//      const int index = DOS.get_index(i, j);
//      const int rindex = result.get_Bin(E);
//      if (std::isfinite(DOS[index])){
//        result[rindex] = addlogwise(result[rindex], DOS[index] - Beta * E);
//        lnZ = addlogwise(lnZ, DOS[index] - Beta * E);
//      }
//    }
//  }
//  for (size_t i = 0; i < result.size(); ++i){
//      result[i] -= lnZ;
//      result[i] = std::exp(result[i]);
//  }
//  return result;
//}
//
//DiscreteAxis reweight_probdens(const DiscreteAxis2D& DOS, const double step,
//                           const double Beta, const double Parameter) {
//  double minE = DOS.get_min_first() + Parameter * DOS.get_min_second();
//  double maxE = DOS.get_max_first() + Parameter * DOS.get_max_second();
//  DiscreteAxis result(minE, step, maxE, LOGZERO);
//  double lnZ = LOGZERO;
//  for (int i = 0; i < DOS.get_NumBins_first(); ++i){
//    const double E1 = DOS.get_Value_first(i);
//    for (int j = 0; j < DOS.get_NumBins_second(); ++j){
//      const double E2 = DOS.get_Value_second(j);
//      const double E = E1 + Parameter * E2;
//      const int index = DOS.get_index(i, j);
//      const int rindex = result.get_Bin(E);
//      if (std::isfinite(DOS[index])){
//        result[rindex] = addlogwise(result[rindex], DOS[index] - Beta * E);
//        lnZ = addlogwise(lnZ, DOS[index] - Beta * E);
//      }
//    }
//  }
//  for (size_t i = 0; i < result.size(); ++i){
//      result[i] -= lnZ;
//      result[i] = std::exp(result[i])*step;
//  }
//  return result;
//}
//
//double reweight_dP(const DiscreteAxis2D& DOS, const DiscreteAxis2D& MicroMean,
//                   const double Beta, const double Parameter) {
//  double min = std::numeric_limits<double>::max();
//  double minE = std::numeric_limits<double>::max();
//  double minOE = std::numeric_limits<double>::max();
//  for (int i = 0; i < DOS.get_NumBins_first(); ++i) {
//    const double E1 = DOS.get_Value_first(i);
//    for (int j = 0; j < DOS.get_NumBins_second(); ++j) {
//      const int index = DOS.get_index(i, j);
//      if (std::isfinite(MicroMean[index])) {
//        min = std::min(min, MicroMean[index]);
//        minE = std::min(minE, E1);
//        minOE = std::min(minOE, MicroMean[index] * (E1));
//      }
//    }
//  }
//  double shiftO = min - 0.1;
//  double shiftE = minE - 0.1;
//  double shiftOE = minOE - 0.1;
//
//  double lnZ = LOGZERO;
//  for (int i = 0; i < DOS.get_NumBins_first(); ++i) {
//    const double E1 = DOS.get_Value_first(i);
//    for (int j = 0; j < DOS.get_NumBins_second(); ++j) {
//      const double E2 = DOS.get_Value_second(j);
//      const int index = DOS.get_index(i, j);
//      if (std::isfinite(DOS[index])) {
//        lnZ = addlogwise(lnZ, DOS[index] - Beta * (E1 + Parameter * E2));
//      }
//    }
//  }
//
//  double O = LOGZERO;
//  double OE = LOGZERO;
//  double E = LOGZERO;
//  for (int i = 0; i < DOS.get_NumBins_first(); ++i) {
//    const double e1 = DOS.get_Value_first(i);
//    for (int j = 0; j < DOS.get_NumBins_second(); ++j) {
//      const double e2 = DOS.get_Value_second(j);
//      const int index = DOS.get_index(i, j);
//      if (std::isfinite(DOS[index]) and std::isfinite(MicroMean[index])) {
//        double tmp = DOS[index] - Beta * (e1 + Parameter * e2) - lnZ;
//        O = addlogwise(O, tmp + std::log(MicroMean[index] - shiftO));
//        OE = addlogwise(OE, tmp + std::log((e1) * MicroMean[index] - shiftOE));
//        E = addlogwise(E, tmp + std::log((e1) - shiftE));
//      }
//    }
//  }
//  O = std::exp(O) + shiftO;
//  OE = std::exp(OE) + shiftOE;
//  E = std::exp(E) + shiftE;
//  return Beta * (OE - O * E);
//}
//
//vec1d reweight_dP(const DiscreteAxis2D& DOS, const DiscreteAxis2D& MicroMean,
//                  const vec1d& Betas, const double Parameter) {
//  double min = std::numeric_limits<double>::max();
//  double minE = std::numeric_limits<double>::max();
//  double minOE = std::numeric_limits<double>::max();
//  for (int i = 0; i < DOS.get_NumBins_first(); ++i) {
//    const double E1 = DOS.get_Value_first(i);
//    for (int j = 0; j < DOS.get_NumBins_second(); ++j) {
//      const int index = DOS.get_index(i, j);
//      if (std::isfinite(MicroMean[index])) {
//        min = std::min(min, MicroMean[index]);
//        minE = std::min(minE, E1);
//        minOE = std::min(minOE, MicroMean[index] * (E1));
//      }
//    }
//  }
//  double shiftO = min - 0.1;
//  double shiftE = minE - 0.1;
//  double shiftOE = minOE - 0.1;
//
//  vec1d result(Betas.size());
//
//  for (size_t k = 0; k < Betas.size(); ++k) {
//    double lnZ = LOGZERO;
//    double O = LOGZERO;
//    double OE = LOGZERO;
//    double E = LOGZERO;
//    for (int i = 0; i < DOS.get_NumBins_first(); ++i) {
//      const double E1 = DOS.get_Value_first(i);
//      for (int j = 0; j < DOS.get_NumBins_second(); ++j) {
//        const double E2 = DOS.get_Value_second(j);
//        const int index = DOS.get_index(i, j);
//        if (std::isfinite(DOS[index])) {
//          lnZ = addlogwise(lnZ, DOS[index] - Betas[k] * (E1 + Parameter * E2));
//        }
//      }
//    }
//
//    for (int i = 0; i < DOS.get_NumBins_first(); ++i) {
//      const double e1 = DOS.get_Value_first(i);
//      for (int j = 0; j < DOS.get_NumBins_second(); ++j) {
//        const double e2 = DOS.get_Value_second(j);
//        const int index = DOS.get_index(i, j);
//        if (std::isfinite(DOS[index]) and std::isfinite(MicroMean[index])) {
//          double tmp = DOS[index] - Betas[k] * (e1 + Parameter * e2) - lnZ;
//          double tmpE = e1;
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
//    result[k] = Betas[k] * Betas[k] * (OE - O * E);
//  }
//  return result;
//}
//
//DiscreteAxis make_hist_obs(const vec2d& Timeseries, const int col,
//                           const Range& range) {
//  DiscreteAxis Histogram(range, 0.0);
//  for (size_t i = 0; i < Timeseries[col].size(); ++i) {
//    int index = Histogram.get_Bin(Timeseries[col][i]);
//    Histogram[index] += 1.0;
//  }
//  return Histogram;
//}
//
//DiscreteAxis add_hist_obs(DiscreteAxis& Hist, const vec2d& Timeseries,
//                          const int col) {
//  for (size_t i = 0; i < Timeseries[col].size(); ++i) {
//    int index = Hist.get_Bin(Timeseries[col][i]);
//    Hist[index] += 1.0;
//  }
//  return Hist;
//}
//
//DiscreteAxis2D make_hist_obs2d(const vec2d& Timeseries, const int col1,
//                               const int col2, const Range& range1,
//                               const Range& range2) {
//  DiscreteAxis2D Histogram(range1, range2);
//  std::fill(Histogram.begin(), Histogram.end(), 0.0);
//  for (size_t i = 0; i < Timeseries[col1].size(); ++i){
//    int index = Histogram.get_Bin( Timeseries[col1][i], Timeseries[col2][i] );
//    Histogram[index] += 1.0;
//  }
//  return Histogram;
//}
//
//DiscreteAxis2D add_hist_obs2d(DiscreteAxis2D& Hist, const vec2d& Timeseries,
//                              const int col1, const int col2) {
//  for (size_t i = 0; i < Timeseries[col1].size(); ++i) {
//    int index = Hist.get_Bin(Timeseries[col1][i], Timeseries[col2][i]);
//    Hist[index] += 1.0;
//  }
//  return Hist;
//}
//
} /* end of namespace WHAM2D */ 
} /* end of namespace consus */ 

