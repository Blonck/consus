#pragma once

#include <vector>
#include <tuple>
#include <algorithm>
#include <chrono>
#include "reweight.hpp"
#include "HistInfo2D.hpp"
#include "Ensemble.hpp"
#include "../DiscreteAxis2D.hpp"
#include "../vector/helper.hpp"
#include "../addlogwise.hpp"

namespace consus
{

/// everthing one need to calculate the 2 dimensional density of states dependent on two variables
namespace WHAM2D
{

namespace detail {

// checks if two rectangles (in the filled entries in histogram) overlaps
inline bool overlaps(const HistInfo2D& h1, const HistInfo2D& h2) {
  if ((h1.ffi_first >= h2.ffi_first and h1.ffi_first <= h2.lfi_first) or
      (h1.lfi_first >= h2.ffi_first and h1.lfi_first <= h2.lfi_first)) {
    if ((h1.ffi_second >= h2.ffi_second and h1.ffi_second <= h2.lfi_second) or
        (h1.lfi_second >= h2.ffi_second and h1.lfi_second <= h2.lfi_second)) {
      return true;
    }
  }
  return false;
}

} /* end of namespace detail */

#pragma omp declare reduction( +: DiscreteAxis2D : \
    std::transform(omp_in.begin(), omp_in.end(),   \
      omp_out.begin(), omp_out.begin(),            \
      std::plus<double>()) )                       \
      initializer (omp_priv(omp_orig))

std::tuple<DiscreteAxis2D, HistInfo2D>
make_histogram2d(const vec2<double>& Timeseries, int col1, int col2,
                 const Range& range1, const Range& range2) {
  DiscreteAxis2D Histogram(range1, range2, 0.0);
  vec1<HistInfo2D> HistInfos;
  int ffi_first = std::numeric_limits<int>::max();
  int lfi_first = std::numeric_limits<int>::lowest();
  int ffi_second = std::numeric_limits<int>::max();
  int lfi_second = std::numeric_limits<int>::lowest();
  //TODO: test reduction before use it
  //#pragma omp parallel for reduction(+: Histogram)
  //                         reduction(min: ffi_first, ffi_second)
  //                         reduction(max: lfi_first, lfi_second)
  for (size_t j = 0; j < Timeseries[col1].size(); ++j) {
    const double E1 = Timeseries[col1][j];
    const double E2 = Timeseries[col2][j];
    ffi_first = std::min(ffi_first, Histogram.get_bin_first(E1));
    lfi_first = std::max(lfi_first, Histogram.get_bin_first(E1));
    ffi_second = std::min(ffi_second, Histogram.get_bin_second(E2));
    lfi_second = std::max(lfi_second, Histogram.get_bin_second(E2));
    int index = Histogram.get_bin(E1, E2);
    Histogram[index] += 1.0;
  }
  //#pragma omp parallel for
  for (int i = 0; i < Histogram.size(); ++i){
    if (Histogram[i] != 0.0){
      Histogram[i] = std::log(Histogram[i]);
    }else{
      Histogram[i] = log_zero<double>();
    }
  }
  double log_length = std::log(Timeseries[col1].size());
  HistInfo2D HistInfo(log_length, ffi_first, lfi_first, ffi_second,
                          lfi_second);
  return std::make_tuple(Histogram, HistInfo);
}

void merge_histograms2d(DiscreteAxis2D& fullHist, HistInfo2D& fullHistInfo,
                      const DiscreteAxis2D& tmpHist, const HistInfo2D& tmpHistInfo){
  assert(fullHist.size() == tmpHist.size());
  assert(fullHist.get_min_first() == tmpHist.get_min_first());
  assert(fullHist.get_min_second() == tmpHist.get_min_second());
  assert(fullHist.get_max_first() == tmpHist.get_max_first());
  assert(fullHist.get_max_second() == tmpHist.get_max_second());
  assert(fullHist.get_step_first() == tmpHist.get_step_first());
  assert(fullHist.get_step_second() == tmpHist.get_step_second());
  fullHistInfo.ffi_first =
      std::min(fullHistInfo.ffi_first, tmpHistInfo.ffi_first);
  fullHistInfo.lfi_first =
      std::max(fullHistInfo.lfi_first, tmpHistInfo.lfi_first);
  fullHistInfo.ffi_second =
      std::min(fullHistInfo.ffi_second, tmpHistInfo.ffi_second);
  fullHistInfo.lfi_second =
      std::max(fullHistInfo.lfi_second, tmpHistInfo.lfi_second);
  fullHistInfo.log_length =
      addlogwise(fullHistInfo.log_length, tmpHistInfo.log_length);

  #pragma omp parallel for
  for (int i = 0; i < tmpHist.size(); ++i){
    if (tmpHist[i] != log_zero<double>()){
      fullHist[i] = addlogwise(fullHist[i], tmpHist[i]);
    }
  }
}

DiscreteAxis2D add_histogram2d(const vec2<double>& Timeseries, int col1,
                               int col2, DiscreteAxis2D& Histogram,
                               HistInfo2D& HistInfo,
                               vec1<HistInfo2D>& HistInfos) {
  int ffi_first = std::numeric_limits<int>::max();
  int lfi_first = std::numeric_limits<int>::lowest();
  int ffi_second = std::numeric_limits<int>::max();
  int lfi_second = std::numeric_limits<int>::lowest();
  DiscreteAxis2D tmpHist(Histogram, 0.0);
  assert(tmpHist.size() == Histogram.size());
  assert(tmpHist.get_min_first() == Histogram.get_min_first());
  assert(tmpHist.get_min_second() == Histogram.get_min_second());
  //TODO: test recuction before use it
  //#pragma omp parallel for reduction(+: tmpHist) reduction(min: ffi_first, ffi_second) reduction(max: lfi_first, lfi_second)
  for (size_t j = 0; j < Timeseries[col1].size(); ++j) {
    const double E1 = Timeseries[col1][j];
    const double E2 = Timeseries[col2][j];
    ffi_first = std::min(ffi_first, Histogram.get_bin_first(E1));
    lfi_first = std::max(lfi_first, Histogram.get_bin_first(E1));
    ffi_second = std::min(ffi_second, Histogram.get_bin_second(E2));
    lfi_second = std::max(lfi_second, Histogram.get_bin_second(E2));
    int index = tmpHist.get_bin(E1, E2);
    tmpHist[index] += 1.0;
  }

  #pragma omp parallel for
  for (int i = 0; i < tmpHist.size(); ++i){
    if (tmpHist[i] != 0){
      tmpHist[i] = std::log(tmpHist[i]);
      Histogram[i] = addlogwise(Histogram[i],tmpHist[i]);
    }
  }
#ifdef VERBOSE
  std::cout << "old Range [" << HistInfo.ffi_first << ", "
            << HistInfo.lfi_first << "] [" << HistInfo.ffi_second << ", "
            << HistInfo.lfi_second << "]\n";
#endif
  HistInfo.ffi_first = std::min(ffi_first, HistInfo.ffi_first);
  HistInfo.lfi_first = std::max(lfi_first, HistInfo.lfi_first);
  HistInfo.ffi_second = std::min(ffi_second, HistInfo.ffi_second);
  HistInfo.lfi_second = std::max(lfi_second, HistInfo.lfi_second);
  HistInfo.log_length =
      addlogwise(HistInfo.log_length, std::log(Timeseries[col1].size()));
  HistInfos.push_back({std::log(Timeseries[col1].size()), ffi_first, lfi_first,
                       ffi_second, lfi_second});
#ifdef VERBOSE
  std::cout << "add Range [" << ffi_first << ", " << lfi_first << "] ["
            << ffi_second << ", " << lfi_second << "]\n";
#endif
  return tmpHist;
}

/// preliminary lnZ for Parameter_Target from histogram
/// measured at Parameter_Target
///
/// Hist - 2d histogram measured at Parameter
/// HistInfo - info object for Hist
/// Parameter - parameters for Hist
/// Parameter_Target - lnZ is estimated for this parameters
template <class TEnsemble>
double estimate_lnZ_from_hist(const DiscreteAxis2D& Hist,
                         const HistInfo2D& HistInfo,
                         const std::pair<double, double>& Parameter,
                         const std::pair<double, double>& Parameter_Target) {
  auto log_Z_ratio = log_zero<double>();
  //TODO: maybe pragma with log_Z_ratio shared
  for (int i = HistInfo.ffi_first; i < HistInfo.lfi_first; ++i) {
    auto E1 = Hist.get_value_first(i);
    for (int j = HistInfo.ffi_second; j < HistInfo.lfi_second; ++j) {
      auto E2 = Hist.get_value_second(j);
      auto index = Hist.get_index(i, j);
      if (Hist[index] != log_zero<double>()) {
        log_Z_ratio = addlogwise(
            log_Z_ratio,
            Hist[index] + TEnsemble::log_weight(Parameter_Target.first, E1,
                                                Parameter_Target.second, E2) -
                TEnsemble::log_weight(Parameter.first, E1, Parameter.second,
                                      E2));
      }
    }
  }
  return log_Z_ratio;
}

void normalize_lnZ(vec1<double>& lnZ){
  for (size_t i = 1; i < lnZ.size(); ++i) {
    lnZ[i] -= lnZ[0];
  }
  lnZ[0] = 0;
}

void normalize_logDOS(DiscreteAxis2D& logDOS, const vec1<double>& lnZ) {
  assert(lnZ.size() > 0);
#pragma omp parallel for
  for (int i = 0; i < logDOS.size(); ++i) {
    logDOS[i] -= lnZ[0];
  }
}

template <class TEnsemble>
double calc_lnZ(const DiscreteAxis2D& DOS, const HistInfo2D& HistInfo,
                 const std::pair<double, double>& Parameter) {
  double lnZ = log_zero<double>();
    for (int i = HistInfo.ffi_first; i < HistInfo.lfi_first; ++i) {
      auto E1 = DOS.get_value_first(i);
      for (int j = HistInfo.ffi_second; j < HistInfo.lfi_second; ++j) {
        auto E2 = DOS.get_value_second(j);
        int index = DOS.get_index(i, j);
        lnZ = addlogwise(
            lnZ, DOS[index] + TEnsemble::log_weight(Parameter.first, E1,
                                                    Parameter.second, E2));
      }
    }
  return lnZ;
}

/// TODO: use HistInfos for each parameter instead of whole area from HistInfo
template <class TEnsemble>
vec1<double> calc_lnZ(const DiscreteAxis2D& DOS, const HistInfo2D& HistInfo,
                      const vec1<std::pair<double, double>>& Parameters) {
  vec1<double> lnZ(Parameters.size(), log_zero<double>());
#pragma omp parallel for
  for (size_t k = 0; k < Parameters.size(); ++k) {
    for (int i = HistInfo.ffi_first; i < HistInfo.lfi_first; ++i) {
      auto E1 = DOS.get_value_first(i);
      for (int j = HistInfo.ffi_second; j < HistInfo.lfi_second; ++j) {
        auto E2 = DOS.get_value_second(j);
        int index = DOS.get_index(i, j);
        lnZ[k] = addlogwise(lnZ[k], DOS[index] + TEnsemble::log_weight(
                                                     Parameters[k].first, E1,
                                                     Parameters[k].second, E2));
      }
    }
  }
  normalize_lnZ(lnZ);
  return lnZ;
}

template <class TEnsemble>
vec1<double> calc_lnZ_reduced(
    const DiscreteAxis2D& DOS, const vec1<HistInfo2D>& HistInfos,
    const vec1<std::pair<double, double>>& Parameters) {
  vec1<double> lnZ(Parameters.size(), log_zero<double>());
  assert(Parameters.size() == HistInfos.size());
#pragma omp parallel for
  for (size_t k = 0; k < Parameters.size(); ++k) {
    for (int i = HistInfos[k].ffi_first; i < HistInfos[k].lfi_first; ++i) {
      auto E1 = DOS.get_value_first(i);
      for (int j = HistInfos[k].ffi_second; j < HistInfos[k].lfi_second; ++j) {
        auto E2 = DOS.get_value_second(j);
        int index = DOS.get_index(i, j);
        lnZ[k] = addlogwise(lnZ[k], DOS[index] + TEnsemble::log_weight(
                                                     Parameters[k].first, E1,
                                                     Parameters[k].second, E2));
      }
    }
  }
  normalize_lnZ(lnZ);
  return lnZ;
}

template <class TEnsemble>
vec1<double> calc_lnZ_reduced(const DiscreteAxis2D& DOS,
                              const vec1<HistInfo2D>& HistInfos,
                              const vec1<std::pair<double, double>>& Parameters,
                              const vec1<int>& Overlap, const int treshold,
                              const vec1<double>& old_lnZ) {
  vec1<double> lnZ(Parameters.size(), log_zero<double>());
  assert(old_lnZ.size() == Parameters.size());
  assert(Overlap.size() > 0);
  assert(Parameters.size() == HistInfos.size());
  assert(Parameters.size()-1 == Overlap.size());
#pragma omp parallel for schedule(dynamic)
  for (size_t k = 0; k < Parameters.size()-1; ++k) {
    if (Overlap[k] >= treshold) {
      for (int i = HistInfos[k].ffi_first; i < HistInfos[k].lfi_first; ++i) {
        auto E1 = DOS.get_value_first(i);
        for (int j = HistInfos[k].ffi_second; j < HistInfos[k].lfi_second;
             ++j) {
          auto E2 = DOS.get_value_second(j);
          int index = DOS.get_index(i, j);
          lnZ[k] = addlogwise(
              lnZ[k],
              DOS[index] + TEnsemble::log_weight(Parameters[k].first, E1,
                                                 Parameters[k].second, E2));
        }
      }
    }else{
     lnZ[k] = old_lnZ[k]; 
    }
  }
  size_t k = Overlap.size();
  for (int i = HistInfos[k].ffi_first; i < HistInfos[k].lfi_first; ++i) {
    auto E1 = DOS.get_value_first(i);
    for (int j = HistInfos[k].ffi_second; j < HistInfos[k].lfi_second; ++j) {
      auto E2 = DOS.get_value_second(j);
      int index = DOS.get_index(i, j);
      lnZ[k] = addlogwise(
          lnZ[k], DOS[index] + TEnsemble::log_weight(Parameters[k].first, E1,
                                                     Parameters[k].second, E2));
    }
  }

  normalize_lnZ(lnZ);
  return lnZ;
}

/// assumption lnZ[0] = 0
double deviation(const vec1<double>& lnZ_old, const vec1<double>& lnZ_new) {
  assert(lnZ_old.size() == lnZ_new.size());
  assert(lnZ_new.size() > 0);
  auto dev = 0.0;
  for (size_t i = 1; i < lnZ_old.size(); ++i) {
    double tmp = (lnZ_new[i] - lnZ_old[i]) / lnZ_new[i];
    dev += tmp * tmp;
  }
  dev = std::sqrt(dev) / lnZ_old.size();
  return dev;
}

template <class TEnsemble>
void iterate_logDOS(const DiscreteAxis2D& Histogram, const HistInfo2D& HistInfo,
                    const vec1<HistInfo2D>& HistInfos,
                    const vec1<std::pair<double, double>>& Parameters,
                    const vec1<double>& lnZ, DiscreteAxis2D& logDOS) {
#pragma omp parallel for
  for (auto i = HistInfo.ffi_first; i < HistInfo.lfi_first; ++i) {
    auto E1 = logDOS.get_value_first(i);
    for (auto j = HistInfo.ffi_second; j < HistInfo.lfi_second; ++j) {
      auto E2 = logDOS.get_value_second(j);
      auto index = logDOS.get_index(i, j);
      auto tmp = log_zero<double>();
      for (size_t k = 0; k < lnZ.size(); ++k) {
        tmp = addlogwise(tmp,
                         HistInfos[k].log_length +
                             TEnsemble::log_weight(Parameters[k].first, E1,
                                                   Parameters[k].second, E2) -
                             lnZ[k]);
      }
      logDOS[index] = Histogram[index] - tmp;
    }
  }
}

template <class TEnsemble>
void calc_logDOS_full(const DiscreteAxis2D& Histogram,
                      const HistInfo2D& HistInfo,
                      const vec1<HistInfo2D>& HistInfos,
                      const vec1<std::pair<double, double>>& Parameters,
                      const double devmax, const int min_iteration,
                      vec1<double>& lnZ, DiscreteAxis2D& logDOS) {
  assert(HistInfos.size() == Parameters.size());
  assert(Parameters.size() == lnZ.size());
  assert(Histogram.size() == logDOS.size());
  assert(Histogram.get_min_first() == logDOS.get_min_first());
  assert(Histogram.get_max_second() == logDOS.get_max_second());
  double dev;
  int count = 1;
  for (int i = 0; i < min_iteration; ++i){
    iterate_logDOS<TEnsemble>(Histogram, HistInfo, HistInfos,  Parameters, lnZ, logDOS);
    auto new_lnZ = calc_lnZ<TEnsemble>(logDOS, HistInfo, Parameters);
    dev = deviation(lnZ, new_lnZ);
    std::swap(lnZ, new_lnZ);
    std::cout << count << ": " << dev << " " << devmax << " - full, forced\n";
    ++count;
  }
  while(dev > devmax){
    iterate_logDOS<TEnsemble>(Histogram, HistInfo, HistInfos, Parameters, lnZ,
                              logDOS);
    auto new_lnZ = calc_lnZ<TEnsemble>(logDOS, HistInfo, Parameters);
    dev = deviation(lnZ, new_lnZ);
    lnZ.swap(new_lnZ);
    std::cout << count << ": "<< dev << " " << devmax << " - full\n";
    ++count;
  };

  normalize_logDOS(logDOS, lnZ);
}

template <class TEnsemble>
void calc_logDOS_reduced(const DiscreteAxis2D& Histogram,
                         const HistInfo2D& HistInfo,
                         const vec1<HistInfo2D>& HistInfos,
                         const vec1<std::pair<double, double>>& Parameters,
                         const double devmax, vec1<double>& lnZ,
                         DiscreteAxis2D& logDOS) {
  double dev;
  int count = 1;
  do{
    iterate_logDOS<TEnsemble>(Histogram, HistInfo, HistInfos, Parameters, lnZ,
                              logDOS);
    auto new_lnZ = calc_lnZ_reduced<TEnsemble>(logDOS, HistInfos, Parameters);
    dev = deviation(lnZ, new_lnZ);
    lnZ.swap(new_lnZ);
    std::cout << count << ": " << dev << " " << devmax << " - reduced\n";
    ++count;
  } while(dev > devmax);
  normalize_logDOS(logDOS, lnZ);
}

template <class TEnsemble>
void calc_logDOS_reduced(const DiscreteAxis2D& Histogram,
                         const HistInfo2D& HistInfo,
                         const vec1<HistInfo2D>& HistInfos,
                         const vec1<std::pair<double, double>>& Parameters,
                         const double devmax, const vec1<int>& Overlap,
                         int treshold, vec1<double>& lnZ,
                         DiscreteAxis2D& logDOS) {
  double dev;
  int count = 1;
  do{
    iterate_logDOS<TEnsemble>(Histogram, HistInfo, HistInfos, Parameters, lnZ,
                              logDOS);
    auto new_lnZ = calc_lnZ_reduced<TEnsemble>(logDOS, HistInfos, Parameters,
                                               Overlap, treshold, lnZ);
    dev = deviation(lnZ, new_lnZ);
    std::cout << count << ": " << dev << " " << devmax << " - reduced, overlap\n";
    //if (not std::isfinite(dev)){
    //  std::cout << "lnZ " << lnZ << "\n";
    //  std::cout << "new lnZ " << new_lnZ << "\n";
    //}
    lnZ.swap(new_lnZ);
    ++count;
  } while(dev > devmax);
  normalize_logDOS(logDOS, lnZ);
}

/// call with first timeseries
///
/// returns in this order:
void logDOS_iteration_start(const vec2<double>& Timeseries, const int col1,
                            const int col2, const Range& range1,
                            const Range& range2, DiscreteAxis2D& Hist,
                            HistInfo2D& HistInfo,
                            vec1<HistInfo2D>& HistInfos, vec1<double>& lnZ) {
  assert(lnZ.size() == 1);
  std::tie(Hist, HistInfo) =
      make_histogram2d(Timeseries, col1, col2, range1, range2);
  assert(HistInfos.size() == 1);
  HistInfos.push_back(HistInfo);
}

template <class TEnsemble>
void logDOS_iteration_second(const vec2<double>& Timeseries, const int col1,
                             const int col2,
                             const vec1<std::pair<double, double>>& Parameters,
                             const double devmax, DiscreteAxis2D& Hist,
                             HistInfo2D& HistInfo,
                             vec1<HistInfo2D>& HistInfos, vec1<double>& lnZ,
                             DiscreteAxis2D& logDOS) {
  assert(Parameters.size() == 2);
  assert(lnZ.size() == 1);
  assert(HistInfos.size() == 1);
  lnZ.push_back(lnZ[0] + estimate_lnZ_from_hist<TEnsemble>(Hist, HistInfo, Parameters[0],
                                                Parameters[1]));
  add_histogram2d(Timeseries, col1, col2, Hist, HistInfo, HistInfos);
  assert(lnZ.size() == 2);
  assert(HistInfos.size() == 2);
  calc_logDOS_reduced<TEnsemble>(Hist, HistInfo, HistInfos, Parameters, devmax,
                                 lnZ, logDOS);
}

template <class TEnsemble>
void logDOS_iteration_next(const vec2<double>& Timeseries, const int col1,
                           const int col2, const double devmax,
                           const vec1<std::pair<double, double>>& Parameters,
                           DiscreteAxis2D& Hist,
                           HistInfo2D& HistInfo, vec1<HistInfo2D>& HistInfos,
                           vec1<double>& lnZ, DiscreteAxis2D& logDOS) {
  assert(Parameters.size() == lnZ.size() + 1);
  assert(lnZ.size() > 1);
  assert(HistInfos.size() > 1);
  lnZ.push_back(calc_lnZ<TEnsemble>(logDOS, HistInfo, Parameters.back()));
  add_histogram2d(Timeseries, col1, col2, Hist, HistInfo, HistInfos);
  assert(lnZ.size() > 2);
  assert(HistInfos.size() > 2);
  calc_logDOS_reduced<TEnsemble>(Hist, HistInfo, HistInfos, Parameters, devmax,
                                 lnZ, logDOS);
}

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

