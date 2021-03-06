#pragma once

#include <vector>
#include <tuple>
#include <algorithm>
#include "HistInfo.hpp"
#include "../Range.hpp"
#include "../DiscreteAxis.hpp"
#include "../vector/helper.hpp"
#include "addlogwise.hpp"

namespace consus {
namespace WHAM2D {

///1d WHAM from 2 dimensional data
namespace OneD {
  
template <class TEnsemble>
std::tuple<DiscreteAxis, HistInfo>
make_histogram1d(const vec2<double>& Timeseries, int col1, int col2,
                 const Range& range, std::pair<double, double>& Parameter) {
  DiscreteAxis Histogram(range, 0.0);
  vec1<HistInfo> HistInfos;
  int ffi = std::numeric_limits<int>::max();
  int lfi = std::numeric_limits<int>::lowest();
  //TODO: test reduction before use it
  //#pragma omp parallel for reduction(+: Histogram)
  //                         reduction(min: ffi_first, ffi_second)
  //                         reduction(max: lfi_first, lfi_second)
  for (size_t j = 0; j < Timeseries[col1].size(); ++j) {
    const double E1 = Timeseries[col1][j];
    const double E2 = Timeseries[col2][j];
    const double E = TEnsemble::ham(E1, Parameter.second, E2);
    ffi = std::min(ffi, Histogram.get_bin(E));
    lfi = std::max(lfi, Histogram.get_bin(E));
    int index = Histogram.get_bin(E);
    Histogram[index] += 1.0;
  }
  //#pragma omp parallel for
  //for (int i = 0; i < Histogram.size(); ++i){
  //  if (Histogram[i] != 0.0){
  //    Histogram[i] = std::log(Histogram[i]);
  //  }else{
  //    Histogram[i] = log_zero<double>();
  //  }
  //}
  double log_length = std::log(Timeseries[col1].size());
  HistInfo HistInfo(log_length, ffi, lfi);
  return std::make_tuple(Histogram, HistInfo);
}

template <class TEnsemble>
std::tuple<DiscreteAxis, HistInfo> make_histogram2d_jk(
    const vec2<double>& Timeseries, int col1, int col2, const Range& range,
    std::pair<double, double>& Parameter,
    const size_t num_bin, const size_t num_bins) {
  size_t length_bins = Timeseries[0].size() / num_bins;
  DiscreteAxis Histogram(range, 0.0);
  vec1<HistInfo> HistInfos;
  int ffi = std::numeric_limits<int>::max();
  int lfi = std::numeric_limits<int>::lowest();
  //TODO: test reduction before use it
  //#pragma omp parallel for reduction(+: Histogram)
  //                         reduction(min: ffi_first, ffi_second)
  //                         reduction(max: lfi_first, lfi_second)
  for (size_t i = 0; i < num_bin; ++i) {
    for (size_t j = 0; j < length_bins; ++j) {
      size_t i_ts = i * length_bins + j;
      const double E1 = Timeseries[col1][i_ts];
      const double E2 = Timeseries[col2][i_ts];
      const double E = TEnsemble::ham(E1, Parameter.second, E2);
      ffi = std::min(ffi, Histogram.get_bin(E));
      lfi = std::max(lfi, Histogram.get_bin(E));
      int index = Histogram.get_bin(E);
      Histogram[index] += 1.0;
    }
  }
  for (size_t i = num_bin+1; i < num_bins; ++i) {
    for (size_t j = 0; j < length_bins; ++j) {
      size_t i_ts = i * length_bins + j;
      const double E1 = Timeseries[col1][i_ts];
      const double E2 = Timeseries[col2][i_ts];
      const double E = TEnsemble::ham(E1, Parameter.second, E2);
      ffi = std::min(ffi, Histogram.get_bin(E));
      lfi = std::max(lfi, Histogram.get_bin(E));
      int index = Histogram.get_bin(E);
      Histogram[index] += 1.0;
    }
  }
  //#pragma omp parallel for
  //for (int i = 0; i < Histogram.size(); ++i){
  //  if (Histogram[i] != 0.0){
  //    Histogram[i] = std::log(Histogram[i]);
  //  }else{
  //    Histogram[i] = log_zero<double>();
  //  }
  //}
  double log_length = std::log(Timeseries[col1].size()-length_bins);
  HistInfo HistInfo(log_length, ffi, lfi);
  return std::make_tuple(Histogram, HistInfo);
}

/// merge two 2d histogram and its corresponding HistInfo2D object
///
/// \param fullHist will contain the data of the second histogram plus its own
/// \param fullHistInfo will contain the merged information of both histograms
/// \param tmpHist will be merged into fullHist
/// \param tmpHistInfo will be merged into fullHistInfo
void merge_histograms2d(DiscreteAxis& fullHist, HistInfo& fullHistInfo,
                        const DiscreteAxis& tmpHist,
                        const HistInfo& tmpHistInfo) {
  assert(fullHist.size() == tmpHist.size());
  fullHistInfo.ffi = std::min(fullHistInfo.ffi, tmpHistInfo.ffi);
  fullHistInfo.lfi = std::max(fullHistInfo.lfi, tmpHistInfo.lfi);
  addlogwise(fullHistInfo.log_length, tmpHistInfo.log_length);

  //#pragma omp parallel for
  for (int i = 0; i < tmpHist.size(); ++i){
    //if (tmpHist[i] != log_zero<double>()){
    //  fullHist[i] = addlogwise(fullHist[i], tmpHist[i]);
    //}
    fullHist[i] += tmpHist[i];
  }
}

template <class TEnsemble>
double estimate_lnZ_from_hist(const DiscreteAxis& Hist,
                         const HistInfo& HistInfo,
                         double Parameter,
                         double Parameter_Target) {
  auto log_Z_ratio = log_zero<double>();
  // Z_t =  \sum_{E_1, E_2} H_i W_t(E_1, E_2)/W_i(E_1, E_2)
  //TODO: maybe pragma with log_Z_ratio shared
  for (int i = HistInfo.ffi; i <= HistInfo.lfi; ++i) {
    auto E = Hist.get_value(i);
      if (Hist[i] != log_zero<double>()) {
        log_Z_ratio = addlogwise(log_Z_ratio,
            Hist[i] + TEnsemble::log_weight(Parameter_Target, E) -
                TEnsemble::log_weight(Parameter, E));
      }
    }
  }
  return log_Z_ratio;
}

/// normalize lnZ
///
/// \param lnZ contains log(Z) for each parameter
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
double calc_lnZ(const DiscreteAxis2D& DOS,
                const std::pair<double, double>& Parameter) {
  double lnZ = log_zero<double>();
    for (int i = 0; i < DOS.get_num_bins_first(); ++i) {
      auto E1 = DOS.get_value_first(i);
      for (int j = 0; j < DOS.get_num_bins_second(); ++j) {
        auto E2 = DOS.get_value_second(j);
        int index = DOS.get_index(i, j);
        lnZ = addlogwise(
            lnZ, DOS[index] + TEnsemble::log_weight(Parameter.first, E1,
                                                    Parameter.second, E2));
      }
    }
  return lnZ;
}

template <class TEnsemble>
double calc_lnZ_reduced(const DiscreteAxis2D& DOS, const HistInfo2D& HistInfo,
                        const std::pair<double, double>& Parameter) {
  double lnZ = log_zero<double>();
    for (int i = HistInfo.ffi_first; i <= HistInfo.lfi_first; ++i) {
      auto E1 = DOS.get_value_first(i);
      for (int j = HistInfo.ffi_second; j <= HistInfo.lfi_second; ++j) {
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
vec1<double> calc_lnZ(const DiscreteAxis2D& DOS,
                      const vec1<std::pair<double, double>>& Parameters) {
  vec1<double> lnZ(Parameters.size(), log_zero<double>());
  #pragma omp parallel for
  for (size_t k = 0; k < Parameters.size(); ++k) {
    for (int i = 0; i < DOS.get_num_bins_first(); ++i) {
      auto E1 = DOS.get_value_first(i);
      for (int j = 0; j < DOS.get_num_bins_second(); ++j) {
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
    for (int i = HistInfos[k].ffi_first; i <= HistInfos[k].lfi_first; ++i) {
      auto E1 = DOS.get_value_first(i);
      for (int j = HistInfos[k].ffi_second; j <= HistInfos[k].lfi_second; ++j) {
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
      for (int i = HistInfos[k].ffi_first; i <= HistInfos[k].lfi_first; ++i) {
        auto E1 = DOS.get_value_first(i);
        for (int j = HistInfos[k].ffi_second; j <= HistInfos[k].lfi_second;
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
  for (int i = HistInfos[k].ffi_first; i <= HistInfos[k].lfi_first; ++i) {
    auto E1 = DOS.get_value_first(i);
    for (int j = HistInfos[k].ffi_second; j <= HistInfos[k].lfi_second; ++j) {
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
  for (auto i = HistInfo.ffi_first; i <= HistInfo.lfi_first; ++i) {
    auto E1 = logDOS.get_value_first(i);
    for (auto j = HistInfo.ffi_second; j <= HistInfo.lfi_second; ++j) {
      auto E2 = logDOS.get_value_second(j);
      auto index = logDOS.get_index(i, j);
      auto tmp = HistInfos[0].log_length +
                 TEnsemble::log_weight(Parameters[0].first, E1,
                                       Parameters[0].second, E2) -
                 lnZ[0];
      for (size_t k = 1; k < lnZ.size(); ++k) {
        tmp = addlogwise(tmp,
                         HistInfos[k].log_length +
                             TEnsemble::log_weight(Parameters[k].first, E1,
                                                   Parameters[k].second, E2) -
                             lnZ[k]);
      }
      if (Histogram[index] != 0){
        logDOS[index] = std::log(Histogram[index]) - tmp;
      }else{
        logDOS[index] = log_zero<double>() - tmp;
      }
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
    iterate_logDOS<TEnsemble>(Histogram, HistInfo, HistInfos, Parameters, lnZ,
                              logDOS);
    auto new_lnZ = calc_lnZ<TEnsemble>(logDOS, Parameters);
    dev = deviation(lnZ, new_lnZ);
    std::swap(lnZ, new_lnZ);
    std::cout << count << ": " << dev << " " << devmax << " - full, forced\n";
    ++count;
  }
  while(dev > devmax){
    iterate_logDOS<TEnsemble>(Histogram, HistInfo, HistInfos,
                                        Parameters, lnZ, logDOS);
    auto new_lnZ = calc_lnZ<TEnsemble>(logDOS, Parameters);
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
    lnZ.swap(new_lnZ);
    ++count;
  } while(dev > devmax);
  normalize_logDOS(logDOS, lnZ);
}

} /* end of namespace OneD */ 

} /* end of namespace WHAM2D */ 

} /* end of namespace consus */ 
