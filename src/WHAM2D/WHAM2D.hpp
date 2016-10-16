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
#include "../optimization.hpp"

#include "dlib/time_this.h"

namespace consus
{

/// everthing one need to calculate the 2 dimensional density of states
/// dependent on two variables
namespace WHAM2D
{

#pragma omp declare reduction( +: DiscreteAxis2D : \
    std::transform(omp_in.begin(), omp_in.end(),   \
      omp_out.begin(), omp_out.begin(),            \
      std::plus<double>()) )                       \
      initializer (omp_priv(omp_orig))

/// constructs a 2D histogram and the corresponding HistInfo2D object
///
/// \param Timeseries - vector of vectors (2dim vector)
///                     each vector contains a single time series 
/// \param col1 - column of the first energy term of the Hamiltonian E1
/// \param col2 - column of the second energy term of the Hamiltonian E2
/// \param range1 - range for discretizing E1
/// \param range2 - range for discretizing E2
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
  //for (int i = 0; i < Histogram.size(); ++i){
  //  if (Histogram[i] != 0.0){
  //    Histogram[i] = std::log(Histogram[i]);
  //  }else{
  //    Histogram[i] = log_zero<double>();
  //  }
  //}
  double log_length = std::log(Timeseries[col1].size());
  HistInfo2D HistInfo(log_length, ffi_first, lfi_first, ffi_second,
                          lfi_second);
  return std::make_tuple(Histogram, HistInfo);
}

/// constructs  2D histogram and the corresponding HistInfo2D object 
/// from a the jackknifed time series
///
/// \param Timeseries - vector of vectors (2dim vector)
///                     each vector contains a single time series 
/// \param col1 column of the first energy term of the Hamiltonian E1
/// \param col2 column of the second energy term of the Hamiltonian E2
/// \param range1 range for discretizing E1
/// \param range2 range for discretizing E2
/// \param num_bin number of the constructed jackknifed histogram
/// \param num_bins total number of used jackknife blocks
/// \return a tuple with two elements, the first is the histogram and the second the
///         HistInfo2D object of the histogram
std::tuple<DiscreteAxis2D, HistInfo2D> make_histogram2d_jk(
    const vec2<double>& Timeseries, int col1, int col2, const Range& range1,
    const Range& range2, const size_t num_bin, const size_t num_bins) {
  size_t length_bins = Timeseries[0].size() / num_bins;
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
  for (size_t i = 0; i < num_bin; ++i) {
    for (size_t j = 0; j < length_bins; ++j) {
      size_t i_ts = i * length_bins + j;
      const double E1 = Timeseries[col1][i_ts];
      const double E2 = Timeseries[col2][i_ts];
      ffi_first = std::min(ffi_first, Histogram.get_bin_first(E1));
      lfi_first = std::max(lfi_first, Histogram.get_bin_first(E1));
      ffi_second = std::min(ffi_second, Histogram.get_bin_second(E2));
      lfi_second = std::max(lfi_second, Histogram.get_bin_second(E2));
      int index = Histogram.get_bin(E1, E2);
      Histogram[index] += 1.0;
    }
  }
  for (size_t i = num_bin+1; i < num_bins; ++i) {
    for (size_t j = 0; j < length_bins; ++j) {
      size_t i_ts = i * length_bins + j;
      const double E1 = Timeseries[col1][i_ts];
      const double E2 = Timeseries[col2][i_ts];
      ffi_first = std::min(ffi_first, Histogram.get_bin_first(E1));
      lfi_first = std::max(lfi_first, Histogram.get_bin_first(E1));
      ffi_second = std::min(ffi_second, Histogram.get_bin_second(E2));
      lfi_second = std::max(lfi_second, Histogram.get_bin_second(E2));
      int index = Histogram.get_bin(E1, E2);
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
  HistInfo2D HistInfo(log_length, ffi_first, lfi_first, ffi_second,
                          lfi_second);
  return std::make_tuple(Histogram, HistInfo);
}

/// merge two 2d histogram and its corresponding HistInfo2D object
///
/// \param fullHist will contain the data of the second histogram plus its own
/// \param fullHistInfo will contain the merged information of both histograms
/// \param tmpHist will be merged into fullHist
/// \param tmpHistInfo will be merged into fullHistInfo
void merge_histograms2d(DiscreteAxis2D& fullHist, HistInfo2D& fullHistInfo,
                        const DiscreteAxis2D& tmpHist,
                        const HistInfo2D& tmpHistInfo) {
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

  //#pragma omp parallel for
  for (int i = 0; i < tmpHist.size(); ++i){
    //if (tmpHist[i] != log_zero<double>()){
    //  fullHist[i] = addlogwise(fullHist[i], tmpHist[i]);
    //}
    fullHist[i] += tmpHist[i];
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

  //#pragma omp parallel for
  for (int i = 0; i < tmpHist.size(); ++i){
    //if (tmpHist[i] != 0){
    //  tmpHist[i] = std::log(tmpHist[i]);
    //  Histogram[i] = addlogwise(Histogram[i],tmpHist[i]);
    //}
    Histogram[i] += tmpHist[i];
  }
  HistInfo.ffi_first = std::min(ffi_first, HistInfo.ffi_first);
  HistInfo.lfi_first = std::max(lfi_first, HistInfo.lfi_first);
  HistInfo.ffi_second = std::min(ffi_second, HistInfo.ffi_second);
  HistInfo.lfi_second = std::max(lfi_second, HistInfo.lfi_second);
  HistInfo.log_length =
      addlogwise(HistInfo.log_length, std::log(Timeseries[col1].size()));
  HistInfos.push_back({std::log(Timeseries[col1].size()), ffi_first, lfi_first,
                       ffi_second, lfi_second});
  return tmpHist;
}

/// preliminary lnZ for Parameter from histogram measured at Parameter_Target
///
/// the function implements the logarithmic version of the following sum
/// Z = \sum( H(E1, E2)*exp(Weight_{Target}(E1, E2))/exp(Weight_{measured}(E1, E2)))
///
/// \param Hist - 2d histogram measured at Parameter
/// \param HistInfo - info object for Hist
/// \param Parameter - parameters for Hist
/// \param Parameter_Target - lnZ is estimated for this parameters
template <class TEnsemble>
double estimate_lnZ_from_hist(const DiscreteAxis2D& Hist,
                         const HistInfo2D& HistInfo,
                         const std::pair<double, double>& Parameter,
                         const std::pair<double, double>& Parameter_Target) {
  auto log_Z_ratio = log_zero<double>();
  // Z_t =  \sum_{E_1, E_2} H_i W_t(E_1, E_2)/W_i(E_1, E_2)
  //TODO: maybe pragma with log_Z_ratio shared
  for (int i = HistInfo.ffi_first; i <= HistInfo.lfi_first; ++i) {
    auto E1 = Hist.get_value_first(i);
    for (int j = HistInfo.ffi_second; j <= HistInfo.lfi_second; ++j) {
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

/// calculate lnZ for a given logathmic density of states (DOS)
/// at a certain Parameter
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

/// same as calc_lnZ(const DiscreteAxis2D& DOS, const std::pair<double,
/// double>& Parameter) but only over the reduced set of histogram entries
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

// TODO: use HistInfos for each parameter instead of whole area from HistInfo
/// same as calc_lnZ(const DiscreteAxis2D& DOS, const std::pair<double,
/// double>& Parameter but for a whole set of parameters at once
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

///same as calc_lnZ_reduced(const DiscreteAxis2D& DOS, const HistInfo2D&
///HistInfo, const std::pair<double, double>& Parameter) { but for a whole set
///of parameters at once
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

/// same as calc_lnZ_reduced(const DiscreteAxis2D& DOS, const vec1<HistInfo2D>&
/// HistInfos, const vec1<std::pair<double, double>>& Parameters)
/// but incoperating the overlap of the histograms
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
        for (int j = HistInfos[k].ffi_second; j <= HistInfos[k].lfi_second; ++j) {
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

/// calculates the deviation between two vectors of lnZ
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

double maxdev(const vec1<double>& lnZ_old, const vec1<double>& lnZ_new){
  double max = std::abs(lnZ_old[0] - lnZ_new[0]);
  for (size_t i = 1; i < lnZ_new.size(); ++i){
    max = std::max(max, std::abs(lnZ_old[i] - lnZ_new[i]));
  }
  return max;
}   

/// one iteration step in the WHAM procedure
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

/// calculates the logathmic density of states (DOS) from 
/// a given set of histograms an lnZ
/// 
/// full means do not use any overlap or histogram boundaries
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
    std::cout << maxdev(lnZ, new_lnZ) << "\n";
    ++count;
  }
  while(dev > devmax){
    iterate_logDOS<TEnsemble>(Histogram, HistInfo, HistInfos,
                                        Parameters, lnZ, logDOS);
    auto new_lnZ = calc_lnZ<TEnsemble>(logDOS, Parameters);
    dev = deviation(lnZ, new_lnZ);
    lnZ.swap(new_lnZ);
    std::cout << count << ": "<< dev << " " << devmax << " - full\n";
    std::cout << maxdev(lnZ, new_lnZ) << "\n";
    ++count;
  };

  normalize_logDOS(logDOS, lnZ);
}

/// same as calc_logDOS_reduced but i
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
    std::cout << maxdev(lnZ, new_lnZ) << "\n";
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
    std::cout << maxdev(lnZ, new_lnZ) << "\n";
    lnZ.swap(new_lnZ);
    ++count;
  } while(dev > devmax);
  normalize_logDOS(logDOS, lnZ);
}

//template <class TEnsemble>
//void calc_logDOS_full_lbfgs(const DiscreteAxis2D& Histogram,
//                            const HistInfo2D& HistInfo,
//                            const vec1<HistInfo2D>& HistInfos,
//                            const vec1<std::pair<double, double>>& Parameters,
//                            const double devmax, const int min_iteration,
//                            vec1<double>& lnZ, DiscreteAxis2D& logDOS) {
//  column_vector starting_point(lnZ.size());
//  std::copy(lnZ.begin(), lnZ.end(), starting_point.begin());
//  calc_J_lnZ<TEnsemble> func = calc_J_lnZ<TEnsemble>{Histogram, HistInfo, HistInfos, Parameters};
//  dlib::find_min_using_approximate_derivatives(
//      dlib::lbfgs_search_strategy(10),
//      dlib::objective_delta_stop_strategy(devmax).be_verbose(), func,
//      starting_point, -1);
//  std::copy(starting_point.begin(), starting_point.end(), lnZ.begin());
//  std::cout << func(starting_point) << "\n";
//  iterate_logDOS<TEnsemble>(Histogram, HistInfo, HistInfos, Parameters, lnZ, logDOS);
//  normalize_logDOS(logDOS, lnZ);
//}

template <class TEnsemble>
void calc_logDOS_full_broydn2(const DiscreteAxis2D& Histogram,
                              const HistInfo2D& HistInfo,
                              const vec1<HistInfo2D>& HistInfos,
                              const vec1<std::pair<double, double>>& Parameters,
                              const double devmax, vec1<double>& lnZ,
                              DiscreteAxis2D& logDOS) {
  std::cout << "start broydn2\n";
  auto func = [&Histogram, &HistInfo, &HistInfos, &Parameters, &logDOS](vec1<double>& lnZ) -> vec1<double>{
    iterate_logDOS<TEnsemble>(Histogram, HistInfo, HistInfos, Parameters, lnZ,
                              logDOS);
    auto diff_lnZ = calc_lnZ<TEnsemble>(logDOS, Parameters);
    for (size_t i = 0; i < diff_lnZ.size(); ++i){
      diff_lnZ[i] -= lnZ[i];
    }
    return diff_lnZ;
  };
  broydn2(lnZ, func, 0.1, devmax, 12);
  iterate_logDOS<TEnsemble>(Histogram, HistInfo, HistInfos, Parameters, lnZ, logDOS);
  normalize_logDOS(logDOS, lnZ);
}

template <class TEnsemble>
void calc_logDOS_reduced_broydn2(const DiscreteAxis2D& Histogram,
                         const HistInfo2D& HistInfo,
                         const vec1<HistInfo2D>& HistInfos,
                         const vec1<std::pair<double, double>>& Parameters,
                         const double devmax, vec1<double>& lnZ,
                         DiscreteAxis2D& logDOS) {
  std::cout << "start broydn2\n";
  auto func = [&Histogram, &HistInfo, &HistInfos, &Parameters,
               &logDOS](vec1<double>& lnZ) -> vec1<double> {
    iterate_logDOS<TEnsemble>(Histogram, HistInfo, HistInfos, Parameters, lnZ,
                              logDOS);
    auto diff_lnZ = calc_lnZ_reduced<TEnsemble>(logDOS, HistInfos, Parameters);
    for (size_t i = 0; i < diff_lnZ.size(); ++i){
      diff_lnZ[i] -= lnZ[i];
    }
    return diff_lnZ;
  };
  broydn2(lnZ, func, 0.1, devmax, 12);
  iterate_logDOS<TEnsemble>(Histogram, HistInfo, HistInfos, Parameters, lnZ, logDOS);
  normalize_logDOS(logDOS, lnZ);
}

} /* end of namespace WHAM2D */ 

} /* end of namespace consus */ 
