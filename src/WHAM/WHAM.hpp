#pragma once

#include "../addlogwise.hpp"
#include "../DiscreteAxis.hpp"
#include "../Range.hpp"
#include "../vector/helper.hpp"
#include "Ensemble.hpp"

namespace consus
{

namespace WHAM
{

/** Construct a historgram from a time series
 *
 * @param timeseries Timeseries from which histogram will be constructed.
 * @param range Range object defining min, max, and step width.
 * @return Pair with histogram and sum of all entries.
 * */
std::pair<DiscreteAxis, double> make_histogram(const vec1<double>& timeseries,
                                               const Range& range) {
  DiscreteAxis hist(range, 0.0);
  for (auto& energy : timeseries){
    const int index = hist.get_bin(energy);
    hist[index] += 1.0;
  }

  return std::make_pair(hist, timeseries.size());
}

/** Construct a 'jackknifed' histogram from a time series.
 *
 * Identical to make_histogram(), but neglects the @a nth jackknife bin.
 * @param timeseries Timeseries from which histogram will be constructed.
 * @param range Range object defining min, max, and step width.
 * @param num_bin Which jackknife block will be neglected.
 * @param num_bins Total number of jackknife bins. 
 * @return Pair with histogram and sum of all entries.
 * */
std::pair<DiscreteAxis, double> make_histogram_jk(
    const vec1<double>& timeseries, const Range& range,
    const size_t num_bin, const size_t num_bins) {
  size_t length_bins = timeseries.size() / num_bins;
  DiscreteAxis hist(range, 0.0);
  double length_hist = 0;
  for (size_t i = 0; i < num_bin; ++i) {
    for (size_t j = 0; j < length_bins; ++j) {
      const size_t i_ts = i * length_bins + j;
      const double energy = timeseries[i_ts];
      const int index = hist.get_bin(energy);
      hist[index] += 1.0;
      ++length_hist;
    }
  }

  for (size_t i = num_bin + 1; i < num_bins; ++i) {
    for (size_t j = 0; j < length_bins; ++j) {
      const size_t i_ts = i * length_bins + j;
      const double energy = timeseries[i_ts];
      const int index = hist.get_bin(energy);
      hist[index] += 1.0;
      ++length_hist;
    }
  }

  return std::make_pair(hist, length_hist);
}

/** Merge one histogram into another.
 *
 * @param fullHist Histogram in which @tmpHist will be merged to. */
void merge_histograms(DiscreteAxis& fullHist, double& full_hist_length,
                      const DiscreteAxis& tmpHist,
                      const double tmp_hist_length) {
  assert(fullHist.size() == tmpHist.size());
  assert(fullHist.get_min() == tmpHist.get_min());
  assert(fullHist.get_max() == tmpHist.get_max());
  assert(fullHist.get_step() == tmpHist.get_step());

  #pragma omp parallel for
  for (int i = 0; i < tmpHist.size(); ++i){
    fullHist[i] += tmpHist[i];
  }

  full_hist_length += tmp_hist_length;
}

/// preliminary lnZ for Parameter_Target from histogram
/// measured at Parameter_Target
///
/// Hist - 2d histogram measured at Parameter
/// Parameter - parameters for Hist
/// Parameter_Target - lnZ is estimated for this parameters
template <class TEnsemble>
double estimate_lnZ_from_hist(const DiscreteAxis& hist, const double Parameter,
                              const double Parameter_Target) {
  auto log_Z_ratio = log_zero<double>();
  // Z_t =  \sum_{E} H_i W_t(E)/W_i(E)

  for (int i = 0; i < hist.get_num_bins(); ++i) {
    auto energy = hist.get_value(i);
    if (hist[i] != log_zero<double>()) {
      log_Z_ratio =
          addlogwise(log_Z_ratio,
                     hist[i] + TEnsemble::log_weight(Parameter_Target, energy) -
                         TEnsemble::log_weight(Parameter, energy));
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

void normalize_logDOS(DiscreteAxis& logDOS, const vec1<double>& lnZ) {
  assert(lnZ.size() > 0);
  #pragma omp parallel for
  for (int i = 0; i < logDOS.size(); ++i) {
    logDOS[i] -= lnZ[0];
  }
}

template <class TEnsemble>
double calc_lnZ(const DiscreteAxis& DOS, const double Parameter) {
  double lnZ = log_zero<double>();
  for (int i = 0; i < DOS.get_num_bins(); ++i) {
    auto energy = DOS.get_value(i);
    lnZ = addlogwise(lnZ, DOS[i] + TEnsemble::log_weight(Parameter, energy));
  }
  return lnZ;
}

template <class TEnsemble>
vec1<double> calc_lnZ(const DiscreteAxis& DOS,
                      const vec1<double>& parameters) {
  vec1<double> lnZ(parameters.size(), log_zero<double>());
#pragma omp parallel for
  for (size_t k = 0; k < parameters.size(); ++k) {
    for (int i = 0; i < DOS.get_num_bins(); ++i) {
      auto energy = DOS.get_value(i);
      lnZ[k] = addlogwise(
          lnZ[k], DOS[i] + TEnsemble::log_weight(parameters[k], energy));
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
void iterate_logDOS(const DiscreteAxis& histogram, const double log_length_hist,
                    const vec1<double>& parameters,
                    const vec1<double>& lnZ, DiscreteAxis& logDOS) {
  #pragma omp parallel for
  for (auto i = 0; i < histogram.get_num_bins(); ++i) {
    const auto energy = logDOS.get_value(i);
    auto tmp =
        log_length_hist + TEnsemble::log_weight(parameters[0], energy) - lnZ[0];
    for (size_t k = 1; k < lnZ.size(); ++k) {
      tmp = addlogwise(tmp, log_length_hist +
                                TEnsemble::log_weight(parameters[k], energy) -
                                lnZ[k]);
    }
    if (histogram[i] != 0) {
      logDOS[i] = std::log(histogram[i]) - tmp;
    } else {
      logDOS[i] = log_zero<double>() - tmp;
    }
  }
}

template <class TEnsemble>
void calc_logDOS_full(const DiscreteAxis& Histogram, double log_length_hist,
                      const vec1<double>& Parameters, const double devmax,
                      const int min_iteration, vec1<double>& lnZ,
                      DiscreteAxis& logDOS) {
  assert(Parameters.size() == lnZ.size());
  assert(Histogram.size() == logDOS.size());
  assert(Histogram.get_min() == logDOS.get_min());
  assert(Histogram.get_max() == logDOS.get_max());

  double dev;
  int count = 1;
  for (int i = 0; i < min_iteration; ++i){
    iterate_logDOS<TEnsemble>(Histogram, log_length_hist, Parameters, lnZ,
                              logDOS);
    auto new_lnZ = calc_lnZ<TEnsemble>(logDOS, Parameters);
    dev = deviation(lnZ, new_lnZ);
    std::swap(lnZ, new_lnZ);
    ++count;
  }
  while(dev > devmax){
    iterate_logDOS<TEnsemble>(Histogram, log_length_hist, Parameters, lnZ,
                              logDOS);
    auto new_lnZ = calc_lnZ<TEnsemble>(logDOS, Parameters);
    dev = deviation(lnZ, new_lnZ);
    lnZ.swap(new_lnZ);
    std::cout << count << ": " << dev << " " << devmax << " - full\n";
    ++count;
  };

  normalize_logDOS(logDOS, lnZ);
}

} /* end of namespace WHAM */ 

}
