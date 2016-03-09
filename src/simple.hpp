#pragma once

#include <algorithm>
#include <tuple>

/*-----------------------------------------------------------------------------
 *  this file contains some standard statistical functions
 *-----------------------------------------------------------------------------*/
namespace consus
{

/// calculates the mean of Data
///
/// just a convinient call to std::accumulate
/// \tparam T   - container with begin() and end()
/// \param Data - container with something convertable to a floating point number
template<class T>
double mean(const T& Data){
  return std::accumulate(Data.begin(), Data.end(), 0.0);
}
 
/// calculates the minimum/maximum/mean value/variance/naive error/binning errr
/// of Data
///
/// ignores numbers which are not finite (NaN/inf/etc)
///
/// \tparam Data - container with random access and size() method
/// \param Data - container with something convertable to a floating point number
/// \param NumBins - number of bins used to calculate binning error
/// \result - tuple with 6 elements (min, max, mean, variance, naive error, binning error)
template <class T>
std::tuple<double, double, double, double, double, double> simple_stats(
    const T& Data, const size_t NumBins) {
  assert(NumBins <= Data.size());
  size_t BinWidth = Data.size() / NumBins;
  double min = std::numeric_limits<double>::max();
  double max = -1.0 * min;
  double mean = 0.0, variance = 0.0, naiveerror = 0.0, binerr = 0.0;
  size_t j, i;
  double binmean;
  double num = 0;
  for (i = 0; i < NumBins; ++i) {
    binmean = 0.0;
    unsigned nums = 0;
    for (j = 0; j < BinWidth; ++j) {
      if (std::isfinite(Data[j + i * BinWidth])) {
        ++nums;
        min = std::min(Data[j + i * BinWidth], min);
        max = std::max(Data[j + i * BinWidth], max);
        mean += Data[j + i * BinWidth];
        variance += Data[j + i * BinWidth] * Data[j + i * BinWidth];
        binmean += Data[j + i * BinWidth];
      }
    }
    num += nums;
    binmean /= nums;
    binerr += binmean * binmean;
  }
  mean /= num;
  variance -= num * mean * mean;
  variance /= (num + 1.0);
  naiveerror = std::sqrt(variance / num);
  binerr -= mean * mean * static_cast<double>(NumBins);
  binerr /= static_cast<double>(NumBins + 1.0);
  binerr = std::sqrt(binerr / (static_cast<double>(NumBins)));
  return std::make_tuple(min, max, mean, variance, naiveerror, binerr);
}

template <class T, class TFunc>
std::tuple<double, double, double, double, double, double> simple_stats_func(
    const T& Data, const size_t NumBins, TFunc func) {
  assert(NumBins <= Data.size());
  size_t BinWidth = Data.size() / NumBins;
  double min = std::numeric_limits<double>::max();
  double max = -1.0 * min;
  double mean = 0.0, variance = 0.0, naiveerror = 0.0, binerr = 0.0;
  size_t j, i;
  double binmean;
  double num = 0;
  for (i = 0; i < NumBins; ++i) {
    binmean = 0.0;
    unsigned nums = 0;
    for (j = 0; j < BinWidth; ++j) {
      auto data = func(Data[j+i*BinWidth]);
      if (std::isfinite(data)) {
        ++nums;
        min = std::min(data, min);
        max = std::max(data, max);
        mean += data;
        variance += data * data;
        binmean += data;
      }
    }
    num += nums;
    binmean /= nums;
    binerr += binmean * binmean;
  }
  mean /= num;
  variance -= num * mean * mean;
  variance /= (num + 1.0);
  naiveerror = std::sqrt(variance / num);
  binerr -= mean * mean * static_cast<double>(NumBins);
  binerr /= static_cast<double>(NumBins + 1.0);
  binerr = std::sqrt(binerr / (static_cast<double>(NumBins)));
  return std::make_tuple(min, max, mean, variance, naiveerror, binerr);
}

template <class T, class TFunc>
std::tuple<double, double, double, double, double, double> simple_stats_func(
    const T& Data1, const T& Data2, const size_t NumBins, TFunc func) {
  assert(NumBins <= Data1.size());
  assert(Data1.size() == Data2.size());
  size_t BinWidth = Data1.size() / NumBins;
  double min = std::numeric_limits<double>::max();
  double max = -1.0 * min;
  double mean = 0.0, variance = 0.0, naiveerror = 0.0, binerr = 0.0;
  size_t j, i;
  double binmean;
  double num = 0;
  for (i = 0; i < NumBins; ++i) {
    binmean = 0.0;
    unsigned nums = 0;
    for (j = 0; j < BinWidth; ++j) {
      auto data = func(Data1[j+i*BinWidth], Data2[j+i*BinWidth]);
      if (std::isfinite(data)) {
        ++nums;
        min = std::min(data, min);
        max = std::max(data, max);
        mean += data;
        variance += data * data;
        binmean += data;
      }
    }
    num += nums;
    binmean /= nums;
    binerr += binmean * binmean;
  }
  mean /= num;
  variance -= num * mean * mean;
  variance /= (num + 1.0);
  naiveerror = std::sqrt(variance / num);
  binerr -= mean * mean * static_cast<double>(NumBins);
  binerr /= static_cast<double>(NumBins + 1.0);
  binerr = std::sqrt(binerr / (static_cast<double>(NumBins)));
  return std::make_tuple(min, max, mean, variance, naiveerror, binerr);
}

template <class T, class TFunc>
std::tuple<double, double, double, double, double, double> simple_stats_func(
    const T& Data1, const T& Data2, const T& Data3, const size_t NumBins, TFunc func) {
  assert(NumBins <= Data1.size());
  assert(Data1.size() == Data2.size());
  size_t BinWidth = Data1.size() / NumBins;
  double min = std::numeric_limits<double>::max();
  double max = -1.0 * min;
  double mean = 0.0, variance = 0.0, naiveerror = 0.0, binerr = 0.0;
  size_t j, i;
  double binmean;
  double num = 0;
  for (i = 0; i < NumBins; ++i) {
    binmean = 0.0;
    unsigned nums = 0;
    for (j = 0; j < BinWidth; ++j) {
      auto data = func(Data1[j+i*BinWidth], Data2[j+i*BinWidth], Data3[j+i*BinWidth]);
      if (std::isfinite(data)) {
        ++nums;
        min = std::min(data, min);
        max = std::max(data, max);
        mean += data;
        variance += data * data;
        binmean += data;
      }
    }
    num += nums;
    binmean /= nums;
    binerr += binmean * binmean;
  }
  mean /= num;
  variance -= num * mean * mean;
  variance /= (num + 1.0);
  naiveerror = std::sqrt(variance / num);
  binerr -= mean * mean * static_cast<double>(NumBins);
  binerr /= static_cast<double>(NumBins + 1.0);
  binerr = std::sqrt(binerr / (static_cast<double>(NumBins)));
  return std::make_tuple(min, max, mean, variance, naiveerror, binerr);
}

} /* end of namespace consus */ 
