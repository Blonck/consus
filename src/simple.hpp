#pragma once

#include <algorithm>
#include <tuple>

namespace consus
{

template<class T>
double mean(const T& Data){
  return std::accumulate(Data.begin(), Data.end(), 0.0);
}
 
/*-----------------------------------------------------------------------------
 *  min, max, mean value, variance, naive error, binning error
 *-----------------------------------------------------------------------------*/
template <class T>
std::tuple<double, double, double, double, double, double> simple_stats(
    const T& Data, const unsigned int NumBins) {
  assert(NumBins <= Data.size());
  size_t BinWidth = Data.size() / NumBins;
  double min = std::numeric_limits<double>::max();
  double max = -1.0 * min;
  double mean = 0.0, variance = 0.0, naiveerror = 0.0, binerr = 0.0;
  size_t j, i;
  double binmean;
  double num = 0;
  {
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
