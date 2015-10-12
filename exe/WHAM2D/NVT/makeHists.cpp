#define VERBOSE 

#include <string>
#include <iostream>
#include <chrono>
#include <regex>
#include "../../src/data/load_file.hpp"
#include "../../src/data/all_matching.hpp"
#include "../../src/vector/helper.hpp"
#include "../../src/DiscreteAxis2D.hpp"
#include "../../src/WHAM2D/WHAM2D.hpp"
#include "../../src/eigenvalues/add_eigenvalues.hpp"

#include "boost/filesystem.hpp"
#include "boost/progress.hpp"

using namespace consus;
using namespace consus::WHAM2D;

typedef std::chrono::high_resolution_clock myclock;

int main(int argc, char const *argv[])
{
  if (argc < 4){
    std::cerr << "arguments are missing\n";
    std::exit(1);
  }

  constexpr int colE1 = 0;
  constexpr int colE2 = 1;
  const double stepE1 = std::stod(argv[1]);
  const double stepE2 = std::stod(argv[2]);

  std::vector<std::string> filenames;
  for (int i = 3; i < argc; ++i){
    filenames.push_back(argv[i]);
  }

  std::regex param_regex(".*timeseries/(.*)/PT_T(.*).dat", std::regex::egrep);
  std::smatch match;

  // sort first in kappa and then in T
  std::sort(filenames.begin(), filenames.end(),
            [&](const std::string& a, const std::string& b) {
              std::regex_search(a, match, param_regex);
              if (match.size() != 3) {
                std::cerr << "could not determine parameters\n";
                std::exit(1);
              }
              double a_t = std::stod(match[2]);
              double a_kappa = std::stod(match[1]);
              std::regex_search(b, match, param_regex);
              if (match.size() != 3) {
                std::cerr << "could not determine parameters\n";
                std::exit(1);
              }
              double b_t = std::stod(match[2]);
              double b_kappa = std::stod(match[1]);
              if (a_kappa == b_kappa) {
                return a_t < b_t;
              } else {
                return a_kappa < b_kappa;
              }
            });

  // reverse every second sequence of temperatures
  bool flip = false;
  std::regex_search(filenames[0], match, param_regex);
  double cur_kappa = std::stod(match[1]);
  size_t start_i = 0;
  for (size_t i = 1; i < filenames.size(); ++i){
    std::regex_search(filenames[i], match, param_regex);
    double kappa = std::stod(match[1]);
    if (cur_kappa != kappa){
      cur_kappa = kappa;
      if (flip){
        std::reverse(filenames.begin() + start_i, filenames.begin() + i);
        flip = false;
      }else{
        flip = true;
        start_i = i;
      }
    }
  }
  if (flip){
    std::reverse(filenames.begin() + start_i, filenames.end());
  }

  double minE1 = std::numeric_limits<double>::max();
  double maxE1 = std::numeric_limits<double>::lowest();
  double minE2 = std::numeric_limits<double>::max();
  double maxE2 = std::numeric_limits<double>::lowest();

  std::cout << "searchin min/max values in data\n";
  boost::progress_display progress(filenames.size());
  #pragma omp parallel for reduction(min:minE1, minE2) reduction(max:maxE1, maxE2)
  for (size_t i = 0; i < filenames.size(); ++i) {
    #pragma omp critical
    ++progress;
    auto timeseries = read_ssv(filenames[i]);
    auto resE1 =
        std::minmax_element(timeseries[colE1].begin(), timeseries[colE1].end());
    auto resE2 =
        std::minmax_element(timeseries[colE2].begin(), timeseries[colE2].end());
    minE1 = std::min(minE1, *(resE1.first));
    maxE1 = std::max(maxE1, *(resE1.second));
    minE2 = std::min(minE2, *(resE2.first));
    maxE2 = std::max(maxE2, *(resE2.second));
  }

  std::vector<std::pair<double, double>> Parameters;
  vec2<int> Overlap;
  vec1<double> length_Hists;
  Range rangeE1 = {minE1, stepE1, maxE1};
  Range rangeE2 = {minE2, stepE2, maxE2};
  
  std::string path("analysis");
  std::string pathh("analysis/Hists");
  boost::filesystem::create_directories(path);
  boost::filesystem::create_directories(pathh);

  for (size_t i = 0; i < filenames.size(); ++i) {
    std::regex_search(filenames[i], match, param_regex);
    Parameters.push_back(
        std::make_pair(1.0 / std::stod(match[2]), std::stod(match[1])));
  }

  
  std::cout << "loading files\n";
  boost::progress_display progress2(filenames.size());
  #pragma omp parallel for
  for (size_t i = 0; i < filenames.size(); ++i) {
      DiscreteAxis2D Hist;
      HistInfo2D HistInfo;
      #pragma omp critical
      ++progress2;
      auto timeseries = read_ssv(filenames[i]);
      std::tie(Hist, HistInfo) =
          make_histogram2d(timeseries, colE1, colE2, rangeE1, rangeE2);
      dlib::serialize(pathh + "/Hist_" + std::to_string(i) + ".obj") << Hist;
      dlib::serialize(pathh + "/HistInfo_" + std::to_string(i) + ".obj")
          << HistInfo;
  }

  assert(Parameters.size() == filenames.size());

  std::cout << "RangeE1 " << rangeE1 << "\n";
  std::cout << "RangeE2 " << rangeE2 << "\n";

  for (size_t i = 0; i < Parameters.size(); ++i){
    dlib::serialize(pathh + "/Parameter_" + std::to_string(i) + ".obj")
        << Parameters[i];
  }

  dlib::serialize(path + "/rangeE1.obj") << rangeE1;
  dlib::serialize(path + "/rangeE2.obj") << rangeE2;
  dlib::serialize(path + "/filenames.obj") << filenames;
  dlib::serialize(path + "/Parameters.obj") << Parameters;
}
