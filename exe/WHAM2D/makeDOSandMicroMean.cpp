#include <string>
#include <iostream>
#include <chrono>
#include <regex>
#include "../../src/load_file.hpp"
#include "../../src/vector/helper.hpp"
#include "../../src/DiscreteAxis2D.hpp"
#include "../../src/WHAM2D/WHAM2D.hpp"
#include "../../src/eigenvalues/add_eigenvalues.hpp"

#include "boost/filesystem.hpp"

using namespace consus;
using namespace consus::WHAM2D;

typedef std::chrono::high_resolution_clock myclock;

int main(int argc, char const *argv[])
{
  if (argc < 5){
    std::cerr << "arguments are missing\n";
    std::exit(1);
  }

  constexpr int colE1 = 0;
  constexpr int colE2 = 1;
  const double stepE1 = std::stod(argv[1]);
  const double stepE2 = std::stod(argv[2]);
  const double devmax = std::stod(argv[3]);

  std::vector<std::string> filenames;
  for (int i = 4; i < argc; ++i){
    filenames.push_back(argv[i]);
    std::cout << argv[i] << "\n";
  }
  std::vector<std::string> header = read_header(filenames[0]);
  bool add_eig = false;
  if (header.size() > 10) {
    if (header[5] == "Rxx" and header[6] == "Ryy" and header[7] ==
        "Rzz" and header[8] == "Rxy" and header[9] == "Rxz" and header[10] ==
        "Ryz") {
      add_eig = true;
    }
  }
  if (add_eig){
    consus::eig::add_eigenvalues_header(header);
  }
  std::cout << header << "\n";

  double minE1 = std::numeric_limits<double>::max();
  double maxE1 = std::numeric_limits<double>::lowest();
  double minE2 = std::numeric_limits<double>::max();
  double maxE2 = std::numeric_limits<double>::lowest();

  #pragma omp parallel for reduction(min:minE1, minE2) reduction(max:maxE1, maxE2)
  for (size_t i = 0; i < filenames.size(); ++i) {
    std::cout << filenames[i] << "\n";
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
  std::cout << minE1 << " " << maxE1 << " " << minE2 << " " << maxE2 << "\n";

  std::vector<std::pair<double, double>> Parameters;

  std::regex param_regex(".*timeseries/(.*)/PT_T(.*).dat", std::regex::egrep);
  std::smatch match;

  DiscreteAxis2D<double> Hist;
  HistInfo2D<double> HistInfo;
  std::vector<HistInfo2D<double>> HistInfos;
  vec1<double> length_Hists;
  std::vector<DiscreteAxis2D<double>> MicroMeans;
  DiscreteAxis2D<double> DOS;
  Range<double> rangeE1 = {minE1, stepE1, maxE1};
  Range<double> rangeE2 = {minE2, stepE2, maxE2};
  std::vector<double> lnZ(1, 0);
  std::vector<std::vector<bool>> Overlaps;
  DiscreteAxis2D<double> logDOS(rangeE1, rangeE2, log_zero<double>());
  for (size_t i = 0; i < filenames.size(); ++i) {
    std::cout << filenames[i] << "\n";
    std::regex_search(filenames[i], match, param_regex);
    if (match.size() != 3) {
      std::cerr << "could not determine parameters"
                << "\n";
      std::exit(1);
    }
    Parameters.push_back(
        std::make_pair(1.0 / std::stod(match[2]), std::stod(match[1])));
    std::cout << "load file\n";
    auto timeseries = read_ssv(filenames[i]);
    if (add_eig){
      consus::eig::add_eigenvalues(timeseries, 5, 6, 7, 8, 9, 10);
    }
    if( i == 0){
      logDOS_iteration_start(timeseries, colE1, colE2, rangeE1, rangeE2, Hist,
                             MicroMeans, HistInfo, HistInfos, lnZ);
      continue;
    }
    if ( i == 1){
      logDOS_iteration_second<NVT<double>>(timeseries, colE1, colE2, Parameters, Hist,
                              MicroMeans, HistInfo, HistInfos, lnZ, logDOS, Overlaps,
                              devmax);
      continue;
    }
    logDOS_iteration_next<NVT<double>>(timeseries, colE1, colE2, Parameters, Hist,
                          MicroMeans, HistInfo, HistInfos, lnZ, logDOS,
                          Overlaps, devmax);
  }
  normalize_MicroMeans(Hist, MicroMeans);
  calc_logDOS_full<NVT<double>>(Hist, HistInfo, HistInfos, Parameters, devmax, lnZ, logDOS);

  std::string path("analysis");
  std::cout << "RangeE1 " << rangeE1 << "\n";
  std::cout << "RangeE2 " << rangeE2 << "\n";
  boost::filesystem::create_directories(path);

  dlib::serialize(path + "/rangeE1.obj") << rangeE1;
  dlib::serialize(path + "/rangeE2.obj") << rangeE2;
  dlib::serialize(path + "/filenames.obj") << filenames;
  dlib::serialize(path + "/Hist.obj") << Hist;
  dlib::serialize(path + "/lnZ.obj") << lnZ;
  dlib::serialize(path + "/logDOS.obj") << DOS;
  dlib::serialize(path + "/MicroMeans.obj") << MicroMeans;
  dlib::serialize(path + "/header.obj") << header;
  dlib::serialize(path + "/HistInfo.obj") << HistInfo;
  dlib::serialize(path + "/HistInfos.obj") << HistInfos;
  dlib::serialize(path + "/Parameters.obj") << Parameters;
}
