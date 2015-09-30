#include <string>
#include <iostream>
#include <chrono>
#include <regex>
#include "../../src/load_file.hpp"
#include "../../src/all_matching.hpp"
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
  if (argc < 3){
    std::cerr << "arguments are missing\n";
    std::exit(1);
  }

  constexpr int colE1 = 0;
  constexpr int colE2 = 1;
  const double devmax = std::stod(argv[1]);
  const int NumBins = std::stoi(argv[2]);


  Range rangeE1;
  Range rangeE2;
  std::vector<std::string> filenames;
  std::vector<std::string> header;
  std::vector<std::pair<double, double>> Parameters;
  std::vector<double> initial_lnZ;
  DiscreteAxis2D initial_logDOS(rangeE1, rangeE2, log_zero<double>());
  bool add_eig;

  std::string path("analysis");
  std::string path_jk("analysis/JK");
  dlib::serialize(path + "/rangeE1.obj") << rangeE1;
  dlib::serialize(path + "/rangeE2.obj") << rangeE2;
  dlib::serialize(path + "/filenames.obj") << filenames;
  dlib::serialize(path + "/lnZ.obj") << initial_lnZ;
  dlib::serialize(path + "/logDOS.obj") << initial_logDOS;
  dlib::serialize(path + "/header.obj") << header;
  dlib::serialize(path + "/Parameters.obj") << Parameters;
  dlib::deserialize(path + "/add_eigenvalues.obj") >> add_eig;

  for (int j = 0; j < NumBins; ++j) {
    DiscreteAxis2D Hist;
    HistInfo2D HistInfo;
    std::vector<HistInfo2D> HistInfos;
    vec1<double> length_Hists;
    std::vector<DiscreteAxis2D> MicroMeans;
    std::vector<double> lnZ(initial_lnZ);
    DiscreteAxis2D logDOS(initial_logDOS);
    for (size_t i = 0; i < filenames.size(); ++i) {
      std::cout << filenames[i] << "\n";
      std::cout << "load file\n";
      auto timeseries = read_ssv_jk(filenames[i], j, NumBins);
      if (add_eig) {
        consus::eig::add_eigenvalues(timeseries, 5, 6, 7, 8, 9, 10);
      }
      if (i == 0) {
        std::tie(Hist, HistInfos, MicroMeans) =
            make_histogram2d(timeseries, colE1, colE2, rangeE1, rangeE2);
        HistInfo = HistInfos[0];
      } else {
        add_histogram2d(timeseries, colE1, colE2, Hist, HistInfo, HistInfos,
                        MicroMeans);
      }
    }

    normalize_MicroMeans(Hist, MicroMeans);
    std::cout << "last full WHAM\n";
    calc_logDOS_full<NVT>(Hist, HistInfo, HistInfos, Parameters, devmax, 10,
                          lnZ, logDOS);
    std::string path("analysis/JK/" + std::to_string(j));
    boost::filesystem::create_directories(path);
    dlib::serialize(path + "/Hist.obj") << Hist;
    dlib::serialize(path + "/lnZ.obj") << lnZ;
    dlib::serialize(path + "/logDOS.obj") << logDOS;
    dlib::serialize(path + "/MicroMeans.obj") << MicroMeans;
    dlib::serialize(path + "/HistInfo.obj") << HistInfo;
    dlib::serialize(path + "/HistInfos.obj") << HistInfos;
  }
}
