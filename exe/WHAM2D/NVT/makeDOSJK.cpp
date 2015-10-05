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

using namespace consus;
using namespace consus::WHAM2D;

typedef std::chrono::high_resolution_clock myclock;

int main(int argc, char const *argv[])
{
  std::cout << "jk calculation\n";
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
  DiscreteAxis2D initial_logDOS;
  std::cout << "read in old state\n";

  std::string path("analysis");
  std::string path_jk("analysis/JK");
  dlib::deserialize(path + "/rangeE1.obj") >> rangeE1;
  dlib::deserialize(path + "/rangeE2.obj") >> rangeE2;
  dlib::deserialize(path + "/filenames.obj") >> filenames;
  dlib::deserialize(path + "/lnZ.obj") >> initial_lnZ;
  dlib::deserialize(path + "/logDOS.obj") >> initial_logDOS;
  dlib::deserialize(path + "/header.obj") >> header;
  dlib::deserialize(path + "/Parameters.obj") >> Parameters;

  std::cout << "start jk calculation\n";
  for (int j = 0; j < NumBins; ++j) {
    std::cout << j << "th jackknife bin\n";
    DiscreteAxis2D Hist;
    HistInfo2D HistInfo;
    std::vector<HistInfo2D> HistInfos;
    vec1<double> length_Hists;
    std::vector<double> lnZ(initial_lnZ);
    DiscreteAxis2D logDOS(initial_logDOS);
    for (size_t i = 0; i < filenames.size(); ++i) {
      std::cout << "load file\n";
      std::cout << filenames[i] << "\n";
      auto timeseries = read_ssv_jk(filenames[i], j, NumBins);
      if (i == 0) {
        std::tie(Hist, HistInfo) =
            make_histogram2d(timeseries, colE1, colE2, rangeE1, rangeE2);
        HistInfos.push_back(HistInfo);
      } else {
        add_histogram2d(timeseries, colE1, colE2, Hist, HistInfo, HistInfos);
      }
    }

    std::cout << "WHAM\n";
    calc_logDOS_full<NVT>(Hist, HistInfo, HistInfos, Parameters, devmax, 10,
                          lnZ, logDOS);
    std::string path("analysis/JK/" + std::to_string(j));
    boost::filesystem::create_directories(path);
    dlib::serialize(path + "/Hist.obj") << Hist;
    dlib::serialize(path + "/lnZ.obj") << lnZ;
    dlib::serialize(path + "/logDOS.obj") << logDOS;
    dlib::serialize(path + "/HistInfo.obj") << HistInfo;
    dlib::serialize(path + "/HistInfos.obj") << HistInfos;
  }
  dlib::serialize(path + "/NumBins.obj") << NumBins;
}
