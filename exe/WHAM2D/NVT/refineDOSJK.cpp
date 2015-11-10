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
  std::cout << "jk calculation\n";
  if (argc < 2){
    std::cerr << "arguments are missing\n";
    std::exit(1);
  }

  const double devmax = std::stod(argv[1]);

  Range rangeE1;
  Range rangeE2;
  std::vector<std::string> filenames;
  std::vector<std::pair<double, double>> Parameters;
  int NumBins;

  std::string path("analysis");
  std::string path_jk("analysis/JK");
  dlib::deserialize(path + "/rangeE1.obj") >> rangeE1;
  dlib::deserialize(path + "/rangeE2.obj") >> rangeE2;
  dlib::deserialize(path + "/filenames.obj") >> filenames;
  dlib::deserialize(path + "/Parameters.obj") >> Parameters;
  dlib::deserialize(path + "/NumBins.obj") >> NumBins;

  boost::progress_display progress2(NumBins);
  #pragma omp parallel for
  for (int j = 0; j < NumBins; ++j) {
    std::string path("analysis/JK/" + std::to_string(j));
    DiscreteAxis2D Hist;
    HistInfo2D HistInfo;
    std::vector<HistInfo2D> HistInfos;
    vec1<double> length_Hists;
    std::vector<double> lnZ;
    DiscreteAxis2D logDOS;
    dlib::deserialize(path + "/lnZ.obj") >> lnZ;
    dlib::deserialize(path + "/logDOS.obj") >> logDOS;
    dlib::deserialize(path + "/Hist.obj") >> Hist;
    dlib::deserialize(path + "/HistInfo.obj") >> HistInfo;
    dlib::deserialize(path + "/HistInfos.obj") >> HistInfos;

    ++progress2;
    calc_logDOS_full<NVT>(Hist, HistInfo, HistInfos,
                          Parameters, devmax, 3, lnZ, logDOS);
    boost::filesystem::create_directories(path);
    dlib::serialize(path + "/lnZ.obj") << lnZ;
    dlib::serialize(path + "/logDOS.obj") << logDOS;
  }
}
