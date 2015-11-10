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
  if (argc < 2){
    std::cerr << "arguments are missing\n";
    std::exit(1);
  }

  const double devmax = std::stod(argv[1]);
  std::string path("analysis");

  std::vector<std::string> filenames;
  std::vector<std::pair<double, double>> Parameters;
  DiscreteAxis2D Hist;
  HistInfo2D HistInfo;
  std::vector<HistInfo2D> HistInfos;
  vec1<double> length_Hists;
  std::vector<double> lnZ;
  DiscreteAxis2D logDOS;
  Range rangeE1;
  Range rangeE2;
  dlib::deserialize(path + "/rangeE1.obj") >> rangeE1;
  dlib::deserialize(path + "/rangeE2.obj") >> rangeE2;
  dlib::deserialize(path + "/filenames.obj") >> filenames;
  dlib::deserialize(path + "/lnZ.obj") >> lnZ;
  dlib::deserialize(path + "/logDOS.obj") >> logDOS;
  dlib::deserialize(path + "/Parameters.obj") >> Parameters;
  dlib::deserialize(path + "/Hist.obj") >> Hist;
  dlib::deserialize(path + "/HistInfo.obj") >> HistInfo;
  dlib::deserialize(path + "/HistInfos.obj") >> HistInfos;

  std::cout << "refine DOS" << "\n";
  calc_logDOS_full<NVT>(Hist, HistInfo, HistInfos, Parameters, devmax, 10, lnZ,
                        logDOS);

  dlib::serialize(path + "/lnZ.obj") << lnZ;
  dlib::serialize(path + "/logDOS.obj") << logDOS;
}
