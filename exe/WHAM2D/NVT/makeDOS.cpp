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
  std::string pathh("analysis/Hists");

  std::vector<std::string> filenames;
  std::vector<std::pair<double, double>> Parameters;
  DiscreteAxis2D Hist;
  HistInfo2D HistInfo;
  std::vector<HistInfo2D> HistInfos;
  vec1<double> length_Hists;
  std::vector<double> lnZ;
  Range rangeE1;
  Range rangeE2;
  dlib::deserialize(path + "/rangeE1.obj") >> rangeE1;
  dlib::deserialize(path + "/rangeE2.obj") >> rangeE2;

  DiscreteAxis2D logDOS(rangeE1, rangeE2, 0.0);

  dlib::deserialize(path + "/filenames.obj") >> filenames;

  for (size_t i = 0; i < filenames.size(); ++i) {
    DiscreteAxis2D tmpHist;
    HistInfo2D tmpHistInfo;
    vec1<int> Overlap;
    std::pair<double, double> Parameter;
    dlib::deserialize(pathh + "/Hist_" + std::to_string(i) + ".obj") >> tmpHist;
    dlib::deserialize(pathh + "/HistInfo_" + std::to_string(i) + ".obj")
        >> tmpHistInfo;
    dlib::deserialize(pathh + "/Overlap_" + std::to_string(i) + ".obj")
        >> Overlap;
    dlib::deserialize(pathh + "/Parameter_" + std::to_string(i) + ".obj")
        >> Parameter;
    HistInfos.push_back(tmpHistInfo);
    Parameters.push_back(Parameter);
    std::cout << "==================================================\n";
    std::cout << "Temperature " << 1.0 / Parameter.first << " Kappa "
              << Parameter.second << "\n";
    std::cout << "==================================================\n";
    if (i == 0) {
      HistInfo = tmpHistInfo;
      Hist = tmpHist;
      lnZ.push_back(0);
      assert(lnZ.size() == 1);
      assert(HistInfos.size() == 1);
      continue;
    }
    if (i == 1) {
      merge_histograms2d(Hist, HistInfo, tmpHist, tmpHistInfo);
      lnZ.push_back(lnZ[0] + estimate_lnZ_from_hist<NVT>(
                                 Hist, HistInfo, Parameters[0], Parameters[1]));
      calc_logDOS_reduced<NVT>(Hist, HistInfo, HistInfos, Parameters, devmax,
                               Overlap, 1, lnZ, logDOS);
      continue;
    }
    assert(Parameters.size() == lnZ.size() + 1);
    assert(lnZ.size() > 1);
    assert(HistInfos.size() > 1);
    merge_histograms2d(Hist, HistInfo, tmpHist, tmpHistInfo);
    lnZ.push_back(calc_lnZ_reduced<NVT>(logDOS, tmpHistInfo, Parameters.back()));
    assert(lnZ.size() > 2);
    assert(HistInfos.size() > 2);
    calc_logDOS_reduced<NVT>(Hist, HistInfo, HistInfos, Parameters, devmax,
                             Overlap, 1, lnZ, logDOS);
    //if (std::accumulate(Overlap.begin(), Overlap.end(), 0) > 3){
    //  calc_logDOS_reduced<NVT>(Hist, HistInfo, HistInfos, Parameters, devmax, Overlap, 1,
    //                           lnZ, logDOS);
    //}else{
    //  std::cout << "small overlap: using full calculation\n";
    //  calc_logDOS_full<NVT>(Hist, HistInfo, HistInfos, Parameters, devmax, 3,
    //                        lnZ, logDOS);
    //}
  }

  std::cout << "last full WHAM" << "\n";
  calc_logDOS_full<NVT>(Hist, HistInfo, HistInfos, Parameters, devmax, 10, lnZ,
                        logDOS);

  dlib::serialize(path + "/Hist.obj") << Hist;
  dlib::serialize(path + "/lnZ.obj") << lnZ;
  dlib::serialize(path + "/logDOS.obj") << logDOS;
  dlib::serialize(path + "/HistInfo.obj") << HistInfo;
  dlib::serialize(path + "/HistInfos.obj") << HistInfos;
}
