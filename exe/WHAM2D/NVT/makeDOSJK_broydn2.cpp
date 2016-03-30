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
  dlib::deserialize(path + "/Parameters.obj") >> Parameters;

  vec1<DiscreteAxis2D> jk_Hist(NumBins);
  vec1<HistInfo2D> jk_HistInfo(NumBins);
  vec2<HistInfo2D> jk_HistInfos(NumBins, vec1<HistInfo2D>(filenames.size()));
  vec1<vec1<double>> jk_length_Hists(NumBins);
  vec1<std::vector<double>> jk_lnZ(NumBins, initial_lnZ);
  vec1<DiscreteAxis2D> jk_logDOS(NumBins, initial_logDOS);

  boost::progress_display progress(filenames.size());
  for (size_t i = 0; i < filenames.size(); ++i) {
    ++progress;
    auto timeseries = read_ssv(filenames[i]);
    #pragma omp parallel for
    for (int j = 0; j < NumBins; ++j) {
      DiscreteAxis2D Hist;
      HistInfo2D HistInfo;
      std::tie(Hist, HistInfo) = make_histogram2d_jk(
          timeseries, colE1, colE2, rangeE1, rangeE2, j, NumBins);
      jk_HistInfos[j][i] = HistInfo;
      if (i == 0) {
        jk_Hist[j] = Hist;
        jk_HistInfo[j] = HistInfo;
      } else {
        merge_histograms2d(jk_Hist[j], jk_HistInfo[j], Hist, HistInfo );
      }
    }
  }

  boost::progress_display progress2(NumBins);
  #pragma omp parallel for
  for (int j = 0; j < NumBins; ++j) {
    ++progress2;
    calc_logDOS_full_broydn2<NVT>(jk_Hist[j], jk_HistInfo[j], jk_HistInfos[j],
                                  Parameters, devmax, jk_lnZ[j],
                                  jk_logDOS[j]);
    std::string path("analysis/JK/" + std::to_string(j));
    boost::filesystem::create_directories(path);
    dlib::serialize(path + "/Hist.obj") << jk_Hist[j];
    dlib::serialize(path + "/lnZ.obj") << jk_lnZ[j];
    dlib::serialize(path + "/logDOS.obj") << jk_logDOS[j];
    dlib::serialize(path + "/HistInfo.obj") << jk_HistInfo[j];
    dlib::serialize(path + "/HistInfos.obj") << jk_HistInfos[j];
  }
  dlib::serialize(path + "/NumBins.obj") << NumBins;
}
