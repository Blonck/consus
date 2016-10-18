#include <string>
#include <iostream>
#include <chrono>
#include <regex>
#include "../../src/data/load_file.hpp"
#include "../../src/data/all_matching.hpp"
#include "../../src/vector/helper.hpp"
#include "../../src/DiscreteAxis.hpp"
#include "../../src/WHAM/WHAM.hpp"
#include "../../src/eigenvalues/add_eigenvalues.hpp"

#include "boost/filesystem.hpp"
#include "boost/progress.hpp"

using namespace consus;
using namespace consus::WHAM;

typedef std::chrono::high_resolution_clock myclock;

int main(int argc, char const *argv[])
{
  std::cout << "jk calculation\n";
  if (argc < 3){
    std::cerr << "arguments are missing\n";
    std::exit(1);
  }

  constexpr int col_energy = 2;
  const double devmax = std::stod(argv[1]);
  const int NumBins = std::stoi(argv[2]);

  Range range;
  std::vector<std::string> filenames;
  std::vector<double> Parameters;
  std::vector<double> initial_lnZ;
  DiscreteAxis initial_logDOS;
  std::cout << "read in old state\n";

  std::string path("analysis");
  std::string path_jk("analysis/JK");
  dlib::deserialize(path + "/range.obj") >> range;
  dlib::deserialize(path + "/filenames.obj") >> filenames;
  dlib::deserialize(path + "/lnZ.obj") >> initial_lnZ;
  dlib::deserialize(path + "/logDOS.obj") >> initial_logDOS;
  dlib::deserialize(path + "/Parameters.obj") >> Parameters;

  vec1<DiscreteAxis> jk_Hist(NumBins);
  vec1<double> jk_length_hist(NumBins);
  vec1<vec1<double>> jk_length_Hists(NumBins, vec1<double>(filenames.size()));
  vec1<std::vector<double>> jk_lnZ(NumBins, initial_lnZ);
  vec1<DiscreteAxis> jk_logDOS(NumBins, initial_logDOS);

  boost::progress_display progress(filenames.size());
  for (size_t i = 0; i < filenames.size(); ++i) {
    ++progress;
    auto timeseries = read_ssv(filenames[i]);
    for (int j = 0; j < NumBins; ++j) {
      DiscreteAxis Hist;
      double length_hist;
      std::tie(Hist, length_hist) =
          make_histogram_jk(timeseries[col_energy], range, j, NumBins);
      jk_length_Hists[j][i] = length_hist;
      if (i == 0) {
        jk_Hist[j] = Hist;
        jk_length_hist[j] = length_hist;
      } else {
        merge_histograms(jk_Hist[j], jk_length_hist[j], Hist, length_hist );
      }
    }
  }

  boost::progress_display progress2(NumBins);
  #pragma omp parallel for
  for (int j = 0; j < NumBins; ++j) {
    ++progress2;
    calc_logDOS_full<NVT>(jk_Hist[j], jk_length_hist[j],
                          Parameters, devmax, 3, jk_lnZ[j], jk_logDOS[j]);
    std::string path("analysis/JK/" + std::to_string(j));
    boost::filesystem::create_directories(path);
    dlib::serialize(path + "/Hist.obj") << jk_Hist[j];
    dlib::serialize(path + "/lnZ.obj") << jk_lnZ[j];
    dlib::serialize(path + "/logDOS.obj") << jk_logDOS[j];
    dlib::serialize(path + "/length_hist.obj") << jk_length_hist[j];
    dlib::serialize(path + "/length_hists.obj") << jk_length_Hists[j];
  }
  dlib::serialize(path + "/NumBins.obj") << NumBins;
}
