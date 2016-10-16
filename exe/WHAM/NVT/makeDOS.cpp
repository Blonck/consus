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

using namespace consus;
using namespace consus::WHAM;

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
  std::vector<double> Betas;
  DiscreteAxis Hist;
  double length_hist;
  std::vector<double> length_hists;
  vec1<double> length_Hists;
  std::vector<double> lnZ;
  Range range;
  dlib::deserialize(path + "/range.obj") >> range;

  DiscreteAxis logDOS(range, 0.0);
  dlib::deserialize(path + "/filenames.obj") >> filenames;

  for (size_t i = 0; i < filenames.size(); ++i) {
    DiscreteAxis tmpHist;
    double tmp_length_hist;
    double Parameter;
    dlib::deserialize(pathh + "/Hist_" + std::to_string(i) + ".obj") >> tmpHist;
    dlib::deserialize(pathh + "/length_hist_" + std::to_string(i) + ".obj")
        >> tmp_length_hist;
    dlib::deserialize(pathh + "/Parameter_" + std::to_string(i) + ".obj")
        >> Parameter;

    length_hists.push_back(tmp_length_hist);
    Betas.push_back(Parameter);
    std::cout << "==================================================\n";
    std::cout << "Temperature " << 1.0 / Parameter << "\n";
    std::cout << "==================================================\n";
    if (i == 0) {
      length_hist = tmp_length_hist;
      Hist = tmpHist;
      lnZ.push_back(0);
      assert(lnZ.size() == 1);
      assert(length_hists.size() == 1);
      continue;
    }
    if (i == 1) {
      merge_histograms(Hist, length_hist, tmpHist, tmp_length_hist);
      lnZ.push_back(lnZ[0] +
                    estimate_lnZ_from_hist<NVT>(Hist, Betas[0], Betas[1]));
      calc_logDOS_full<NVT>(Hist, std::log(length_hist), Betas, devmax, 1, lnZ, logDOS);
      continue;
    }
    assert(Betas.size() == lnZ.size() + 1);
    assert(lnZ.size() > 1);
    assert(length_hists.size() > 1);
    merge_histograms(Hist, length_hist, tmpHist, tmp_length_hist);
    lnZ.push_back(calc_lnZ<NVT>(logDOS, Betas.back()));
    assert(lnZ.size() > 2);
    assert(length_hists.size() > 2);
    calc_logDOS_full<NVT>(Hist, std::log(length_hist), Betas, devmax, 1,
                          lnZ, logDOS);
  }

  std::cout << "last full WHAM" << "\n";
  calc_logDOS_full<NVT>(Hist, std::log(length_hist), Betas, devmax, 10,
                        lnZ, logDOS);

  dlib::serialize(path + "/Hist.obj") << Hist;
  dlib::serialize(path + "/lnZ.obj") << lnZ;
  dlib::serialize(path + "/logDOS.obj") << logDOS;
  dlib::serialize(path + "/length_hist.obj") << length_hist;
  dlib::serialize(path + "/length_hists.obj") << length_hists;
}
