
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

  std::vector<std::string> filenames;
  std::vector<double> Parameters;
  DiscreteAxis Hist;
  double length_hist;
  std::vector<double> length_hists;
  vec1<double> length_Hists;
  std::vector<double> lnZ;
  DiscreteAxis logDOS;
  Range range;
  dlib::deserialize(path + "/range.obj") >> range;
  dlib::deserialize(path + "/filenames.obj") >> filenames;
  dlib::deserialize(path + "/lnZ.obj") >> lnZ;
  dlib::deserialize(path + "/logDOS.obj") >> logDOS;
  dlib::deserialize(path + "/Parameters.obj") >> Parameters;
  dlib::deserialize(path + "/Hist.obj") >> Hist;
  dlib::deserialize(path + "/length_hist.obj") >> length_hist;
  dlib::deserialize(path + "/length_hists.obj") >> length_hists;

  std::cout << "refine DOS" << "\n";
  calc_logDOS_full<NVT>(Hist, std::log(length_hist), Parameters, devmax, 10, lnZ,
                        logDOS);

  dlib::serialize(path + "/lnZ.obj") << lnZ;
  dlib::serialize(path + "/logDOS.obj") << logDOS;
}
