#include <string>
#include <iostream>
#include <chrono>
#include "../../src/data/load_file.hpp"
#include "../../src/data/all_matching.hpp"
#include "../../src/vector/helper.hpp"
#include "../../src/DiscreteAxis2D.hpp"
#include "../../src/WHAM2D/WHAM2D.hpp"

#include "boost/filesystem.hpp"
#include "boost/regex.hpp"
#include "boost/progress.hpp"

using namespace consus;
using namespace consus::WHAM2D;
typedef std::chrono::high_resolution_clock myclock;

int main(int argc, char const *argv[])
{
  if (argc != 6){
    std::cerr << "wrong number of arguments\n";
    std::exit(1);
  }

  std::string path("analysis");

  const double Kappa = std::stod(argv[1]);
  const double step_E = std::stod(argv[2]);
  const double beta_min = std::stod(argv[3]);
  const double beta_max = std::stod(argv[4]);
  const double step_b = std::stod(argv[5]);
  
  DiscreteAxis2D DOS;

  dlib::deserialize(path + "/logDOS.obj") >> DOS;

  std::string pathr = "results/logWeights";
  boost::filesystem::create_directories(pathr);

  std::string path_kappa = pathr + "/" + std::to_string(Kappa);
  boost::filesystem::create_directories(path_kappa);
  
  std::ofstream file(path_kappa + "/InvBetaMicro.dat", std::ios::trunc);
  file << "#beta_m E_max\n";
  for (double beta = beta_min; beta <= beta_max; beta += step_b) {
    auto Hist = reweight_hist<NVT>(DOS, step_E, {beta, Kappa});
    auto max = std::max_element(Hist.begin(), Hist.end());
    auto bin = std::distance(Hist.begin(), max);
    file << beta << " " << Hist.get_value(bin) << "\n";
  }

}
