#include <string>
#include <iostream>
#include <chrono>
#include "../../src/data/load_file.hpp"
#include "../../src/data/all_matching.hpp"
#include "../../src/vector/helper.hpp"
#include "../../src/DiscreteAxis2D.hpp"
#include "../../src/WHAM2D/WHAM2D.hpp"
#include "../../src/WHAM2D/reweight.hpp"

#include "boost/filesystem.hpp"
#include "boost/regex.hpp"
#include "boost/progress.hpp"

using namespace consus;
using namespace consus::WHAM2D;
typedef std::chrono::high_resolution_clock myclock;

int main(int argc, char const *argv[])
{
  if (argc != 4){
    std::cerr << "wrong number of arguments\n";
    std::exit(1);
  }

  std::string path("analysis");
  
  DiscreteAxis2D DOS;
  std::vector<double> length_Hists;
  std::vector<DiscreteAxis2D> MicroMeans;
  std::vector<std::string> header;

  dlib::deserialize(path + "/logDOS.obj") >> DOS;
  dlib::deserialize(path + "/header.obj") >> header;
  dlib::deserialize(path + "/MicroMeans.obj") >> MicroMeans;


  const double Temp= std::stod(argv[1]);
  const double Kappa = std::stod(argv[2]);
  const double step = std::stod(argv[3]);
  auto Parameter = std::make_pair(1.0/Temp, Kappa);
  
  std::cout << "Temperature " << Temp << " Kappa " << Kappa << "\n";

  std::vector<double> bounds;
  auto Hist = reweight_hist<NVT>(DOS, step, Parameter);
  std::string pathr = "results/histograms/2D";
  boost::filesystem::create_directories(pathr);

  std::string filename = pathr + "/T" + std::to_string(Temp) + "_K" + std::to_string(Kappa) + ".dat";
  std::ofstream out(filename, std::ios::trunc);
  for (int i = 0; i < Hist.get_num_bins(); ++i){
      out << Hist.get_value(i) << " " << Hist[i] << "\n";
  }
}
