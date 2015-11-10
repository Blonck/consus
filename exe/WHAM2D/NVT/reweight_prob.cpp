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
  if (argc < 5){
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


  const int column = std::stoi(argv[1]);
  const double Temp= std::stod(argv[2]);
  const double Kappa = std::stod(argv[3]);
  auto Parameter = std::make_pair(1.0/Temp, Kappa);
  
  std::cout << "Temperature " << Temp << " Kappa " << Kappa << "\n";

  std::vector<double> bounds;

  if (argv[4] == std::string("-e")){
    for (int i = 5; i < argc; ++i){
      bounds.push_back(std::stod(argv[i]));
    }
  }else if(argv[4] == std::string("-l")){
    auto min = std::stod(argv[5]);
    auto max = std::stod(argv[6]);
    auto step = std::stod(argv[7]);
    for (auto value = min; value < max; value += step){
      bounds.push_back(value);
    }
  }else{
    std::cerr << "4th argument either -e (explicit) or -l (linspace) \n";
    std::exit(1);
  }


  std::cout << header << "\n";
  std::cout << "reweighting " << header[column] << "\n";
  auto results = reweight_prob<NVT>(DOS, MicroMeans[column], Parameter, bounds);

  auto replacediv = [](const std::string name) {
    std::string ret = name;
    auto pos = name.find('/');
    if (pos != std::string::npos) {
      ret = name.substr(0, pos);
      ret = ret + "div";
      ret = ret + name.substr(pos + 1);
    }
    return ret;
  };

  std::string pathr = "results/Probs/" + replacediv(header[column]);
  boost::filesystem::create_directories(pathr);

  std::cout << results << "\n";
  std::string filename = pathr + "/T" + std::to_string(Temp) + "_K" + std::to_string(Kappa) + ".dat";
  std::ofstream out(filename, std::ios::trunc);
  out << "#value p\n";
  for (size_t i = 0; i < results.size()-1; ++i){
    out << bounds[i] << " " << results[i] << "\n";
  }
  out.close();
}
