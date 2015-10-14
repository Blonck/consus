#include <string>
#include <iostream>
#include <chrono>
#include <regex>
#include "../../src/data/load_file.hpp"
#include "../../src/data/all_matching.hpp"
#include "../../src/vector/helper.hpp"
#include "../../src/DiscreteAxis2D.hpp"
#include "../../src/WHAM2D/WHAM2D.hpp"

#include "boost/filesystem.hpp"

using namespace consus;
using namespace consus::WHAM2D;

int main(int argc, char const *argv[])
{
  if (argc != 8){
    std::cerr << "wrong number of arguments\n";
    std::exit(1);
  }

  std::string path("analysis");
  
  DiscreteAxis2D DOS;
  DiscreteAxis2D logHist;
  vec1d length_Hists;
  std::vector<DiscreteAxis2D> MicroMeans;
  std::vector<std::string> header;

  dlib::deserialize(path + "/logDOS.obj") >> DOS;
  dlib::deserialize(path + "/header.obj") >> header;
  dlib::deserialize(path + "/MicroMeans.obj") >> MicroMeans;

  const size_t column = std::stoi(argv[1]);
  if (column > header.size()) {
    std::cerr << "column doesn't exists"
              << "\n";
    std::exit(1);
  }
  std::cout << "reweighting " << header[column] << "\n";

  const double TempStart = std::stod(argv[2]);
  const double TempEnd = std::stod(argv[3]);
  const int    NumTemp = std::stoi(argv[4]);
  const double KappaStart = std::stod(argv[5]);
  const double KappaEnd = std::stod(argv[6]);
  const double KappaStep = std::stod(argv[7]);

  std::cout << "reweight from Temp " << TempStart << " to " << TempEnd << "\n";
  std::cout << "reweight from Kappa " << KappaStart << " to " << KappaEnd << "\n";

  vec1d Temperatures;
  vec1<std::pair<double, double>> Parameters;
  size_t NumKappas = 0;

  double logMin = std::log(TempStart);
  double logMax = std::log(TempEnd);
  double factor = exp((logMax - logMin) / NumTemp);
  for (double Kappa = KappaStart; Kappa <= KappaEnd; Kappa += KappaStep) {
    ++NumKappas;
    Temperatures.push_back(TempStart);
    Parameters.push_back({1.0 / TempStart, Kappa});
    for (int i = 1; i < NumTemp; ++i) {
      double Temp = Temperatures.back() * factor;
      Temperatures.push_back(Temp);
      Parameters.push_back({1.0 / Temp, Kappa});
    }
  }
  vec1d results(Parameters.size(), 0);

  results = reweight_dT2<NVT>(DOS, MicroMeans[column], Parameters);

  path = "results";
  boost::filesystem::create_directories(path);
  

  std::string path_kappa = path + "/Kappa";
  boost::filesystem::create_directories(path_kappa);

  std::string path_temp = path + "/Temperature";
  boost::filesystem::create_directories(path_temp);

  auto replacediv = [](const std::string name){
    std::string ret = name;
    auto pos = name.find( '/' );
    if (pos != std::string::npos){
      ret = name.substr(0, pos);
      ret = ret + "div";
      ret = ret + name.substr(pos+1);
    }
    return ret;
  };
  
  std::string path_contour = path + "/contour";
  boost::filesystem::create_directories(path_contour);

  std::string filename = path_contour + "/" + replacediv(header[column]) + "_dT.dat";
  std::ofstream out(filename, std::ios::trunc);
  out << "#Kappa Temperature " << header[column] << "_dT\n";
  for (size_t k = 0; k < Parameters.size(); ++k) {
    out << Parameters[k].second << " " << Temperatures[k] << " " << results[k]
        << "\n";
    if ((k % NumTemp) == (NumTemp - 1)) {
      out << "\n";
    }
  }
}
