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
  
  DiscreteAxis2D DOS;
  std::vector<DiscreteAxis2D> MicroMeans;
  std::vector<std::string> header;

  dlib::deserialize(path + "/logDOS.obj") >> DOS;
  dlib::deserialize(path + "/header.obj") >> header;
  dlib::deserialize(path + "/MicroMeans.obj") >> MicroMeans;


  const double TempStart = std::stod(argv[1]);
  const double TempEnd = std::stod(argv[2]);
  const double TempStep = std::stod(argv[3]);
  const double Kappa = std::stod(argv[4]);

  const size_t column = std::stod(argv[5]);

  std::cout << "reweight from Temp " << TempStart << " to " << TempEnd << "\n";

  vec1d Temperatures;
  vec1d Betas;
  vec1<std::pair<double, double>> Parameters;

  for (double Temp = TempStart; Temp <= TempEnd; Temp += TempStep){
    Temperatures.push_back(Temp);
    Betas.push_back(1.0/Temp);
    Parameters.push_back(std::make_pair(1.0/Temp, Kappa));
  }

  std::cout << header << "\n";

  vec1d results(Temperatures.size(), 0.0);
  std::cout << "reweighting " << header[column] << "\n";
  results = reweight_dT<NVT>(DOS, MicroMeans[column], Parameters);

  boost::regex filter(".*[0-9]+");
  auto jk_paths = all_matching_paths(path + "/JK/", filter);
  std::sort(jk_paths.begin(), jk_paths.end());
  size_t num_bins = jk_paths.size();
  std::cout << "calculation of jk blocks\n";
  boost::progress_display progress2(num_bins);
  vec2d results_jk(num_bins);
  #pragma omp parallel for
  for (size_t i = 0; i < num_bins; ++i) {
    DiscreteAxis2D jk_DOS;
    std::vector<DiscreteAxis2D> jk_MicroMeans;
    dlib::deserialize(jk_paths[i] + "/logDOS.obj") >> jk_DOS;
    dlib::deserialize(jk_paths[i] + "/MicroMeans.obj") >> jk_MicroMeans;
    results_jk[i] = reweight_dT<NVT>(jk_DOS, jk_MicroMeans[column], Parameters);
    ++progress2;
  }


  std::cout << "error calculation\n";
  vec1d jk_means(Temperatures.size(), 0.0);
  vec1d jk_error(Temperatures.size(), 0.0);
  for (size_t i = 0; i < Betas.size(); ++i) {
    for (size_t j = 0; j < num_bins; ++j) {
      jk_means[i] += results_jk[j][i];
      jk_error[i] += results_jk[j][i] * results_jk[j][i];
    }
  }
  for (size_t i = 0; i < jk_means.size(); ++i) {
    jk_means[i] = jk_means[i] / static_cast<double>(num_bins);
    jk_error[i] -=
        static_cast<double>(num_bins) * jk_means[i] * jk_means[i];
    jk_error[i] *=
        static_cast<double>(num_bins - 1) / static_cast<double>(num_bins);
    jk_error[i] = std::sqrt(jk_error[i]);
  }

  std::string pathr = "results";
  boost::filesystem::create_directories(pathr);

  std::string path_kappa = pathr + "/Kappa";
  boost::filesystem::create_directories(path_kappa);

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

  std::string path_files = path_kappa + "/" + std::to_string(Kappa);
  boost::filesystem::create_directories(path_files);
  std::string filename =
      path_files + "/" + replacediv(header[column]) + "_dT.dat";
  std::cout << "write " << filename << "\n";
  std::ofstream out(filename, std::ios::trunc);
  out << "#Temperature " << header[column] << "_dT "
      << "error(" << header[column] << "_dT)\n";
  for (size_t j = 0; j < Temperatures.size(); ++j) {
    out << Temperatures[j] << " " << results[j] << " " << jk_error[j] << "\n";
  }
  out.close();
}
