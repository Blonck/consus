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
  std::vector<double> length_Hists;
  std::vector<DiscreteAxis2D> MicroMeans;
  std::vector<std::string> header;

  dlib::deserialize(path + "/logDOS.obj") >> DOS;
  dlib::deserialize(path + "/header.obj") >> header;
  dlib::deserialize(path + "/MicroMeans.obj") >> MicroMeans;


  const int column = std::stoi(argv[1]);
  const double KappaStart = std::stod(argv[2]);
  const double KappaEnd = std::stod(argv[3]);
  const double KappaStep = std::stod(argv[4]);
  const double Temp = std::stod(argv[5]);
  const double Beta = 1.0/Temp;

  std::cout << "reweight from Kappa " << KappaStart << " to " << KappaEnd << "\n";

  vec1d Kappas;
  vec1<std::pair<double, double>> Parameters;

  for (double Kappa = KappaStart; Kappa <= KappaEnd; Kappa += KappaStep) {
    Kappas.push_back(Kappa);
    Parameters.push_back(std::make_pair(1.0/Temp, Kappa));
  }

  std::cout << header << "\n";

  std::cout << "reweighting " << header[column] << "\n";
  vec1d results(Kappas.size());
  results = reweight_dT<NVT>(DOS, MicroMeans[column], Parameters);

  vec1d jk_means(Kappas.size(), 0.0);
  vec1d jk_error(Kappas.size(), 0.0);
  boost::regex filter(".*[0-9]+");
  auto jk_paths = all_matching_paths(path + "/JK/", filter);
  std::sort(jk_paths.begin(), jk_paths.end());
  size_t num_bins = jk_paths.size();
  vec2d results_jk(Kappas.size(), vec1d(num_bins, 0));
  std::cout << "error calculation\n";
  boost::progress_display progress2(num_bins);
  #pragma omp parallel for
  for (size_t i = 0; i < num_bins; ++i) {
    ++progress2;
    DiscreteAxis2D jk_DOS;
    std::vector<DiscreteAxis2D> jk_MicroMeans;
    dlib::deserialize(jk_paths[i] + "/logDOS.obj") >> jk_DOS;
    dlib::deserialize(jk_paths[i] + "/MicroMeans.obj") >> jk_MicroMeans;
    for (size_t j = 0; j < Kappas.size(); ++j) {
      results_jk[j][i] =
          reweight_dT<NVT>(jk_DOS, jk_MicroMeans[column], {Beta, Kappas[j]});
    }
  }

  for (size_t i = 0; i < results_jk.size(); ++i) {
    for (size_t j = 0; j < results_jk[i].size(); ++j) {
      jk_means[i] += results_jk[i][j];
      jk_error[i] += results_jk[i][j] * results_jk[i][j];
    }
  }
  for (size_t i = 0; i < jk_means.size(); ++i) {
    jk_means[i] = jk_means[i] / static_cast<double>(num_bins);
    jk_error[i] -= static_cast<double>(num_bins) * jk_means[i] * jk_means[i];
    jk_error[i] *=
        static_cast<double>(num_bins - 1) / static_cast<double>(num_bins);
    jk_error[i] = std::sqrt(jk_error[i]);
  }

  std::string pathr = "results";
  boost::filesystem::create_directories(pathr);

  std::string path_T = pathr + "/Temperature";
  boost::filesystem::create_directories(path_T);

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

  std::string path_files = path_T + "/" + std::to_string(Temp);
  boost::filesystem::create_directories(path_files);
  std::string filename = path_files + "/" + replacediv(header[column]) + "_dT.dat";
  std::ofstream out(filename, std::ios::trunc);
  out << "#Kappa " << header[column] << " "
      << "error(" << header[column] << ")\n";
  for (size_t j = 0; j < Kappas.size(); ++j) {
    out << Kappas[j] << " " << results[j] << " " << jk_error[j] << "\n";
  }
  out.close();
}
