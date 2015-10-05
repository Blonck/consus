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

using namespace consus;
using namespace consus::WHAM2D;
typedef std::chrono::high_resolution_clock myclock;

int main(int argc, char const *argv[])
{
  if (argc != 5){
    std::cerr << "wrong number of arguments\n";
    std::exit(1);
  }

  std::string path("analysis");
  
  DiscreteAxis2D DOS;
  DiscreteAxis2D Hist;
  std::vector<double> length_Hists;
  std::vector<DiscreteAxis2D> MicroMeans;
  std::vector<std::string> header;

  dlib::deserialize(path + "/logDOS.obj") >> DOS;
  dlib::deserialize(path + "/header.obj") >> header;
  dlib::deserialize(path + "/MicroMeans.obj") >> MicroMeans;

  const double TempStart = std::stod(argv[1]);
  const double TempEnd = std::stod(argv[2]);
  const double TempStep = std::stod(argv[3]);
  const double Kappa = std::stod(argv[4]);

  std::cout << "reweight from Temp " << TempStart << " to " << TempEnd << "\n";

  vec1d Temperatures;
  vec1d Betas;

  for (double Temp = TempStart; Temp <= TempEnd; Temp += TempStep){
    Temperatures.push_back(Temp);
    Betas.push_back(1.0/Temp);
  }

  std::cout << header << "\n";

  for (size_t column = 0; column < header.size(); ++column) {
    std::cout << "reweighting " << header[column] << "\n";
    vec1d results(Temperatures.size());
    #pragma omp parallel for
    for (size_t j = 0; j < Betas.size(); ++j) {
      results[j] = reweight_dT2<NVT>(DOS, MicroMeans[column], {Betas[j], Kappa});
    }

    vec1d jk_means(Temperatures.size());
    vec1d jk_error(Temperatures.size());
    if (boost::filesystem::is_directory(path + "/JK")) {
      boost::regex filter(".*[0-9]+");
      auto jk_paths = all_matching_paths(path + "/JK/", filter);
      std::sort(jk_paths.begin(), jk_paths.end());
      size_t num_bins = jk_paths.size();
      vec2d results_jk(Temperatures.size(), vec1d(num_bins, 0));
      std::cout << jk_paths << "\n";
      if (jk_paths.size() > 0) {
        for (size_t i = 0; i < num_bins; ++i) {
          DiscreteAxis2D jk_DOS;
          std::vector<DiscreteAxis2D> jk_MicroMeans;
          dlib::deserialize(jk_paths[i] + "/logDOS.obj") >> jk_DOS;
          dlib::deserialize(jk_paths[i] + "/MicroMeans.obj") >> jk_MicroMeans;
          #pragma omp parallel for
          for (size_t j = 0; j < Betas.size(); ++j) {
            results_jk[j][i] =
                reweight_dT2<NVT>(jk_DOS, jk_MicroMeans[column], {Betas[j], Kappa});
          }
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
        jk_error[i] -=
            static_cast<double>(num_bins) * jk_means[i] * jk_means[i];
        jk_error[i] *=
            static_cast<double>(num_bins - 1) / static_cast<double>(num_bins);
        jk_error[i] = std::sqrt(jk_error[i]);
      }
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
        path_files + "/" + replacediv(header[column]) + "_dT2.dat";
    std::ofstream out(filename, std::ios::trunc);
    out << "#Temperature " << header[column] << "_dT2 "
        << "error(" << header[column] << "_dT)\n";
    for (size_t j = 0; j < Temperatures.size(); ++j) {
      out << Temperatures[j] << " " << results[j] << " " << jk_error[j] << "\n";
    }
    out.close();
  }
}
