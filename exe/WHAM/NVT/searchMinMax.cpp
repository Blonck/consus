
#include <string>
#include <iostream>
#include <chrono>
#include <regex>
#include "../../../src/search/where_diff.hpp"
#include "../../../src/search/find_zero_brent.hpp"
#include "../../../src/data/load_file.hpp"
#include "../../../src/data/all_matching.hpp"
#include "../../../src/vector/helper.hpp"
#include "../../../src/DiscreteAxis.hpp"
#include "../../../src/WHAM/WHAM.hpp"
#include "../../../src/WHAM/reweight.hpp"

#include <omp.h>
#include "boost/filesystem.hpp"
#include "boost/progress.hpp"

using namespace consus;
using namespace consus::WHAM;
typedef std::chrono::high_resolution_clock myclock;

int main(int argc, char const *argv[])
{
  if (argc != 5){
    std::cerr << "wrong number of arguments\n";
    std::exit(1);
  }

  std::string path("analysis");
  
  DiscreteAxis DOS;
  vec1d length_Hists;
  std::vector<DiscreteAxis> MicroMeans;
  std::vector<std::string> header;

  dlib::deserialize(path + "/logDOS.obj") >> DOS;
  dlib::deserialize(path + "/header.obj") >> header;
  dlib::deserialize(path + "/MicroMeans.obj") >> MicroMeans;

  const size_t column =std::stoi(argv[1]);
  if (column > header.size()){
    std::cerr << "column doesn't exists" << "\n";
    std::exit(1);
  }
  std::cout << "reweighting " << header[column] << "\n";

  const double TempStart = std::stod(argv[2]);
  const double TempEnd = std::stod(argv[3]);
  const double TempStep = std::stod(argv[4]);

  std::cout << "intital Temp range " << TempStart << " to " << TempEnd << "\n";

  vec1d Temperatures;
  vec1d Betas;

  for (double Temp = TempStart; Temp <= TempEnd; Temp += TempStep){
    Temperatures.push_back(Temp);
    Betas.push_back(1.0/Temp);
  }

  auto results = reweight_dT2<NVT>(DOS, MicroMeans[column], Betas);
  for (auto& r: results){
    std::cout << r << " ";
  }
  std::cout << "\n";

  vec1d zeros;
  vec1d dT2;
  vec1d dT;
  vec1d value;
  std::vector<std::string> kind;

  auto change_sign = [](double a, double b) -> bool {
    return (std::signbit(a) != std::signbit(b));
  };
  auto pre_zeros = where_diff(results, change_sign);
  for (auto& z: pre_zeros){
    std::cout << z << " ";
  }
  std::cout << "\n";

  auto calc_dT2 = [&](const double beta) {
    return reweight_dT2<NVT>(DOS, MicroMeans[column], beta);
  };
  auto calc_dT = [&](const double beta) {
    return reweight_dT<NVT>(DOS, MicroMeans[column], beta);
  };
  auto calc_value = [&](const double beta) {
    return reweight<NVT>(DOS, MicroMeans[column], beta);
  };
  for (size_t j = 0; j < pre_zeros.size(); ++j) {
    auto ret =
        find_zero_brent(calc_dT2, Betas[pre_zeros[j]], Betas[pre_zeros[j] + 1]);
    if (ret.valid) {
      double tzero = 1.0 / ret.x;
      std::cout << "found " << tzero << "\n";
      zeros.push_back(tzero);
      dT2.push_back(ret.fx);
      value.push_back(calc_value(ret.x));
      auto dtval = calc_dT(ret.x);
      auto before = calc_dT(1.0 / (tzero - 0.1 * TempStep));
      auto after = calc_dT(1.0 / (tzero + 0.1 * TempStep));
      dT.push_back(dtval);
      if ((dtval > 0) and (before < dtval) and (after < dtval)) {
        kind.push_back("MAX");
      } else {
        if ((dtval < 0) and (before > dtval) and (after > dtval)) {
          kind.push_back("MIN");
        } else {
          kind.push_back("NONE");
        }
      }
    } 
  }

  vec1d means_zero(zeros.size(), 0.0);
  vec1d error_zero(zeros.size(), 0.0);
  vec1d means_value(zeros.size(), 0.0);
  vec1d error_value(zeros.size(), 0.0);
  vec1d means_dT(zeros.size(), 0.0);
  vec1d error_dT(zeros.size(), 0.0);
  vec1d means_dT2(zeros.size(), 0.0);
  vec1d error_dT2(zeros.size(), 0.0);
  std::vector<size_t> num_bins_zero(zeros.size(), 0.0);

  boost::regex filter(".*[0-9]+");
  auto jk_paths = all_matching_paths(path + "/JK/", filter);
  std::sort(jk_paths.begin(), jk_paths.end());
  size_t num_bins = jk_paths.size();
  std::cout << "error analysis"
            << "\n";
  if (jk_paths.size() > 0) {
#pragma omp parallel for
    for (size_t i = 0; i < zeros.size(); ++i) {
      std::cout << zeros[i] << "\n";
      for (size_t j = 0; j < num_bins; ++j) {
        //++progress3;
        DiscreteAxis DOS;
        std::vector<DiscreteAxis> MicroMeans;
        std::vector<std::string> header;
        dlib::deserialize(jk_paths[j] + "/logDOS.obj") >> DOS;
        dlib::deserialize(jk_paths[j] + "/MicroMeans.obj") >> MicroMeans;
        auto calc_dT2 = [&](const double beta) {
          return reweight_dT2<NVT>(DOS, MicroMeans[column], beta);
        };
        auto calc_dT = [&](const double beta) {
          return reweight_dT<NVT>(DOS, MicroMeans[column], beta);
        };
        auto calc_value = [&](const double beta) {
          return reweight<NVT>(DOS, MicroMeans[column], beta);
        };
        double lower_bound = zeros[i] - TempStep;
        double min_lower_bound = 0.0;
        if (i != 0) {
          min_lower_bound = zeros[i - 1];
        }
        lower_bound = std::max(lower_bound, min_lower_bound);
        double upper_bound = zeros[i] + TempStep;
        double max_upper_bound = TempEnd + TempStep;
        if (i != zeros.size() - 1) {
          max_upper_bound = zeros[i + 1];
        }
        upper_bound = std::min(upper_bound, max_upper_bound);
        // std::cout << "starting brent with [" << lower_bound << ", " <<
        // upper_bound << "]\n";
        // std::cout << "bound limits [" << min_lower_bound << ", " <<
        // max_upper_bound << "]\n";
        auto ret =
            find_zero_brent(calc_dT2, 1.0 / lower_bound, 1.0 / upper_bound,
                            1.0 / min_lower_bound, 1.0 / max_upper_bound);
        auto tzero = 1.0 / ret.x;
        // std::cout << "zero at " << tzero << " with " << ret.fx << " is "
        //         << ((ret.valid == true) ? "valid" : "not valid") << "\n";
        if (ret.valid) {
          // std::cout << "jk " << j << " " << tzero << "\n";
          auto val = calc_value(ret.x);
          auto dtval = calc_dT(ret.x);
          means_zero[i] += tzero;
          error_zero[i] += tzero * tzero;
          means_value[i] += val;
          error_value[i] += val * val;
          means_dT[i] += dtval;
          error_dT[i] += dtval * dtval;
          means_dT2[i] += ret.fx;
          error_dT2[i] += ret.fx * ret.fx;
          num_bins_zero[i]++;
        }
      }
    }
  }

  for (size_t i = 0; i < zeros.size(); ++i) {
    means_zero[i] = means_zero[i] / static_cast<double>(num_bins_zero[i]);
    error_zero[i] -=
        static_cast<double>(num_bins_zero[i]) * means_zero[i] * means_zero[i];
    error_zero[i] *= static_cast<double>(num_bins_zero[i] - 1) /
                     static_cast<double>(num_bins_zero[i]);
    error_zero[i] = std::sqrt(error_zero[i]);
    means_value[i] = means_value[i] / static_cast<double>(num_bins_zero[i]);
    error_value[i] -=
        static_cast<double>(num_bins_zero[i]) * means_value[i] * means_value[i];
    error_value[i] *= static_cast<double>(num_bins_zero[i] - 1) /
                      static_cast<double>(num_bins_zero[i]);
    error_value[i] = std::sqrt(error_value[i]);
    means_dT[i] = means_dT[i] / static_cast<double>(num_bins_zero[i]);
    error_dT[i] -=
        static_cast<double>(num_bins_zero[i]) * means_dT[i] * means_dT[i];
    error_dT[i] *= static_cast<double>(num_bins_zero[i] - 1) /
                   static_cast<double>(num_bins_zero[i]);
    error_dT[i] = std::sqrt(error_dT[i]);
    means_dT2[i] = means_dT2[i] / static_cast<double>(num_bins_zero[i]);
    error_dT2[i] -=
        static_cast<double>(num_bins_zero[i]) * means_dT2[i] * means_dT2[i];
    error_dT2[i] *= static_cast<double>(num_bins_zero[i] - 1) /
                    static_cast<double>(num_bins_zero[i]);
    error_dT2[i] = std::sqrt(error_dT2[i]);
  }

  path = "results/Maximas";
  boost::filesystem::create_directories(path);

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

    std::string name = path + "/" + replacediv(header[column]) + "_dT.dat";
    std::ofstream file(name, std::ios::out | std::ios::trunc);
    file << "#Maximas of " << header[column] << "_dT\n";
    file << "#Temperature err(Temperature) " << header[column] << " err("
         << header[column] << ") " << header[column] << "_dT err("
         << header[column] << "_dT) " << header[column] << "_dT2 err("
         << header[column] << "_dT) kind\n";
    for (size_t j = 0; j < zeros.size(); ++j) {
      file << zeros[j] << " " << error_zero[j] << " " << value[j] << " "
           << error_value[j] << " " << dT[j] << " " << error_dT[j] << " "
           << dT2[j] << " " << error_dT2[j] << " " << kind[j] << "\n";
    }

}
