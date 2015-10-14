#include <string>
#include <iostream>
#include <chrono>
#include <regex>
#include "../../src/search/where_diff.hpp"
#include "../../src/search/find_zero_brent.hpp"
#include "../../src/data/load_file.hpp"
#include "../../src/data/all_matching.hpp"
#include "../../src/vector/helper.hpp"
#include "../../src/DiscreteAxis2D.hpp"
#include "../../src/WHAM2D/WHAM2D.hpp"

#include <omp.h>
#include "boost/filesystem.hpp"
#include "boost/progress.hpp"

using namespace consus;
using namespace consus::WHAM2D;
typedef std::chrono::high_resolution_clock myclock;

int main(int argc, char const *argv[])
{
  if (argc != 8){
    std::cerr << "wrong number of arguments\n";
    std::exit(1);
  }

  std::string path("analysis");
  
  DiscreteAxis2D DOS;
  vec1d length_Hists;
  std::vector<DiscreteAxis2D> MicroMeans;
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
  const double KappaStart = std::stod(argv[5]);
  const double KappaEnd = std::stod(argv[6]);
  const double KappaStep = std::stod(argv[7]);

  std::cout << "intital Temp range " << TempStart << " to " << TempEnd << "\n";
  std::cout << "search in Kappa from " << KappaStart << " to " << KappaEnd << "\n";

  vec1d Temperatures;
  vec1d Betas;
  vec1d Kappas;

  for (double Temp = TempStart; Temp <= TempEnd; Temp += TempStep){
    Temperatures.push_back(Temp);
    Betas.push_back(1.0/Temp);
  }
  for (double Kappa = KappaStart; Kappa <= KappaEnd; Kappa += KappaStep){
    Kappas.push_back(Kappa);
  }
  vec2d results(Kappas.size(), vec1d(Temperatures.size()));

  std::cout << header << "\n";

  for (size_t i = 0; i < Kappas.size(); ++i) {
    for (size_t j = 0; j < Betas.size(); ++j){
      double beta = Betas[i];
      double kappa = Kappas[j];
      results[i][j] = reweight_dT2<NVT>(DOS, MicroMeans[column], {beta, kappa});
    }
  }

  vec2d zeros(Kappas.size());
  vec2d dT2(Kappas.size());
  vec2d dT(Kappas.size());
  vec2d value(Kappas.size());
  std::vector< std::vector<std::string> > kind(Kappas.size());

#pragma omp parallel for
  for (size_t i = 0; i < Kappas.size(); ++i) {
    std::cout << Kappas[i] << "\n";
    auto change_sign = [](double a, double b)->bool {
      return (std::signbit(a) != std::signbit(b));
    };
    auto pre_zeros = where_diff(results[i], change_sign);
    auto calc_dT2 = [&](const double beta) {
      return reweight_dT2<NVT>(DOS, MicroMeans[column], {beta, Kappas[i]});
    };
    auto calc_dT = [&](const double beta) {
      return reweight_dT<NVT>(DOS, MicroMeans[column], {beta, Kappas[i]});
    };
    auto calc_value = [&](const double beta) {
      return reweight<NVT>(DOS, MicroMeans[column], {beta, Kappas[i]});
    };
    for (size_t j = 0; j < pre_zeros.size(); ++j) {
      auto ret = find_zero_brent(calc_dT2, Betas[pre_zeros[j]],
                                 Betas[pre_zeros[j] + 1]);
      if (ret.valid) {
        double tzero = 1.0 / ret.x;
        std::cout << "found " << tzero << "\n";
        zeros[i].push_back(tzero);
        dT2[i].push_back(ret.fx);
        value[i].push_back(calc_value(ret.x));
        auto dtval = calc_dT(ret.x);
        auto before = calc_dT(1.0 / (tzero - 0.1*TempStep));
        auto after = calc_dT(1.0 / (tzero + 0.1*TempStep));
        dT[i].push_back(dtval);
        if ( (dtval > 0) and (before < dtval) and (after < dtval)){
          kind[i].push_back("MAX");
        }else{
          if ( (dtval < 0) and (before > dtval) and (after > dtval)){
            kind[i].push_back("MIN");
          }else{
            kind[i].push_back("NONE");
          }
        }
      } 
    }
  }

  vec2d means_zero(Kappas.size());
  vec2d error_zero(Kappas.size());
  vec2d means_value(Kappas.size());
  vec2d error_value(Kappas.size());
  vec2d means_dT(Kappas.size());
  vec2d error_dT(Kappas.size());
  vec2d means_dT2(Kappas.size());
  vec2d error_dT2(Kappas.size());
  std::vector< std::vector<size_t> > num_bins_zero(Kappas.size());
  for (size_t i = 0; i < Kappas.size(); ++i) {
    means_zero[i].resize(zeros[i].size(), 0.0);
    error_zero[i].resize(zeros[i].size(), 0.0);
    means_value[i].resize(zeros[i].size(), 0.0);
    error_value[i].resize(zeros[i].size(), 0.0);
    means_dT[i].resize(zeros[i].size(), 0.0);
    error_dT[i].resize(zeros[i].size(), 0.0);
    means_dT2[i].resize(zeros[i].size(), 0.0);
    num_bins_zero[i].resize(zeros[i].size(), 0.0);
    error_dT2[i].resize(zeros[i].size(), 0.0);
  }

    boost::regex filter(".*[0-9]+");
    auto jk_paths = all_matching_paths(path + "/JK/", filter);
    std::sort(jk_paths.begin(), jk_paths.end());
    size_t num_bins = jk_paths.size();
    std::cout << "error analysis" << "\n";
    if (jk_paths.size() > 0) {
      #pragma omp parallel for
      for (size_t k = 0; k < Kappas.size(); ++k) {
        std::cout << Kappas[k] << "\n";
        for (size_t i = 0; i < zeros[k].size(); ++i) {
          std::cout << zeros[k][i] << "\n";
          for (size_t j = 0; j < num_bins; ++j) {
            //std::cout << "bin " << j << "\n";
            //++progress3;
            DiscreteAxis2D DOS;
            std::vector<DiscreteAxis2D> MicroMeans;
            std::vector<std::string> header;
            dlib::deserialize(jk_paths[j] + "/logDOS.obj") >> DOS;
            dlib::deserialize(jk_paths[j] + "/header.obj") >> header;
            dlib::deserialize(jk_paths[j] + "/MicroMeans.obj") >> MicroMeans;
            auto calc_dT2 = [&](const double beta) {
              return reweight_dT2<NVT>(DOS, MicroMeans[column], {beta, Kappas[k]});
            };
            auto calc_dT = [&](const double beta) {
              return reweight_dT<NVT>(DOS, MicroMeans[column], {beta, Kappas[k]});
            };
            auto calc_value = [&](const double beta) {
              return reweight<NVT>(DOS, MicroMeans[column], {beta, Kappas[k]});
            };
            double lower_bound = zeros[k][i] - TempStep;
            double min_lower_bound = 0.0;
            if (i != 0) {
              min_lower_bound = zeros[k][i - 1];
            }
            lower_bound = std::max(lower_bound, min_lower_bound);
            double upper_bound = zeros[k][i] + TempStep;
            double max_upper_bound = TempEnd + TempStep;
            if (i != zeros[k].size() - 1) {
              max_upper_bound = zeros[k][i + 1];
            }
            upper_bound = std::min(upper_bound, max_upper_bound);
            //std::cout << "starting brent with [" << lower_bound << ", " <<
            //upper_bound << "]\n";
            //std::cout << "bound limits [" << min_lower_bound << ", " <<
            //max_upper_bound << "]\n";
            auto ret =
                find_zero_brent(calc_dT2, 1.0 / lower_bound, 1.0 / upper_bound,
                                1.0 / min_lower_bound, 1.0 / max_upper_bound);
            auto tzero = 1.0 / ret.x;
            //std::cout << "zero at " << tzero << " with " << ret.fx << " is "
            //         << ((ret.valid == true) ? "valid" : "not valid") << "\n";
            if (ret.valid) {
              //std::cout << "jk " << j << " " << tzero << "\n";
              auto val = calc_value(ret.x);
              auto dtval = calc_dT(ret.x);
              means_zero[k][i] += tzero;
              error_zero[k][i] += tzero * tzero;
              means_value[k][i] += val;
              error_value[k][i] += val * val;
              means_dT[k][i] += dtval;
              error_dT[k][i] += dtval * dtval;
              means_dT2[k][i] += ret.fx;
              error_dT2[k][i] += ret.fx * ret.fx;
              num_bins_zero[k][i]++;
            }
          }
        }
      }
    }

    for (size_t k = 0; k < Kappas.size(); ++k) {
      for (size_t i = 0; i < zeros[k].size(); ++i) {
        means_zero[k][i] =
            means_zero[k][i] / static_cast<double>(num_bins_zero[k][i]);
        error_zero[k][i] -= static_cast<double>(num_bins_zero[k][i]) *
                            means_zero[k][i] * means_zero[k][i];
        error_zero[k][i] *= static_cast<double>(num_bins_zero[k][i] - 1) /
                            static_cast<double>(num_bins_zero[k][i]);
        error_zero[k][i] = std::sqrt(error_zero[k][i]);
        means_value[k][i] =
            means_value[k][i] / static_cast<double>(num_bins_zero[k][i]);
        error_value[k][i] -= static_cast<double>(num_bins_zero[k][i]) *
                             means_value[k][i] * means_value[k][i];
        error_value[k][i] *= static_cast<double>(num_bins_zero[k][i] - 1) /
                             static_cast<double>(num_bins_zero[k][i]);
        error_value[k][i] = std::sqrt(error_value[k][i]);
        means_dT[k][i] =
            means_dT[k][i] / static_cast<double>(num_bins_zero[k][i]);
        error_dT[k][i] -= static_cast<double>(num_bins_zero[k][i]) *
                          means_dT[k][i] * means_dT[k][i];
        error_dT[k][i] *= static_cast<double>(num_bins_zero[k][i] - 1) /
                          static_cast<double>(num_bins_zero[k][i]);
        error_dT[k][i] = std::sqrt(error_dT[k][i]);
        means_dT2[k][i] =
            means_dT2[k][i] / static_cast<double>(num_bins_zero[k][i]);
        error_dT2[k][i] -= static_cast<double>(num_bins_zero[k][i]) *
                           means_dT2[k][i] * means_dT2[k][i];
        error_dT2[k][i] *= static_cast<double>(num_bins_zero[k][i] - 1) /
                           static_cast<double>(num_bins_zero[k][i]);
        error_dT2[k][i] = std::sqrt(error_dT2[k][i]);
      }
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

  for (size_t i = 0; i < Kappas.size(); ++i) {
    std::string path_files = path + "/" + std::to_string(Kappas[i]);
    boost::filesystem::create_directories(path_files);
    std::string name = path_files + "/" + replacediv(header[column]) + "_dT.dat";
    std::ofstream file(name, std::ios::out | std::ios::trunc);
    file << "#Maximas of " << header[column] << "_dT\n";
    file << "#Temperature err(Temperature) " << header[column] << " err("
         << header[column] << ") " << header[column] << "_dT err("
         << header[column] << "_dT) " << header[column] << "_dT2 err("
         << header[column] << "_dT) kind\n";
    for (size_t j = 0; j < zeros[i].size(); ++j) {
      file << zeros[i][j] << " " << error_zero[i][j] << " "
           << value[i][j] << " " << error_value[i][j] << " " << dT[i][j] << " "
           << error_dT[i][j] << " " << dT2[i][j] << " " << error_dT2[i][j]
           << " " << kind[i][j] << "\n";
    }
  }

  std::string name = path + "/" + replacediv(header[column]) + "_dT.dat";
  std::ofstream file(name, std::ios::out | std::ios::trunc);
  file << "#Maximas of " << header[column] << "_dT\n";
  file << "#Kappa Temperature err(Temperature) " << header[column] << " err("
       << header[column] << ") " << header[column] << "_dT err("
       << header[column] << "_dT) " << header[column] << "_dT2 err("
       << header[column] << "_dT) kind\n";
  for (size_t i = 0; i < Kappas.size(); ++i) {
    for (size_t j = 0; j < zeros[i].size(); ++j) {
      file << Kappas[i] << " " << zeros[i][j] << " " << error_zero[i][j] << " "
           << value[i][j] << " " << error_value[i][j] << " " << dT[i][j] << " "
           << error_dT[i][j] << " " << dT2[i][j] << " " << error_dT2[i][j]
           << " " << kind[i][j] << "\n";
    }
  }
}
