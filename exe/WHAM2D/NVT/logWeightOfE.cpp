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
  if (argc != 5){
    std::cerr << "wrong number of arguments\n";
    std::exit(1);
  }

  std::string path("analysis");
  std::string path_jk("analysis/JK");

  const double Kappa = std::stod(argv[1]);
  const double step = std::stod(argv[2]);
  const int dh1 = std::stoi(argv[3]);
  const int dh2 = std::stoi(argv[4]);
  
  DiscreteAxis2D DOS;
  std::vector<double> length_Hists;
  int NumBins;

  dlib::deserialize(path + "/logDOS.obj") >> DOS;
  dlib::deserialize(path + "/NumBins.obj") >> NumBins;
  double minE = DOS.get_value_first(0) + Kappa * DOS.get_value_second(0);
  double maxE = DOS.get_value_first(DOS.get_num_bins_first() - 1) +
                Kappa * DOS.get_value_second(DOS.get_num_bins_second() - 1);

  DiscreteAxis logW(minE, maxE, step, log_zero<double>());
  DiscreteAxis norm(minE, maxE, step, 0.0);
  for (int i = 0; i < DOS.get_num_bins_first(); ++i){
    const double E1 = DOS.get_value_first(i);
    for (int j = 0; j < DOS.get_num_bins_second(); ++j){
      const double E2 = DOS.get_value_second(j);
      const double E = NVT::ham(E1, Kappa, E2);
      const int index = DOS.get_index(i, j);
      const int rindex = logW.get_bin(E);
      if (std::isfinite(DOS[index] and DOS[index] > log_zero<double>() + 1.0e10)){
        logW[rindex] = addlogwise(logW[rindex], DOS[index]);
        norm[rindex] += 1;
      }
    }
  }

  for (int i = 0; i < logW.size(); ++i){
    if (norm[i] != 0.0 and logW[i] > -1.0e100){
      logW[i] -= std::log(norm[i]);
    }else{
      logW[i] = std::numeric_limits<double>::quiet_NaN();
    }
  }

  std::function<double(int, int, const DiscreteAxis&)> df = [&](
      int h, int x, const DiscreteAxis& LogW) {
    double f_;
    if ( 2*h <= x && x <(int)LogW.size()-2*h){
      f_ = (-1.0 * LogW[x + 2 * h] + 8 * LogW[x + h] - 8 * LogW[x - h] +
            LogW[x - 2 * h]) /
           (12.0 * h * LogW.get_step());
    }else if (x < 2*h){
      f_ = df(h, 2 * h, LogW);
    }else{
      f_ = df(h, (int)LogW.size() - 2 * h - 1, LogW);
    }
    return f_;
  };

  std::function<double(int, int, const DiscreteAxis&)> df2 = [&](
      int h, int x, const DiscreteAxis& LogW) {
    double f_;
    if ( 2*h <= x && x <(int)LogW.size()-2*h){
      f_ = (-1.0 * LogW[x + 2 * h] + 16.0 * LogW[x + h] - 30.0 * LogW[x] +
            16.0 * LogW[x - h] - LogW[x - 2 * h]) /
           (12.0 * h * h * LogW.get_step() * LogW.get_step());
    }else if (2*h > x){
      f_ = df2(h, 2 * h, LogW);
    }else{
      f_ = df2(h, (int)LogW.size() - 2 * h - 1, LogW);
    }
    return f_;
  };

  std::vector<double> S_dE(logW.size());
  std::vector<double> S_dE2(logW.size());
  for (int i = 0; i < logW.size(); ++i){
    S_dE[i] = df(dh1, i, logW);
    S_dE2[i] = df2(dh2, i, logW);
  }

  std::vector<DiscreteAxis> jk_logW(
      NumBins, DiscreteAxis(minE, maxE, step, log_zero<double>()));

  for (int k = 0; k < NumBins; ++k){
    DiscreteAxis2D jk_DOS;
    dlib::deserialize(path_jk + "/" + std::to_string(k) + "/logDOS.obj") >> jk_DOS;
    DiscreteAxis norm(minE, maxE, step, 0.0);
    for (int i = 0; i < jk_DOS.get_num_bins_first(); ++i) {
      const double E1 = jk_DOS.get_value_first(i);
      for (int j = 0; j < jk_DOS.get_num_bins_second(); ++j) {
        const double E2 = jk_DOS.get_value_second(j);
        const double E = NVT::ham(E1, Kappa, E2);
        const int index = jk_DOS.get_index(i, j);
        const int rindex = jk_logW[k].get_bin(E);
        if (std::isfinite(jk_DOS[index] and
                          jk_DOS[index] > log_zero<double>() + 1.0e10)) {
          jk_logW[k][rindex] = addlogwise(jk_logW[k][rindex], jk_DOS[index]);
          norm[rindex] += 1;
        }
      }
    }

    for (int i = 0; i < jk_logW[k].size(); ++i) {
      if (norm[i] != 0.0 and jk_logW[k][i] > -1.0e100) {
        jk_logW[k][i] -= std::log(norm[i]);
      } else {
        jk_logW[k][i] = std::numeric_limits<double>::quiet_NaN();
      }
    }
  }
  
  std::vector<double> logW_mean(logW.size(), 0.0);
  std::vector<double> logW_err(logW.size(), 0.0);
  std::vector<double> S_mean(logW.size(), 0.0);
  std::vector<double> S_err(logW.size(), 0.0);
  std::vector<double> S_dE_mean(logW.size(), 0.0);
  std::vector<double> S_dE_err(logW.size(), 0.0);
  std::vector<double> S_dE2_mean(logW.size(), 0.0);
  std::vector<double> S_dE2_err(logW.size(), 0.0);

  for (size_t i = 0; i < jk_logW.size(); ++i) {
    for (int j = 0; j < jk_logW[i].size(); ++j) {
      logW_mean[j] += jk_logW[i][j];
      logW_err[j] += jk_logW[i][j] * jk_logW[i][j];
      double tmp = -1.0 * jk_logW[i][j];
      S_mean[j] += tmp;
      S_err[j] += tmp*tmp;
      tmp = df(dh1, j, jk_logW[i]);
      S_dE_mean[j] += tmp;
      S_dE_err[j] += tmp*tmp;
      tmp = df2(dh2, j, jk_logW[i]);
      S_dE2_mean[j] += tmp;
      S_dE2_err[j] += tmp*tmp;
    }
  }
  for (size_t i = 0; i < logW_mean.size(); ++i) {
    logW_mean[i] = logW_mean[i] / static_cast<double>(NumBins);
    logW_err[i] -= static_cast<double>(NumBins) * logW_mean[i] * logW_mean[i];
    logW_err[i] *=
        static_cast<double>(NumBins - 1) / static_cast<double>(NumBins);
    logW_err[i] = std::sqrt(logW_err[i]);
    S_mean[i] = S_mean[i] / static_cast<double>(NumBins);
    S_err[i] -= static_cast<double>(NumBins) * S_mean[i] * S_mean[i];
    S_err[i] *=
        static_cast<double>(NumBins - 1) / static_cast<double>(NumBins);
    S_err[i] = std::sqrt(S_err[i]);
    S_dE_mean[i] = S_dE_mean[i] / static_cast<double>(NumBins);
    S_dE_err[i] -= static_cast<double>(NumBins) * S_dE_mean[i] * S_dE_mean[i];
    S_dE_err[i] *=
        static_cast<double>(NumBins - 1) / static_cast<double>(NumBins);
    S_dE_err[i] = std::sqrt(S_dE_err[i]);
    S_dE2_mean[i] = S_dE2_mean[i] / static_cast<double>(NumBins);
    S_dE2_err[i] -= static_cast<double>(NumBins) * S_dE2_mean[i] * S_dE2_mean[i];
    S_dE2_err[i] *=
        static_cast<double>(NumBins - 1) / static_cast<double>(NumBins);
    S_dE2_err[i] = std::sqrt(S_dE2_err[i]);
  }

  std::string pathr = "results/logWeights";
  boost::filesystem::create_directories(pathr);

  std::string path_kappa = pathr + "/" + std::to_string(Kappa);
  boost::filesystem::create_directories(path_kappa);
  
  std::ofstream file(path_kappa + "/logWeight.dat", std::ios::trunc);
  file << "#E lnW(E) err\n";
  std::ofstream fileS(path_kappa + "/logS.dat", std::ios::trunc);
  fileS << "#E lnS(E) err\n";
  std::ofstream filedE(path_kappa + "/logS_dE.dat", std::ios::trunc);
  filedE << "#E lnS_dE(E) err\n";
  std::ofstream filedE2(path_kappa + "/logS_dE2.dat", std::ios::trunc);
  filedE2 << "#E lnS_dT2(E) err\n";
  for (int i = 0; i < logW.size(); ++i){
    double E = logW.get_value(i);
    file << E << " " << logW[i] << " " << logW_err[i] << "\n";
    fileS << E << " " << -1.0 * logW[i] << " " << S_err[i] << "\n";
    filedE << E << " " << S_dE[i] << " " << S_dE_err[i] << "\n";
    filedE2 << E << " " << S_dE2[i] << " " << S_dE2_err[i] << "\n";
  }


}
