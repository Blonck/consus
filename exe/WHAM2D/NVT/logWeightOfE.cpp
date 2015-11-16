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

  const double Kappa = std::stod(argv[1]);
  const double step = std::stod(argv[2]);
  const int dh1 = std::stoi(argv[3]);
  const int dh2 = std::stoi(argv[4]);
  
  DiscreteAxis2D DOS;
  std::vector<double> length_Hists;
  std::vector<DiscreteAxis2D> MicroMeans;
  std::vector<std::string> header;

  dlib::deserialize(path + "/logDOS.obj") >> DOS;
  dlib::deserialize(path + "/header.obj") >> header;
  dlib::deserialize(path + "/MicroMeans.obj") >> MicroMeans;

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

  std::function<double(int, int)>  df = [&](int h, int x){
    double f_;
    if ( 2*h <= x && x <(int)logW.size()-2*h){
      f_ = (-1.0 * logW[x + 2 * h] + 8 * logW[x + h] - 8 * logW[x - h] +
            logW[x - 2 * h]) /
           (12.0 * h * logW.get_step());
    }else if (x < 2*h){
      f_ = df(h,2*h);
    }else{
      f_ = df(h,(int)logW.size() - 2*h-1);
    }
    return f_;
  };

  std::function<double(int, int)> df2 = [&](int h, int x){
    double f_;
    if ( 2*h <= x && x <(int)logW.size()-2*h){
      f_ = (-1.0 * logW[x + 2 * h] + 16.0 * logW[x + h] - 30.0 * logW[x] +
            16.0 * logW[x - h] - logW[x - 2 * h]) /
           (12.0 * h * h * logW.get_step() * logW.get_step());
    }else if (2*h > x){
      f_ = df2(h,2*h);
    }else{
      f_ = df2(h,(int)logW.size() - 2*h-1);
    }
    return f_;
  };


  std::string pathr = "results/logWeights";
  boost::filesystem::create_directories(pathr);

  std::string path_kappa = pathr + "/" + std::to_string(Kappa);
  boost::filesystem::create_directories(path_kappa);
  
  std::ofstream file(path_kappa + "/logWeight.dat", std::ios::trunc);
  file << "#E lnW(E)\n";
  std::ofstream fileS(path_kappa + "/logS.dat", std::ios::trunc);
  fileS << "#E lnS(E)\n";
  std::ofstream filedE(path_kappa + "/logS_dE.dat", std::ios::trunc);
  filedE << "#E lnS_dE(E)\n";
  std::ofstream filedE2(path_kappa + "/logS_dE2.dat", std::ios::trunc);
  filedE2 << "#E lnS_dT2(E)\n";
  for (int i = 0; i < logW.size(); ++i){
    double E = logW.get_value(i);
    file << E << " " << logW[i] << "\n";
    fileS << E << " " << -1.0 * logW[i] << "\n";
    filedE << E << " " << df(dh1, i) << "\n";
    filedE2 << E << " " << df2(dh2, i) << "\n";
  }


}
