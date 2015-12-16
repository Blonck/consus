#include <string>
#include <iostream>
#include <algorithm>
#include <chrono>
#include "../../src/data/load_file.hpp"
#include "../../src/data/all_matching.hpp"
#include "../../src/vector/helper.hpp"
#include "../../src/DiscreteAxis2D.hpp"
#include "../../src/WHAM2D/WHAM2D.hpp"
#include "../../src/simple.hpp"
#include "../../src/eigenvalues/add_eigenvalues.hpp"

#include "boost/filesystem.hpp"
#include "boost/regex.hpp"
#include "boost/progress.hpp"

using namespace consus;
using namespace consus::WHAM2D;
typedef std::chrono::high_resolution_clock myclock;

int main()
{
  std::string path("analysis");
  
  std::vector<std::string> header;
  std::vector<std::string> filenames;
  std::vector<std::pair<double, double>> Parameters;

  dlib::deserialize(path + "/header.obj") >> header;
  dlib::deserialize(path + "/filenames.obj") >> filenames;
  dlib::deserialize(path + "/Parameters.obj") >> Parameters;

  bool add_eig = false;
  int first_gyr = 0;
  if (header.size() > 10) {
    auto cur =  std::find(header.begin(), header.end(), "Rxx");
    first_gyr = std::distance(header.begin(), cur);
    if (cur != header.end() and (cur+5) < header.end()){
      add_eig = true;
    }
  }

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
  if (filenames.size() != Parameters.size()){
    std::cerr << "Size of filenames and Parameters differs\n";
    std::exit(1);
  }

  vec1<double> Kappas;
  vec1<double> Temps;
  for (auto& a: Parameters){
    a.first = 1.0/a.first;
    Temps.push_back(a.first);
    Kappas.push_back(a.second);
  }


  std::sort(Temps.begin(), Temps.end());
  Temps.erase(std::unique(Temps.begin(), Temps.end()), Temps.end());
  std::sort(Kappas.begin(), Kappas.end());
  Kappas.erase(std::unique(Kappas.begin(), Kappas.end()), Kappas.end());

  std::for_each(Temps.begin(), Temps.end(), [](double e){ std::cout << e << " ";});
  std::cout << "\n";
  std::for_each(Kappas.begin(), Kappas.end(), [](double e){ std::cout << e << " ";});
  std::cout << "\n";

  vec2<double> results(Parameters.size(), vec1<double>(header.size(), 0));
  vec2<double> errors(Parameters.size(), vec1<double>(header.size(), 0));

  boost::progress_display progress(filenames.size() * header.size());
  for (size_t i = 0; i < filenames.size(); ++i){
    auto timeseries = read_ssv(filenames[i]);
    if (add_eig){
      consus::eig::add_eigenvalues(timeseries, first_gyr, first_gyr + 1,
                                   first_gyr + 2, first_gyr + 3, first_gyr + 4,
                                   first_gyr + 5);
    }
    for (size_t j = 0; j < header.size(); j++){
      auto res = simple_stats(timeseries[j], 100);
      results[i][j] = std::get<2>(res);
      errors[i][j] = std::get<5>(res);
      ++progress;
    }
  }

  std::string pathr = "results/canonical";
  boost::filesystem::create_directories(pathr);

  std::string path_kappa = pathr + "/Kappa";
  boost::filesystem::create_directories(path_kappa);
  std::string path_temp = pathr + "/Temperatures";
  boost::filesystem::create_directories(path_kappa);

  for (auto& Kappa : Kappas) {
    std::string path_files = path_kappa + "/" + std::to_string(Kappa);
    boost::filesystem::create_directories(path_files);
    for (size_t column = 0; column < header.size(); ++column) {
      std::string filename =
          path_files + "/" + replacediv(header[column]) + ".dat";
      std::ofstream out(filename, std::ios::trunc);
      out << "#Temperature " << header[column] << " "
          << "error(" << header[column] << ")\n";
      for (size_t j = 0; j < Temps.size(); ++j) {
        for (size_t i = 0; i < Parameters.size(); ++i) {
          if (Temps[j] == std::get<0>(Parameters[i]) and
              Kappa == std::get<1>(Parameters[i])) {
            out << Temps[j] << " " << results[i][column] << " "
                << errors[i][column] << "\n";
          }
        }
      }
      out.close();
    }
  }
}
