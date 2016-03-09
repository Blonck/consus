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
  
  std::vector<std::string> header;
  std::vector<std::string> filenames;
  std::vector<std::pair<double, double>> Parameters;

  dlib::deserialize(path + "/header.obj") >> header;
  dlib::deserialize(path + "/filenames.obj") >> filenames;
  dlib::deserialize(path + "/Parameters.obj") >> filenames;

  std::string pathr = "results/canonical";
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
  if (filenames.size() != Parameters.size()){
    std::cerr << "Size of filenames and Parameters differs\n";
    std::exit(1);
  }

  vec2<double> results(

  for (size_t i = 0; i < filenames.size(); ++i){
    auto timeseries = read_ssv(filenames[i]);
  }

  std::string path_files = path_kappa + "/" + std::to_string(Kappa);
  boost::filesystem::create_directories(path_files);
  std::string filename = path_files + "/" + replacediv(header[column]) + ".dat";
  std::ofstream out(filename, std::ios::trunc);
  out << "#Temperature " << header[column] << " "
      << "error(" << header[column] << ")\n";
  for (size_t j = 0; j < Temperatures.size(); ++j) {
    out << Temperatures[j] << " " << results[j] << " " << jk_error[j] << "\n";
  }
  out.close();
}
