#include <chrono>
#include <iostream>
#include <string>

#include "../../src/DiscreteAxis.hpp"
#include "../../src/WHAM/WHAM.hpp"
#include "../../src/data/all_matching.hpp"
#include "../../src/data/load_file.hpp"
#include "../../src/eigenvalues/add_eigenvalues.hpp"
#include "../../src/vector/helper.hpp"

#include "boost/filesystem.hpp"
#include "boost/progress.hpp"

using namespace consus;
using namespace consus::WHAM;

int main(int argc, char const *argv[])
{
  if (argc < 3){
    std::cerr << "arguments are missing\n";
    std::exit(1);
  }

  constexpr int col_energy = 2;
  const double num_bins = std::stod(argv[1]);
  assert(num_bins > 9);

  std::vector<std::string> filenames;
  for (int i = 2; i < argc; ++i){
    filenames.push_back(argv[i]);
  }

  auto get_temp = [] (const std::string& name) -> double {
    auto pos1 = name.find("PT_T");
    pos1 += 4;
    auto pos2 = name.find(".dat", pos1);
    if (std::string::npos != pos1 && std::string::npos != pos2) {
      return std::stod(name.substr(pos1, pos2 - pos1));
    } else {
      std::cerr << "Could not determine temperature from " << name << "\n";
      std::exit(1);
    }

    return 0.0;
  };

  std::sort(filenames.begin(), filenames.end(), [&](const std::string & a,
                                                    const std::string & b) {
    double a_t = get_temp(a);
    double b_t = get_temp(b);
    return a_t < b_t;
  });


  double minE = std::numeric_limits<double>::max();
  double maxE = std::numeric_limits<double>::lowest();

  std::cout << "searchin min/max values in data\n";
  boost::progress_display progress(filenames.size());
  #pragma omp parallel for reduction(min:minE) reduction(max:maxE)
  for (size_t i = 0; i < filenames.size(); ++i) {
    #pragma omp critical
    ++progress;
    auto timeseries = read_ssv(filenames[i]);
    auto resE =
        std::minmax_element(timeseries[col_energy].begin(), timeseries[col_energy].end());
    minE = std::min(minE, *(resE.first));
    maxE = std::max(maxE, *(resE.second));
  }
  const double step = (maxE - minE)/num_bins;

  std::vector<double> Parameters;
  std::vector<double> length_Hists;
  Range range= {minE, step, maxE};
  
  std::string path("analysis");
  std::string pathh("analysis/Hists");
  boost::filesystem::create_directories(path);
  boost::filesystem::create_directories(pathh);

  for (size_t i = 0; i < filenames.size(); ++i) {
    Parameters.push_back(1.0 / get_temp(filenames[i]));
  }
  
  std::cout << "loading files\n";
  boost::progress_display progress2(filenames.size());
  #pragma omp parallel for
  for (size_t i = 0; i < filenames.size(); ++i) {
      DiscreteAxis Hist;
      double length_hist;
      #pragma omp critical
      ++progress2;
      auto timeseries = read_ssv(filenames[i]);
      std::tie(Hist, length_hist) = make_histogram(timeseries[col_energy], range);
      dlib::serialize(pathh + "/Hist_" + std::to_string(i) + ".obj") << Hist;
      dlib::serialize(pathh + "/length_hist_" + std::to_string(i) + ".obj") << length_hist;
  }

  assert(Parameters.size() == filenames.size());

  std::cout << "Range " << range << "\n";

  for (size_t i = 0; i < Parameters.size(); ++i){
    dlib::serialize(pathh + "/Parameter_" + std::to_string(i) + ".obj")
        << Parameters[i];
  }

  dlib::serialize(path + "/range.obj") << range;
  dlib::serialize(path + "/filenames.obj") << filenames;
  dlib::serialize(path + "/Parameters.obj") << Parameters;
}
