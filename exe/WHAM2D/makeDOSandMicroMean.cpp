#include <string>
#include <iostream>
#include <chrono>
#include <regex>
#include "../../src/load_file.hpp"
#include "../../src/vector/helper.hpp"
#include "../../src/DiscreteAxis.hpp"
#include "../../src/WHAM2d/make_histogram2d.hpp"
#include "../../src/eigenvalues/add_eigenvalues.hpp"

#include "boost/filesystem.hpp"

using namespace consus;
using namespace consus::WHAM2D;

typedef std::chrono::high_resolution_clock myclock;

int main(int argc, char const *argv[])
{
  if (argc < 5){
    std::cerr << "arguments are missing\n";
    std::exit(1);
  }

  constexpr int rowE1 = 0;
  constexpr int rowE2 = 1;
  const double stepE1 = std::stod(argv[1]);
  const double stepE2 = std::stod(argv[2]);
  const double devmax = std::stod(argv[3]);

  std::vector<std::string> filenames;
  for (int i = 4; i < argc; ++i){
    filenames.push_back(argv[i]);
    std::cout << argv[i] << "\n";
  }
  auto header = read_header(filenames[0]);
  bool replace_eig = false;
  if (header.size() > 10) {
    if (header[5] == "Rxx" and header[6] == "Ryy" and header[7] ==
        "Rzz" and header[8] == "Rxy" and header[9] == "Rxz" and header[10] ==
        "Ryz") {
      replace_eig = true;
    }
  }
  if (replace_eig){
    consus::eig::replace_eigenvalues_header(header, 5, 6, 7, 8, 9, 10);
  }
  std::cout << header << "\n";

  double minE1 = std::numeric_limits<double>::max();
  double maxE1 = std::numeric_limits<double>::lowest();
  double minE2 = std::numeric_limits<double>::max();
  double maxE2 = std::numeric_limits<double>::lowest();

#pragma omp parallel for reduction(min:minE1, minE2) reduction(max:maxE1, maxE2)
  for (size_t i = 0; i < filenames.size(); ++i) {
    std::cout << filenames[i] << "\n";
    auto timeseries = read_ssv(filenames[i]);
    auto resE1 =
        std::minmax_element(timeseries[rowE1].begin(), timeseries[rowE1].end());
    auto resE2 =
        std::minmax_element(timeseries[rowE2].begin(), timeseries[rowE2].end());
    minE1 = std::min(minE1, *(resE1.first));
    maxE1 = std::max(maxE1, *(resE1.second));
    minE2 = std::min(minE2, *(resE2.first));
    maxE2 = std::max(maxE2, *(resE2.second));
  }

  vec1d Betas;
  vec1d Kappas;

  std::regex param_regex(".*timeseries/(.*)/PT_T(.*).dat", std::regex::egrep);
  std::smatch match;

  DiscreteAxis2D Hist;
  vec1d length_Hists;
  std::vector<DiscreteAxis2D> MicroMeans;
  DiscreteAxis2D DOS;
  Range rangeE1 = {minE1, stepE1, maxE1};
  Range rangeE2 = {minE2, stepE2, maxE2};
  for (size_t i = 0; i < filenames.size(); ++i) {
    std::cout << filenames[i] << "\n";
    std::regex_search(filenames[i], match, param_regex);
    if (match.size() != 3) {
      std::cerr << "could not determine parameters"
                << "\n";
      std::exit(1);
    }
    Betas.push_back(1.0 / std::stod(match[2]));
    Kappas.push_back(std::stod(match[1]));
    std::cout << "load file\n";
    auto timeseries = read_ssv(filenames[i]);
    if (replace_eig){
      consus::eig::replace_eigenvalues(timeseries, 5, 6, 7, 8, 9, 10);
    }
    std::cout << "make histogram" << "\n";
    if( i == 0){
      std::tie(Hist, length_Hists, MicroMeans) =
          make_histogram2d(timeseries, rowE1, rowE2, rangeE1, rangeE2);
      DOS = logDOS(Hist, length_Hists, Betas, Kappas, devmax);
    } else {
      add_histogram2d(timeseries, rowE1, rowE2, Hist, length_Hists,
                      MicroMeans);
      logDOS(DOS, Hist, length_Hists, Betas, Kappas, devmax);
    }
  }
  DOS = logDOS(Hist, length_Hists, Betas, Kappas, devmax);
  finalize(Hist, MicroMeans);

  std::string path("analysis");
  std::cout << "RangeE1 " << rangeE1 << "\n";
  std::cout << "RangeE2 " << rangeE2 << "\n";
  boost::filesystem::create_directories(path);

  dlib::serialize(path + "/rangeE1.obj") << rangeE1;
  dlib::serialize(path + "/rangeE2.obj") << rangeE2;
  dlib::serialize(path + "/filenames.obj") << filenames;
  dlib::serialize(path + "/Hist.obj") << Hist;
  dlib::serialize(path + "/logDOS.obj") << DOS;
  dlib::serialize(path + "/header.obj") << header;
  dlib::serialize(path + "/MicroMeans.obj") << MicroMeans;
  dlib::serialize(path + "/length_Hists.obg") << length_Hists;
  dlib::serialize(path + "/Betas.obj") << Betas;
  dlib::serialize(path + "/Kappas.obj") << Kappas;

}
