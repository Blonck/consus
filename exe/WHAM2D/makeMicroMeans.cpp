#define VERBOSE 

#include <string>
#include <iostream>
#include <chrono>
#include <regex>
#include "../../src/data/load_file.hpp"
#include "../../src/data/all_matching.hpp"
#include "../../src/vector/helper.hpp"
#include "../../src/DiscreteAxis2D.hpp"
#include "../../src/WHAM2D/WHAM2D.hpp"
#include "../../src/eigenvalues/add_eigenvalues.hpp"

#include "boost/filesystem.hpp"

using namespace consus;
using namespace consus::WHAM2D;

typedef std::chrono::high_resolution_clock myclock;

int main(int argc, char const *argv[])
{
  if (argc < 2){
    std::cerr << "arguments are missing\n";
    std::exit(1);
  }
  int NumBins = std::stod(argv[1]);
  const int col1 = 0;
  const int col2 = 1;

  std::vector<std::string> filenames;
  std::vector<std::string> header;
  Range rangeE1;
  Range rangeE2;
  DiscreteAxis2D Hist;
  std::string path("analysis");
  std::string path_jk("analysis/JK");
  dlib::deserialize(path + "/Hist.obj") >> Hist;
  dlib::deserialize(path + "/rangeE1.obj") >> rangeE1;
  dlib::deserialize(path + "/rangeE2.obj") >> rangeE2;
  dlib::deserialize(path + "/header.obj") >> header;
  dlib::deserialize(path + "/filenames.obj") >> filenames;

  bool add_eig = false;
  if (header.size() > 10) {
    if (header[5] == "Rxx" and header[6] == "Ryy" and header[7] ==
        "Rzz" and header[8] == "Rxy" and header[9] == "Rxz" and header[10] ==
        "Ryz") {
      add_eig = true;
    }
  }
  if (add_eig){
    consus::eig::add_eigenvalues_header(header);
  }
  std::cout << header << "\n";

  vec1<DiscreteAxis2D> MicroMeans(header.size(), DiscreteAxis2D(rangeE1, rangeE2));
  vec2<DiscreteAxis2D> jk_MicroMeans(
      NumBins,
      vec1<DiscreteAxis2D>(header.size(), DiscreteAxis2D(rangeE1, rangeE2)));
  for (size_t i = 0; i < filenames.size(); ++i) {
    std::cout << filenames[i] << "\n";
    auto timeseries = read_ssv(filenames[i]);
    int length_bins = timeseries[0].size() / NumBins;
    if (timeseries[0].size() % NumBins != 0) {
      std::cerr << "bin size and time series size doesn't fit\n";
      std::exit(1);
    }
    if (add_eig){
      consus::eig::add_eigenvalues(timeseries, 5, 6, 7, 8, 9, 10);
    }
    for (size_t j = 0; j < timeseries.size(); ++j) {
      const double E1 = timeseries[col1][i];
      const double E2 = timeseries[col2][i];
      int index = MicroMeans[0].get_bin(E1, E2);
      for (size_t k = 0; k < MicroMeans.size(); ++k) {
        MicroMeans[k][index] += timeseries[k][i];
      }
    }
    std::cout << "collecting jackknife blocks\n";
    for (int j = 0; j < NumBins; ++j) {
      std::cout << j << " ";
      for (int k = 0; k < j; ++k) {
        for (int l = 0; l < length_bins; ++l) {
          int i_ts = k * length_bins + l;
          const double E1 = timeseries[col1][i_ts];
          const double E2 = timeseries[col2][i_ts];
          int index = MicroMeans[0].get_bin(E1, E2);
          for (size_t h = 0; h < header.size(); ++h) {
            jk_MicroMeans[j][h][index] = timeseries[h][i_ts];
          }
        }
      }
      for (int k = j + 1; k < NumBins; ++k) {
        for (int l = 0; l < length_bins; ++l) {
          int i_ts = k * length_bins + l;
          const double E1 = timeseries[col1][i_ts];
          const double E2 = timeseries[col2][i_ts];
          int index = MicroMeans[0].get_bin(E1, E2);
          for (size_t h = 0; h < header.size(); ++h) {
            jk_MicroMeans[j][h][index] = timeseries[h][i_ts];
          }
        }
      }
    }
    std::cout << "\n";
  }

  for (size_t i = 0; i < MicroMeans.size(); ++i){
    for (int j = 0; j < MicroMeans[i].size(); ++j){
      MicroMeans[i][j] /= std::exp(Hist[j]);
    }
  }
  dlib::serialize(path + "/MicroMeans.obj") << MicroMeans;
  for (int i = 0; i < NumBins; ++i) {
    DiscreteAxis2D jk_Hist;
    dlib::deserialize(path_jk + "/" + std::to_string(i)) >> jk_Hist;
    for (size_t j = 0; j < jk_MicroMeans[i].size(); ++j) {
      for (int k = 0; k < jk_MicroMeans[i][j].size(); ++k) {
        jk_MicroMeans[i][j][k] /= std::exp(jk_Hist[k]);
      }
    }
    dlib::serialize(path_jk + "/" + std::to_string(i) + "/MicroMeans.obj")
        << jk_MicroMeans[i];
  }

  dlib::serialize(path + "/header.obj") << header;
  dlib::serialize(path + "/add_eigenvalues.obj") << add_eig;
}
