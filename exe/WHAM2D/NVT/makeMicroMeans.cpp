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
#include "boost/progress.hpp"

using namespace consus;
using namespace consus::WHAM2D;

int main()
{
  int NumBins;
  const int col1 = 0;
  const int col2 = 1;

  std::vector<std::string> filenames;
  Range rangeE1;
  Range rangeE2;
  DiscreteAxis2D Hist;
  std::string path("analysis");
  std::string path_jk("analysis/JK");
  dlib::deserialize(path + "/Hist.obj") >> Hist;
  dlib::deserialize(path + "/rangeE1.obj") >> rangeE1;
  dlib::deserialize(path + "/rangeE2.obj") >> rangeE2;
  dlib::deserialize(path + "/filenames.obj") >> filenames;
  dlib::deserialize(path + "/NumBins.obj") >> NumBins;

  std::vector<std::string> header = read_header(filenames[0]);
  std::cout << header << "\n";

  bool add_eig = false;
  int first_gyr = 0;
  if (header.size() > 10) {
    auto cur =  std::find(header.begin(), header.end(), "Rxx");
    first_gyr = std::distance(header.begin(), cur);
    if (cur != header.end() and (cur+5) < header.end()){
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
  
  boost::progress_display progress(filenames.size());
  for (size_t i = 0; i < filenames.size(); ++i) {
    ++progress;
    auto timeseries = read_ssv(filenames[i]);
    int length_bins = timeseries[0].size() / NumBins;
    if (timeseries[0].size() % NumBins != 0) {
      std::cerr << "bin size and time series size doesn't fit\n";
      std::exit(1);
    }
    if (add_eig){
      consus::eig::add_eigenvalues(timeseries, first_gyr, first_gyr + 1,
                                   first_gyr + 2, first_gyr + 3, first_gyr + 4,
                                   first_gyr + 5);
    }
    for (size_t j = 0; j < timeseries[col1].size(); ++j) {
      const double E1 = timeseries[col1][j];
      const double E2 = timeseries[col2][j];
      int index = MicroMeans[0].get_bin(E1, E2);
      for (size_t k = 0; k < MicroMeans.size(); ++k) {
        MicroMeans[k][index] += timeseries[k][j];
      }
    }
    #pragma omp parallel for
    for (int j = 0; j < NumBins; ++j) {
      for (int k = 0; k < j; ++k) {
        for (int l = 0; l < length_bins; ++l) {
          int i_ts = k * length_bins + l;
          const double E1 = timeseries[col1][i_ts];
          const double E2 = timeseries[col2][i_ts];
          int index = MicroMeans[0].get_bin(E1, E2);
          for (size_t h = 0; h < header.size(); ++h) {
            jk_MicroMeans[j][h][index] += timeseries[h][i_ts];
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
            jk_MicroMeans[j][h][index] += timeseries[h][i_ts];
          }
        }
      }
    }
  }

  std::cout << "normalize\n";
  for (size_t i = 0; i < MicroMeans.size(); ++i){
    for (int j = 0; j < MicroMeans[i].size(); ++j){
      if (Hist[j] != 0.0){
        MicroMeans[i][j] /= Hist[j];
      }
    }
  }
  dlib::serialize(path + "/MicroMeans.obj") << MicroMeans;
  for (int i = 0; i < NumBins; ++i) {
    DiscreteAxis2D jk_Hist;
    dlib::deserialize(path_jk + "/" + std::to_string(i) + "/Hist.obj") >> jk_Hist;
    for (size_t j = 0; j < jk_MicroMeans[i].size(); ++j) {
      for (int k = 0; k < jk_MicroMeans[i][j].size(); ++k) {
        if (Hist[k] != 0.0) {
          jk_MicroMeans[i][j][k] /= jk_Hist[k];
        }
      }
    }
    dlib::serialize(path_jk + "/" + std::to_string(i) + "/MicroMeans.obj")
        << jk_MicroMeans[i];
  }

  dlib::serialize(path + "/header.obj") << header;
  dlib::serialize(path + "/add_eigenvalues.obj") << add_eig;
}
