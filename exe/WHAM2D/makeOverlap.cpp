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

int main()
{

  std::vector<std::string> filenames;
  
  std::string path("analysis");
  std::string pathh("analysis/Hists");

  dlib::deserialize(path + "/filenames.obj") >> filenames;

  size_t NH = filenames.size();
  std::cout << "calculating overlap\n";
  std::cout << NH << "\n";

  #pragma omp parallel for
  for (size_t i = 0; i < NH; ++i){
    std::cout << "Overlap for " << i << "th histogram\n";
    vec1<int> Overlap(i);
    DiscreteAxis2D H1;
    dlib::deserialize(pathh + "/Hist_"  + std::to_string(i) + ".obj") >> H1;
    for (size_t j = 0; j < i; ++j){
      DiscreteAxis2D H2;
      dlib::deserialize(pathh + "/Hist_" + std::to_string(j) + ".obj") >> H2;
      for (int k = 0; k < H1.size(); ++k){
        if (std::exp(H1[k]) > 1.1 and
            std::exp(H2[k]) > 1.1) {
          ++Overlap[j];
        }
      }
    }
    assert(Overlap.size() == i);
    dlib::serialize(pathh + "/Overlap_" + std::to_string(i) + ".obj") << Overlap;
  }
}
