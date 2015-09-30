#include <string>
#include <iostream>
#include <fstream>
#include <chrono>
#include <regex>
#include "../../src/vector/helper.hpp"
#include "../../src/DiscreteAxis.hpp"
#include "../../src/WHAM2D/WHAM2D.hpp"

#include "boost/filesystem.hpp"

using namespace consus;
using namespace consus::WHAM2D;
typedef std::chrono::high_resolution_clock myclock;

int main()
{
  std::string path("analysis");
  
  DiscreteAxis2D DOS;
  std::vector<DiscreteAxis2D> MicroMeans;
  std::vector<std::string> header;

  dlib::deserialize(path + "/logDOS.obj") >> DOS;
  dlib::deserialize(path + "/header.obj") >> header;
  dlib::deserialize(path + "/MicroMeans.obj") >> MicroMeans;

  path = "results/logWeights";

  boost::filesystem::create_directories(path);

  std::cout << "size of logDOS " << DOS.size() << "\n";
  std::cout << DOS.get_num_bins_first() << " " << DOS.get_num_bins_second()
            << "\n";

  std::ofstream file(path + "/logWeight.dat", std::ios::trunc);
  file << "#E1 E2 logW(E1, E2)\n";

  for (int i = 0; i < DOS.get_num_bins_first(); ++i){
    double E1 = DOS.get_value_first(i);
    for (int j = 0; j < DOS.get_num_bins_second(); ++j){
      double E2 = DOS.get_value_second(j);
      file << E1 << " " << E2 << " " << DOS[DOS.get_index(i, j)] << "\n";
    }
    file << "\n";
  }
}
