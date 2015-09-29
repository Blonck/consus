#pragma once

#include <array>
#include <string>
#include <vector>
#include <cmath>
#include <tuple>
#include <algorithm>
#include "../vector/helper.hpp"

namespace consus
{

namespace eig {

/// calculates the determinant of a symmetric  3x3 matrix
inline double det_sym(const double Rxx, const double Ryy, const double Rzz,
                      const double Rxy, const double Rxz, const double Ryz) {
  return Rxx * Ryy * Rzz + 2.0 * (Rxy * Ryz * Rxz) - Rxz * Ryy * Rxz -
         Rxy * Rxy * Rzz - Rxx * Ryz * Ryz;
}

/// sort the 3 values in decreasing order
void sort_eig(double& a, double& b, double& c) {
  std::array<double, 3> tmp{ {a, b, c} };
  std::sort(tmp.begin(), tmp.end());
  a = tmp[2]; b = tmp[1]; c = tmp[0];
}

/// replace the entries of a header with ones corresponding to the
/// eigenvalues of the gyration tensor and invariants of iit
void replace_eigenvalues_header(std::vector<std::string>& header,
                                const size_t Rxx_index, const size_t Ryy_index,
                                const size_t Rzz_index, const size_t Rxy_index,
                                const size_t Rxz_index,
                                const size_t Ryz_index) {
  header[Rxx_index] = "l1";
  header[Ryy_index] = "l2";
  header[Rzz_index] = "l3";
  header[Rxy_index] = "b";
  header[Rxz_index] = "c";
  header[Ryz_index] = "k2";
}

/// calculate the 3 eigenvalues and 3 shape descriptors of the gyration tensor
/// 
/// http://en.wikipedia.org/w/index.php?title=Eigenvalue_algorithm&oldid=590404829
inline std::array<double, 6> eigenvalues_plus_impl(
    const double Rxx, const double Ryy, const double Rzz, const double Rxy,
    const double Rxz, const double Ryz) {

  std::array<double, 6> res;
  double l1, l2, l3;
  double p1 = Rxy * Rxy + Rxz * Rxz + Ryz * Ryz;
  double phi;
  if (p1 == 0.0) {
    l1 = Rxx;
    l2 = Ryy;
    l3 = Rzz;
    sort_eig(l1, l2, l3);
  } else {
    double q = (Rxx + Ryy + Rzz) / 3.0;
    double p2 = (Rxx - q) * (Rxx - q) + (Ryy - q) * (Ryy - q) +
                (Rzz - q) * (Rzz - q) + 2.0 * p1;
    double p = std::sqrt(p2 / 6.0);
    double tq = 1.0 / p;
    double r = 0.5 * det_sym(tq * (Rxx - q), tq * (Ryy - q), tq * (Rzz - q),
                             tq * Rxy, tq * Rxz, tq * Ryz);
    if (r <= -1.0) {
      phi = M_PI / 3.0;
    } else {
      if (r >= 1) {
        phi = 0;
      } else {
        phi = std::acos(r) / 3.0;
      }
    }
    l1 = q + 2.0 * p * std::cos(phi);
    l3 = q + 2.0 * p * std::cos(phi + (2.0 * M_PI / 3.0));
    l2 = 3.0 * q - l1 - l3;
  }
#ifndef NDEBUG
  if (l1 < l2 or l1 < l3 or l2 < l3) {
    std::cerr << "l1 " << l1 << " l2 " << l2 << " l3 " << l3 << "\n";
    std::abort();
  }
#endif
  res[0] = l1;
  res[1] = l2;
  res[2] = l3;
  res[3] = l1 - 0.5 * (l2 + l3);
  res[4] = l2 - l3;
  res[5] = 1.5 * (std::pow(l1, 2) + std::pow(l2, 2) + std::pow(l3, 2)) /
               (std::pow(l1 + l2 + l3, 2)) - 0.5;
  if (res[5] > 1.0 or res[5] < 0){
    std::cerr << "Error on calculating relative shape anisotropy\n";
    std::terminate();
  }
#ifndef NDEBUG
  if (res[5] > 1.0 or res[5] < 0){
    std::cerr << Rxx << " " << Ryy << " " << Rzz << "\n";
    std::cerr << Rxy << " " << Rxz << " " << Ryz << "\n";
    std::cerr << phi << "\n";
    std::cerr << res[5] << " " << l1 << " " << l2 << " " << l3 << "\n";
    std::terminate();
  }
#endif

  return res;
}

/// replace entries of the gyration tensor with the eigenvalues and invariants
/// of it
///
/// \param data - 2 dim vector every 1 dim vector contains a time series
/// the entries of the elements of the gyration tensor will be changed to
/// the 3 eigenvalues and b, c and k2
/// \param Rxx_index - index of the first diagonal element of gyration tensor
/// \param Ryy_index - index of the second diagonal element of gyration tensor
/// \param Rzz_index - index of the third diagonal element of gyration tensor
/// \param Rxy_index - index of the (1,2) and (2,1) element of gyration tensor
/// \param Rxz_index - index of the (1,3) and (3,1) element of gyration tensor
/// \param Ryz_index - index of the (2,3) ans (3,2) element of gyration tensor
template <class T = double>
void replace_eigenvalues(vec2<T>& data, const size_t Rxx_index,
                         const size_t Ryy_index, const size_t Rzz_index,
                         const size_t Rxy_index, const size_t Rxz_index,
                         const size_t Ryz_index) {
  for (size_t i = 0; i < data[0].size(); ++i) {
    auto tmp = eigenvalues_plus_impl(data[Rxx_index][i], data[Ryy_index][i],
                                     data[Rzz_index][i], data[Rxy_index][i],
                                     data[Rxz_index][i], data[Ryz_index][i]);
    data[Rxx_index][i] = tmp[0];
    data[Ryy_index][i] = tmp[1];
    data[Rzz_index][i] = tmp[2];
    data[Rxy_index][i] = tmp[3];
    data[Rxz_index][i] = tmp[4];
    data[Ryz_index][i] = tmp[5];
  }
}

template <class T = double>
void add_eigenvalues(vec2<T>& data, const size_t Rxx_index,
                     const size_t Ryy_index, const size_t Rzz_index,
                     const size_t Rxy_index, const size_t Rxz_index,
                     const size_t Ryz_index) {
  size_t last_index = data.size();
  for (size_t i = 0; i < 6; ++i){
    data.push_back( vec1<T>(data[0].size()) );
  }
  for (size_t i = 0; i < data[0].size(); ++i) {
    auto tmp = eigenvalues_plus_impl(data[Rxx_index][i], data[Ryy_index][i],
                                     data[Rzz_index][i], data[Rxy_index][i],
                                     data[Rxz_index][i], data[Ryz_index][i]);
    data[last_index][i] = tmp[0];
    data[last_index + 1][i] = tmp[1];
    data[last_index + 2][i] = tmp[2];
    data[last_index + 3][i] = tmp[3];
    data[last_index + 4][i] = tmp[4];
    data[last_index + 5][i] = tmp[5];
  }
}

} /* end of namespace eig */

template <class T = double>
void replace_gyr_tensor(vec2<T>& data, std::vector<std::string>& header,
                        const std::vector<size_t>& col_gyr) {

  if (col_gyr.size() != 6) {
    std::cerr << "col_gyr should get 6 values"
              << "\n";
    std::terminate();
  }
  auto less_then = [&](const size_t& col) { return col > (header.size() - 1); };
  if (std::find_if(col_gyr.begin(), col_gyr.end(), less_then) !=
      col_gyr.end()) {
    std::cerr << "col_gyr should not get values greater then the number of "
                 "time series" << header.size() << "\n";
  }
  if (not((header[col_gyr[0]] == "Rxx")and(header[col_gyr[1]] == "Ryy")
          and(header[col_gyr[2]] == "Rzz") and(header[col_gyr[3]] == "Rxy")
          and(header[col_gyr[4]] == "Rxz") and(header[col_gyr[5]] == "Ryz"))) {
    std::cerr << "header of elements of Gyration Tensor does not fit\n";
    std::exit(1);
  }
  std::cout << "add eigenvalues to time series\n";
  eig::replace_eigenvalues_header(header, col_gyr[0], col_gyr[1], col_gyr[2],
                                  col_gyr[3], col_gyr[4], col_gyr[5]);
  eig::replace_eigenvalues(data, col_gyr[0], col_gyr[1], col_gyr[2], col_gyr[3],
                           col_gyr[4], col_gyr[5]);
}

template <class T = double>
void add_gyr_tensor(vec2<T>& data, std::vector<std::string>& header,
                    const std::vector<size_t>& col_gyr) {
  if (col_gyr.size() != 6) {
    std::cerr << "col_gyr should have 6 values"
              << "\n";
    std::exit(1);
  }
  auto less_then = [&](const size_t& col) { return col > (header.size() - 1); };
  if (std::find_if(col_gyr.begin(), col_gyr.end(), less_then) !=
      col_gyr.end()) {
    std::cerr << "col_gyr should not get values greater then the number of "
                 "time series" << header.size() << "\n";
  }
  if (not((header[col_gyr[0]] == "Rxx")and(header[col_gyr[1]] == "Ryy")
          and(header[col_gyr[2]] == "Rzz") and(header[col_gyr[3]] == "Rxy")
          and(header[col_gyr[4]] == "Rxz") and(header[col_gyr[5]] == "Ryz"))) {
    std::cerr << "header of elements of Gyration Tensor does not fit\n";
    std::exit(1);
  }
  std::cout << "add eigenvalues to time series\n";
  eig::replace_eigenvalues_header(header, col_gyr[0], col_gyr[1], col_gyr[2],
                                  col_gyr[3], col_gyr[4], col_gyr[5]);
  eig::replace_eigenvalues(data, col_gyr[0], col_gyr[1], col_gyr[2], col_gyr[3],
                           col_gyr[4], col_gyr[5]);
}

} /* end of namespace consus */
