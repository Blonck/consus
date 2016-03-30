#include <limits>
#include <algorithm>
#include <cmath>
#include <vector>
#include "vector/helper.hpp"

namespace consus{

// optimization routines taken from numerical recipies

// computes forward-difference approximation to Jacobian
//template <class T>
//struct NRfdjac {
//  const double eps_ = 1.0e-8; // std::sqrt(std::numeric_limits<double>::epsilon());
//  T& func_;
//
//  NRfdjac(T& func) : func_(func) {};
//
//  vec2<double> operator()(vec1<double>& x, vec1<double>& fvec){;
//    int n = x.size();
//    vec2<double> df(n, vec1<double>(n, 0.0));
//    vec1<double> xh = x;
//    for (int j = 0; j < n; ++j){
//      double tmp = xh[j];
//      double h = eps_ * std::abs(tmp);
//      if (h == 0.0){
//        h = eps_;
//      }
//      xh[j] = tmp + h;
//      h = xh[j] - tmp;
//      vec1<double> f = func_(xh);
//      xh[j] = tmp;
//      for (int i = 0; i < n; ++i){
//        df[i][j] = (f[i] - fvec[i])/h;
//      }
//    }
//    return df;
//  }
//};
//
//// returns 1/2 * F * F and stores F in fvec_
//template <class T>
//struct NRfmin {
//  vec1<double> fvec_;
//  T& func_;
//
//  NRfmin(T& func) : func_(func) {}
//
//  double operator() (vec1<double>& x){
//    int n = x.size();
//    double sum = 0;
//    fvec_ = func_(x);
//    for (int i = 0; i < n; ++i){
//      sum += std::sqrt(fvec_[i]);
//    }
//    return 0.5*sum;
//  }
//};

// vec1 = std::vector<double>
// vec2 = std::vector<std::vector<double>>

template <class T>
void broydn2(vec1<double>& x, T& func, const double eta,
             const double TOLF = 1e-8, int clear = 0) {
  typedef vec1<double> vec;
  const int n = x.size();

  vec F = func(x);
  if ((std::abs(*std::max_element(F.begin(), F.end(), [](double a, double b) {
        return std::abs(a) < std::abs(b);
      }))) <= 0.01 * TOLF) {
    std::cout << "Given initial guess is already a minimum\n";
  return;
  }

  vec xold = x;
  vec Fold = F;
  for (int i = 0; i < n; ++i){
    x[i] = xold[i] + eta * F[i];
  }

  //std::cout << eta << "\n";
  //std::cout << F << "\n";
  //std::cout << xold << x << "\n";


  int i = 0;
  bool check = false;
  vec2<double> u;
  vec2<double> v;

  while (not check) {
    F = func(x);
    vec Fdiff(F.size());
    for (int j = 0; j < n; ++j){
      Fdiff[j] = F[j] - Fold[j];
    }

    vec ui = x;
    for (int j = 0; j < n; ++j) {
      ui[j] -= xold[j] + eta * Fdiff[j];
    }

    for (int j = 0; j < i; ++j) {
      double scalar = 0;
      for (int l = 0; l < n; ++l){
        scalar += v[j][l] * Fdiff[l];
      }
      for (int l = 0; l < n; ++l) {
        ui[l] -=  scalar * u[j][l];
      }
    }
    u.push_back(ui);
    vec vi = Fdiff;
    double scalar = 0;
    for (int j = 0; j < n; ++j){
      scalar += Fdiff[j] * Fdiff[j];
    }
    if (scalar == 0){
      std::cerr << "norm in broydn2 is 0.0\n";
      std::exit(1);
    }
    for (int j = 0; j < n; ++j){
      vi[j] = Fdiff[j] / scalar;
    }
    v.push_back(vi);

    xold = x;
    Fold = F;
    for (int j = 0; j < n; ++j){
      x[j] = xold[j] + eta * F[j];
    }

    for (int j = 0; j < i+1; ++j){
      double scalar = 0;
      for (int l = 0; l < n; ++l){
        scalar += v[j][l] * F[l];
      }
      for (int l = 0; l < n; ++l){
        x[l] -=  scalar * u[j][l];
      }
    }
    //std::cout << F << "\n";
    double test =
        std::abs(*std::max_element(F.begin(), F.end(), [](double a, double b) {
          return std::abs(a) < std::abs(b);
        }));
    std::cout << i << ": " << test << "\n";
    if (test < TOLF){
      check = true;
    }

    //std::cout << x << "\n";
    i++;
    if (i == clear){
      u.clear();
      v.clear();
      i = 0;
      clear = 0;
    }
  }

  return;
}

}
