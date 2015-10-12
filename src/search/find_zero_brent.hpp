#pragma once

#include <cassert>
#include <limits>
#include <cmath>
#include <algorithm>
#include <tuple>

namespace consus
{

struct search_result {
  double x;
  double fx;
  bool valid;

  search_result(const double x_value_, const double fx_value_,
                const bool valid_)
      : x(x_value_), fx(fx_value_), valid(valid_) {}
};

template <class T>
search_result find_zero_brent_impl(T& func, const double lower_limit,
                                   const double upper_limit,
                                   const double error_tol) {
  double a = lower_limit;
  double b = upper_limit;
  double d = std::numeric_limits<double>::max();
  double fa = func(a);
  double fb = func(b);
  double s = 0;
  double fs = 0;
  if (std::abs(fa) < std::abs(fb)){
    std::swap(a, b);
    std::swap(fa, fb);
  }

  double c = a;
  double fc = fa;
  double mflag = true;
  
  size_t i =0;
  while( (not(fb == 0.0)) and std::abs(a - b) > error_tol) {
    if ((fa != fc) and (fb != fc)){
      // Inverse quadratic interpolation
      s = a * fb * fc / (fa - fb) / (fa - fc) +
          b * fa * fc / (fb - fa) / (fb - fc) +
          c * fa * fb / (fc - fa) / (fc - fb);
    }else{
      // Secant Rule
      s = b - fb * (b - a) / (fb - fa);
    }
    double tmp2 = (3 * a + b) / 4.0;
    if (((not(((s > tmp2)and(s < b))or((s < tmp2)and(s > b)))))or(
            mflag and(std::abs(s - b) >= (std::abs(b - c) / 2.0)))
        or(not(mflag) and(std::abs(s - b) >= (std::abs(c - d) / 2.0)))) {
      s = (a + b) / 2.0;
      mflag = true;
    } else {
      if (((mflag and(std::abs(b - c) < error_tol)))or(
              not(mflag) and(std::abs(c - d) < error_tol))) {
        s = (a + b) / 2.0;
        mflag = true;
      } else {
        mflag = true;
      }
    }
    fs = func(s);
    d = c;
    c = b;
    fc = fb;
    if (fa * fs < 0) {
      b = s;
      fb = fs;
    } else {
      a = s;
      fa = fs;
    }

    if (std::abs(fa) < std::abs(fb)) {
      std::swap(a, b);
      std::swap(fa, fb);
      i = i + 1;
    }
    if (i > 1000) {
      //std::cerr << "to much tries in brent algorithm\n";
      return search_result(b, fb, false);
    }
  }
  return search_result(b, fb, true);
}

template <class T>
search_result find_zero_brent(T& func, double lower_bound, double upper_bound,
                              const double error_tol = 1e-9) {
  if (lower_bound > upper_bound) {
    std::swap(lower_bound, upper_bound);
  }
  double old_left = func(lower_bound);
  double old_right = func(upper_bound);
  if (std::signbit(old_left) == std::signbit(old_right)) {
    //std::cerr << "bracket have no zero crossing\n";
    return search_result(lower_bound, old_left, false);
  }
  assert(std::signbit(old_left) != std::signbit(old_right));
  auto u = find_zero_brent_impl(func, lower_bound, upper_bound, error_tol);
  return u;
}

template <class T>
search_result find_zero_brent(T& func, double lower_bound, double upper_bound,
                              double min_lower_bound, double max_upper_bound,
                              const double error_tol = 1e-9,
                              const double attempt_step = 0.005,
                              const size_t attempts = 50) {
    if (lower_bound > upper_bound){
      std::swap(lower_bound, upper_bound);
      std::swap(min_lower_bound, max_upper_bound);
    }
    assert(min_lower_bound <= lower_bound);
    assert(max_upper_bound >= upper_bound);
    double old_left = func(lower_bound);
    double old_right = func(upper_bound);
    if (std::signbit(old_left) == std::signbit(old_right)) {
      //std::cout << "bracket ]" << lower_bound << ", " << upper_bound
      //          << "[ don't have a crossing, try to adjust\n";
      bool min_reached = false;
      bool max_reached = false;
      for (size_t i = 0; i < attempts; ++i) {
        if (min_lower_bound <= lower_bound - attempt_step){
          lower_bound -= attempt_step;
          old_left = func(lower_bound);
        }else{
          min_reached = true;
        }
        if (max_upper_bound >= upper_bound + attempt_step){
          upper_bound += attempt_step;
          old_right = func(upper_bound);
        }else{
          max_reached = true;
        }
        if (std::signbit(old_left) != std::signbit(old_right)
            or(min_reached and max_reached)) {
          break;
        }
      }
    }
    if (std::signbit(old_left) == std::signbit(old_right)) {
      //std::cerr << "bracket have no zero crossing\n";
      return search_result(lower_bound, old_left, false);
    }
    assert(std::signbit(old_left) != std::signbit(old_right));
    auto u = find_zero_brent_impl(func, lower_bound, upper_bound, error_tol);
    return u;
}

} /* end of namespace consus */ 
