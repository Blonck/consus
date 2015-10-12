#pragma once

#include "../vector/helper.hpp"

namespace consus
{

/// returns all indicies where data[i] and data[i+1] have a difference
/// according to a predicate
///
/// returns a std::vector<size_t>
/// \param data - 1 dim vector
/// \param func - predicate which returns true if data[i] and data[i+1] differs
template <class T>
std::vector<size_t> where_diff(const vec1d& data, const T& func) {
  std::vector<size_t> res;
  for (size_t i = 0; i < data.size() - 1; ++i) {
    if (func(data[i], data[i + 1])) {
      res.push_back(i);
    }
  }
  return res;
}

/// returns all indicies where data[j][i] and data[j][i+1] have a difference
/// according to a predicate
///
/// returns a std::vector<std::vector<size_t>>
/// \param data - 2 dim vector
/// \param func - predicate which returns true if data[i] and data[i+1] differs
template <class T>
std::vector<std::vector<size_t>> where_diff(const vec2d& data, const T& func) {
  std::vector<std::vector<size_t>> res(data.size());
  for (size_t j = 0; j < data.size(); ++j) {
    for (size_t i = 0; i < data[j].size() - 1; ++i) {
      if (func(data[j][i], data[j][i + 1])) {
        res[j].push_back(i);
      }
    }
  }
  return res;
}

//template <class T>
//class FWhere_diff : public FBase {
// public:
//  typedef std::vector<int> result_type;
//  const std::vector<int> init_element;
//
//  T func_;
//  std::vector<size_t> res;
//
// public:
//  FWhere_diff(const T& func): func_(func) {};
//
//  inline void operator()(const vec1d& data, const size_t i) {
//    if (i != 0) {
//      if (func(data[i - 1], data[i])) {
//        res.push_back(i - 1);
//      }
//    }
//  }
//
//  inline std::vector<size_t> result() {
//    return res;
//  }
//
//  inline void reset() {
//    res.resize(0);
//  }
//};
//
//template <class T>
//class SWhere_diff : public SBase {
// public:
//  typedef std::vector<std::vector<int>> result_type;
//  const std::vector<std::vector<int>> init_element;
//
//  T func_;
//  std::vector<std::vector<size_t>> res;
//
// public:
//  SWhere_diff(const vec2d& data, const T& func): init_element(data.size), (func_(func) {};
//
//  inline void operator()(const vec1d& data, const size_t i) {
//    for (size_t j = 0; j < data.size(); ++j) {
//      if (i != 0) {
//        if (func(data[j][i - 1], data[j][i])) {
//          res[j].push_back(i - 1);
//        }
//      }
//    }
//  }
//
//  inline std::vector<std::vector<size_t>> result() {
//    return res;
//  }
//
//  inline void reset() {
//    for (auto& elem : res) {
//      elem.resize(0);
//    }
//  }
//};

} /* end of namespace consus */ 
