#pragma once

#include <vector>
#include <iostream>
#include <sstream>

#include "../definitions.hpp"

namespace consus
{
  

typedef std::vector<double> vec1d;
typedef std::vector<std::vector<double>> vec2d;
typedef std::vector<std::vector<std::vector<double>>> vec3d;

template <class T>
using vec1 = std::vector<T>;

template <class T>
using vec2 = std::vector<std::vector<T>>;

template <class T>
using vec3 = std::vector<std::vector<std::vector<T>>>;

constexpr size_t inner_size(double f){
  return 0;
}

inline size_t inner_size(const vec1d& data){
  return data.size();
}

inline size_t inner_size(const vec2d& data){
  return data[0].size();
}

inline size_t inner_size(const vec3d& data){
  return data[0][0].size();
}

inline size_t inner_p1_size(const vec1d& data){
  return 0;
}

inline size_t inner_p1_size(const vec2d& data){
  return data.size();
}

inline size_t inner_p1_size(const vec3d& data){
  return data[0].size();
}

template <class T>
struct dim{
  static constexpr size_t value = 0;
};

template <>
struct dim<vec1d>{
  static constexpr size_t value = 1;
};

template <>
struct dim<vec2d>{
  static constexpr size_t value = 2;
};

template <>
struct dim<vec3d>{
  static constexpr size_t value = 3;
};

template <class T>
struct inner_type{
  typedef typename T::value_type value_type;
};

template <>
struct inner_type<vec1d>{
  typedef typename vec1d::value_type value_type;
};

template <>
struct inner_type<vec2d>{
  typedef typename vec2d::value_type value_type;
};

template <class T>
std::ostream& operator<<(std::ostream& out, const vec1<T>& data){
  out << "[ ";
  for (auto &elem: data){
    out << elem << " ";
  }
  out << "]";
  return out;
}

template <class T>
std::string print_1d(const vec1<T>& vec){
  std::stringstream ss;
  ss << "[ ";
  for (size_t i = 0; i < vec.size()-1; ++i){
    ss << vec[i] << ", ";
  }
  ss << vec.back() << "]";
  return ss.str();
}

template <class T>
std::ostream& operator<<(std::ostream& out, const vec2<T>& data){
  out << "[";
  for (size_t i = 0; i < data.size(); ++i){
    if ( not (i == 0)){
      out << " ";
    }
    out << "[ ";
    for (size_t j = 0; j < data[i].size(); ++j){
      out << data[i][j] << " ";
    }
    out << "]";
    if ( i != (data.size()-1)){
      out << "\n";
    }
  }
  out << "]";
  return out;
}

} /* end of namespace consus */ 
