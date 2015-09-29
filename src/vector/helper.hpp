#pragma once
#include <vector>
#include <iostream>
#include <sstream>

#include "../definitions.hpp"

typedef std::vector<prec_t> vec1d;
typedef std::vector<std::vector<prec_t>> vec2d;
typedef std::vector<std::vector<std::vector<prec_t>>> vec3d;

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
std::ostream& operator<<(std::ostream& out, const std::vector<T>& data){
  std::cout << "[ ";
  for (auto &elem: data){
    std::cout << elem << " ";
  }
  std::cout << "]";
  return out;
}

template <class T>
std::ostream& operator<<(std::ostream& out, const std::vector<std::vector<T>>& data){
  std::cout << "[";
  for (size_t i = 0; i < data.size(); ++i){
    if ( not (i == 0)){
      std::cout << " ";
    }
    std::cout << "[ ";
    for (size_t j = 0; j < data[i].size(); ++j){
      std::cout << data[i][j] << " ";
    }
    std::cout << "]";
    if ( i != (data.size()-1)){
      std::cout << "\n";
    }
  }
  std::cout << "]";
  return out;
}
