#pragma once

#include "../definitions.hpp"

#include <limits>
#include <iostream>

namespace consus
{

namespace WHAM2D
{

/// simple class containing some information about DiscreteAxis2D objects
struct HistInfo2D{
  /// logarithm of sum of all entries
  double log_length;
  /// first filled index for first dimension
  int ffi_first = std::numeric_limits<int>::max();
  /// last filled index for first dimension
  int lfi_first = std::numeric_limits<int>::lowest();
  /// first filled index for second dimension
  int ffi_second = std::numeric_limits<int>::max();
  /// last filled index for second dimension
  int lfi_second = std::numeric_limits<int>::lowest();

  /// empty ctor
  HistInfo2D(){}
  HistInfo2D(double log_length_, int ffi_first_, int lfi_first_,
             int ffi_second_, int lfi_second_)
      : log_length(log_length_),
        ffi_first(ffi_first_),
        lfi_first(lfi_first_),
        ffi_second(ffi_second_),
        lfi_second(lfi_second_) {}
};


inline void serialize(const WHAM2D::HistInfo2D& t, std::ostream& out){
  dlib::serialize(t.log_length, out);
  dlib::serialize(t.ffi_first, out);
  dlib::serialize(t.lfi_first, out);
  dlib::serialize(t.ffi_second, out);
  dlib::serialize(t.lfi_second, out);
}


inline void deserialize(WHAM2D::HistInfo2D& t, std::istream& in){
  dlib::deserialize(t.log_length, in);
  dlib::deserialize(t.ffi_first, in);
  dlib::deserialize(t.lfi_first, in);
  dlib::deserialize(t.ffi_second, in);
  dlib::deserialize(t.lfi_second, in);
}

} /* end of namespace WHAM2D */ 
  
} /* end of namespace consus */ 
