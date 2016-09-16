#pragma once

#include "../definitions.hpp"

namespace consus
{

namespace WHAM2D
{

/// simple class containing some information about DiscreteAxis2D objects
struct HistInfo{
  /// logarithm of sum of all entries
  double log_length;
  /// first filled index
  int ffi = std::numeric_limits<int>::max();
  /// last filled index
  int lfi = std::numeric_limits<int>::lowest();

  /// empty ctor
  HistInfo(){}
  HistInfo(double log_length_, int ffi_, int lfi_)
      : log_length(log_length_),
        ffi(ffi_),
        lfi(lfi_)
        {}
};


inline void serialize(const WHAM2D::HistInfo& t, std::ostream& out){
  dlib::serialize(t.log_length, out);
  dlib::serialize(t.ffi, out);
  dlib::serialize(t.lfi, out);
}

inline void deserialize(WHAM2D::HistInfo& t, std::istream& in){
  dlib::deserialize(t.log_length, in);
  dlib::deserialize(t.ffi, in);
  dlib::deserialize(t.lfi, in);
}

} /* end of namespace WHAM2D */ 
  
} /* end of namespace consus */ 
