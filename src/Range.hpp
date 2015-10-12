#pragma once

#include <iostream>

#include "definitions.hpp"

namespace consus
{

/// simple struct describing a floating point range [start, end[
///
/// used to construct DiscreteAxis objects
struct Range{
  double start;
  double step;
  double end;
};

/// ostream operator for Range
std::ostream& operator<<(std::ostream& out, const Range& range){
  out << range.start << ":" << range.step << ":" << range.end;
  return out;
}

/// serialize Range using dlib
inline void serialize(const Range& t, std::ostream& out){
    serialize(t.start, out);
    serialize(t.step, out);
    serialize(t.end, out);
}

/// deserialize Range using dlib
inline void deserialize(Range& t, std::istream& in){
    deserialize(t.start, in);
    deserialize(t.step, in);
    deserialize(t.end, in);
}

} /* end of namespace betamc */ 
