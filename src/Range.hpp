#pragma once

#include <iostream>

#include "definitions.hpp"

namespace consus
{

struct Range{
  double start;
  double step;
  double end;
};

std::ostream& operator<<(std::ostream& out, const Range& range){
  out << range.start << ":" << range.step << ":" << range.end;
  return out;
}

inline void serialize(const Range& t, std::ostream& out){
    serialize(t.start, out);
    serialize(t.step, out);
    serialize(t.end, out);
}

inline void deserialize(Range& t, std::istream& in){
    deserialize(t.start, in);
    deserialize(t.step, in);
    deserialize(t.end, in);
}

} /* end of namespace betamc */ 
