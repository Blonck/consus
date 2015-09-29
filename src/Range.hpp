#pragma once

#include <iostream>

#include "definitions.hpp"

namespace consus
{

template <class T = double>
struct Range{
  T start;
  T step;
  T end;
};

template <class T>
std::ostream& operator<<(std::ostream& out, const Range<T>& range){
  out << range.start << ":" << range.step << ":" << range.end;
  return out;
}

template <class T>
inline void serialize(const Range<T>& t, std::ostream& out){
    serialize(t.start, out);
    serialize(t.step, out);
    serialize(t.end, out);
}

template <class T>
inline void deserialize(Range<T>& t, std::istream& in){
    deserialize(t.start, in);
    deserialize(t.step, in);
    deserialize(t.end, in);
}

} /* end of namespace betamc */ 
