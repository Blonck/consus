#pragma once

#include "../extern/dlib/dlib/serialize.h"

namespace consus{

//import the serialize/deserialize function of dlib into betamc namespace
using dlib::serialize;
using dlib::deserialize;

typedef double prec_t;
}
