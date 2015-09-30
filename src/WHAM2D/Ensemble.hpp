#pragma once

#include <cmath>

namespace consus
{

namespace WHAM2D
{

/// logW = - \beta ( H_0 + \lambda H_1 )
struct NVT{
  static double log_weight(const double Beta, const double E0,
                           const double Lambda, const double E1) {
    return -1.0 * Beta * (E0 + Lambda * E1);
  }

  static double ham(const double E0, const double Lambda, const double E1) {
    return E0 + Lambda * E1;
  }
};

// logW = - \beta ( E + p logV)
struct NPT{
  static double log_weight(const double Beta, const double E,
                           const double p, const double logV) {
    return -1.0 * Beta * (E + p * std::exp(logV));
  }

  static double ham(const double E, const double p, const double logV){
    return E + p * std::exp(logV);
  }
};
  
} /* end of namespace WHAM2D */ 
  
} /* end of namespace consus */ 
