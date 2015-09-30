#pragma once

namespace consus
{

namespace WHAM2D
{

/// H = \beta ( H_0 + \lambda H_1 )
struct NVT{
  static double log_weight(const double Beta, const double E0,
                           const double Lambda, const double E1) {
    return -1.0 * Beta * (E0 + Lambda * E1);
  }
};
  
} /* end of namespace WHAM2D */ 
  
} /* end of namespace consus */ 
