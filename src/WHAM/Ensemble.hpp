#pragma once

#include <cmath>
#include <utility>

namespace consus
{

namespace WHAM
{

/// Ensemble class for NVT (canonical) ensemble
///
/// functions calculating the DOS and reweighting the data use
/// this class to adept WHAM to the canonical ensemble
struct NVT{

  /// logarithm of the weight of the canonical ensemble
  ///
  /// \param Beta inverse temperature
  /// \param E energy
  static double log_weight(const double Beta, const double E) {
    return -1.0 * Beta * E;
  }
};

} /* end of namespace WHAM */ 
  
}
