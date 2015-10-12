#pragma once

#include <cmath>

namespace consus
{

namespace WHAM2D
{


/// Ensemble class for NVT (canonical) ensemble
///
/// functions calculating the DOS and reweighting the data use
/// this class to adept WHAM to the canonical ensemble
/// Hamiltonian looks like H = H_0 + \lamda H_1
struct NVT{

  /// logarithm of the weight of the canonical ensemble
  ///
  /// \param Beta inverse temperature
  /// \param E0 first energy term of the Hamiltonian
  /// \param Lambda prefactor of the second energy term of the Hamiltonian
  /// \param E1 second energy term of the Hamiltonian
  static double log_weight(const double Beta, const double E0,
                           const double Lambda, const double E1) {
    return -1.0 * Beta * (E0 + Lambda * E1);
  }

  /// the Hamiltonian of the system
  ///
  /// \param E0 first energy term of the Hamiltonian
  /// \param Lambda prefactor of the second energy term of the Hamiltonian
  /// \param E1 second energy term of the Hamiltonian
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

// logW = (3*N - 2)/2 * log(E_t - (E_p + p V))
struct NPH{
  static double N;
  static double log_weight(const double E_t, const double E_p,
                           const double p, const double logV){
    return (1.5*N - 1.0) * std::log(E_t - (E_p + p * std::exp(logV)));
  }
};

double NPH::N = 1.0;

} /* end of namespace WHAM2D */ 
  
} /* end of namespace consus */ 
