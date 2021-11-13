#include "Tools.h"

namespace thermal_scattering_rates {

double fermi_dirac_der(double p, double T) {
  const double e = exp(-p / T);
  return (-1.0 / T) * e / pow(1.0 + e, 2);
}

double bose_einstein_deriv(double p, double T) {
  const double e = exp(-p / T);
  return (-1.0 / T) * e / pow(1.0 - e, 2);
}

} // namespace thermal_scattering_rates
