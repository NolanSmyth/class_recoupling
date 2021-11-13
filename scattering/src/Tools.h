#ifndef THERMAL_SCATTERING_RATES_TOOLS_H
#define THERMAL_SCATTERING_RATES_TOOLS_H

#include <cmath>

namespace thermal_scattering_rates {

template <class T> constexpr auto sqr(T x) -> decltype(x * x) { return x * x; }

/**
 * Compute the derivative of the Fermi-Dirac distribution w.r.t. momentum.
 *
 * @param p The momentum of the dark-radiation
 * @param T The temperature of the equilibrium bath
 */
auto fermi_dirac_der(double p, double T) -> double;

/**
 * Compute the derivative of the Bose-Einstein distribution w.r.t. momentum with
 * factor of 1/2T removed.
 *
 * @param p The momentum of the dark-radiation
 * @param T The temperature of the equilibrium bath
 */
auto bose_einstein_deriv(double p, double T) -> double;

} // namespace thermal_scattering_rates

#endif // THERMAL_SCATTERING_RATES_TOOLS_H
