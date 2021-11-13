// This file contains code for computing the thermal scattering cross section
// for a model with a Majorana DM fermion, a Majorana DR fermion and a scalar
// mediator. The interaction term looks like: L > g_{ij} * \xi_i \xi_j S + h.c.

#include "thermal_scattering_rates/ThermalScatteringRate.h"
#include "thermal_scattering_rates/Interface.h"
#include <boost/math/quadrature/gauss_kronrod.hpp>
#include <gsl/gsl_sf_dilog.h>

// ===========================================================================
// ---- Tools ----------------------------------------------------------------
// ===========================================================================

/**
 * Compute the derivative of the Fermi-Dirac distribution w.r.t. momentum.
 *
 * @param p The momentum of the dark-radiation
 * @param T The temperature of the equilibrium bath
 */
double fermi_dirac_der(double p, double T) {
  const double e = exp(-p / T);
  return (-1.0 / T) * e / pow(1.0 + e, 2);
}

/**
 * Compute the derivative of the Bose-Einstein distribution w.r.t. momentum with
 * factor of 1/2T removed.
 *
 * @param p The momentum of the dark-radiation
 * @param T The temperature of the equilibrium bath
 */
double bose_einstein_deriv(double p, double T) {
  const double e = exp(-p / T);
  return (-1.0 / T) * e / pow(1.0 - e, 2);
}

// ===========================================================================
// ---- t-Averaged Squared Matrix Elements -----------------------------------
// ===========================================================================

// The t-Average of a squared matrix element for the 2->2 process
// A(p') + B(k') -> A(p) + B(k) with A non-relativistic is given by:
//  <M2>_t = 1 / (8 * |k|^4) * Integrate[-t * M2(t), {t, -4 kcm2, 0}]
// where:
//  s = (p1 + p2)^2 = mA^2 + mB^2 + 2 * mA * EB
//  t = (p1 - p3)^2
//  kcm2 = (s-(mA+mB)^2) * (s-(mA-mB)^2) / (2s)

// S-channel only contribution to <|M|^2>
auto Model::msqrd_ss_tavg(double q) const -> double {
  const double r = p_r;
  const double w = p_w;
  const double d = p_delta;

  const double num =
      4 * pow(q, 6) * pow(r, 4) * pow(1 - 2 * (-1 + q) * r, 2) * pow(d, 4);
  const double den =
      pow(2 * q + r, 2) *
      pow(1 + pow(w, 2) + pow(r, 2) * (4 - 8 * q + 4 * pow(q, 2) + pow(w, 2)) +
              r * (4 - 4 * q + 2 * pow(w, 2)),
          2);

  return num / den;
}

// U-channel only contribution to <|M|^2>
auto Model::msqrd_uu_tavg(double q) const -> double {
  const double r = p_r;
  const double d = p_delta;

  const double num = -(
      pow(d, 3) * ((-8 * pow(q, 2) * r *
                    (2 * q + r + q * (4 + q) * r + 2 * pow(r, 2)) * d) /
                       pow(2 * q + r, 2) -
                   (4 * r * pow(q + 2 * q * r, 2) * d) /
                       (r + 2 * (pow(r, 2) + q * pow(1 + r, 2))) +
                   (3 + (6 + 4 * q) * r) * (d + 2 * r * d) *
                       log(1 + (4 * pow(q, 2) * r) /
                                   (r + 2 * (pow(r, 2) + q * pow(1 + r, 2))))));
  const double den = 8;

  return num / den;
}

// s-u interference contribution to <|M|^2>
auto Model::msqrd_su_tavg(double q) const -> double {
  const double r = p_r;
  const double w = p_w;
  const double d = p_delta;

  const double num =
      r * (1 - 2 * (-1 + q) * r) * pow(d, 4) *
      ((-4 * pow(q, 2) * r *
        (2 * pow(q, 2) + q * (3 + 2 * q * (2 + q)) * r +
         (1 + 6 * q * (1 + q)) * pow(r, 2) + 2 * (1 + q) * pow(r, 3))) /
           pow(2 * q + r, 2) +
       (1 + 2 * (1 + q) * r) * (q + r + 2 * q * r + 2 * (1 + q) * pow(r, 2)) *
           log(1 + (4 * pow(q, 2) * r) /
                       (r + 2 * (pow(r, 2) + q * pow(1 + r, 2)))));
  const double den =
      4 * (pow(1 - 2 * (-1 + q) * r, 2) + pow(1 + r, 2) * pow(w, 2));

  return num / den;
}

auto Model::msqrd_tavg(double q) const -> double {
  return msqrd_ss_tavg(q) + msqrd_uu_tavg(q) + msqrd_su_tavg(q);
}

// ===========================================================================
// ---- Thermal Scattering Rate ----------------------------------------------
// ===========================================================================

auto Model::thermal_scattering_rate_integrand(double q, double T) const
    -> double {
  const double msqrd = msqrd_tavg(q);
  const double ker = fermi_dirac_der(q * p_delta, T);
  return -p_delta * msqrd * ker;
}

// Perform integration from 0 to peak
auto Model::thermal_scattering_rate_small_q(double T) const -> double {
  using boost::math::quadrature::gauss_kronrod;

  const double breakpt_peak = p_delta * (1.0 + pow(p_w / 2.0, 2));

  auto integrand = [&T, this](double q) {
    return thermal_scattering_rate_integrand(q, T);
  };

  return p_pre * gauss_kronrod<double, 15>::integrate(integrand, 0.0,
                                                      breakpt_peak, 5, 1e-9);
}

// Perform integral from peak to a cut-off for large r
auto Model::thermal_scattering_rate_mid_q(double T) const -> double {
  using boost::math::quadrature::gauss_kronrod;

  const double breakpt_peak = p_delta * (1.0 + pow(p_w / 2.0, 2));
  const double breakpt_large_r = p_switch_factor * p_r;

  auto integrand = [&T, this](double q) {
    return thermal_scattering_rate_integrand(q, T);
  };
  return p_pre * gauss_kronrod<double, 15>::integrate(integrand, breakpt_peak,
                                                      breakpt_large_r, 5, 1e-9);
}

// Contribution to thermal-scattering-rate valid for q >> r
auto Model::thermal_scattering_rate_large_q(double T) const -> double {
  const double r = p_r;
  const double g = p_g;
  const double delta = p_delta;
  const double lam = p_lam;
  const double d = delta / T;
  const double x = r * d;
  const double xp = p_switch_factor * x;
  const double ex = exp(-xp);
  return 3 * p_pre * lam * sqr(g * x * sqr(delta)) / (8 * T * pow(d, 5)) *
         (-4 * gsl_sf_dilog(-ex) +
          xp * (xp + 4 * log(1 + ex) - xp * tanh(xp / 2)));
}

// Contribution to thermal-scattering-rate from the resonance.
auto Model::thermal_scattering_rate_resonant(double T) const -> double {
  const double dt = p_delta / T;
  return 1 / (48.0 * M_PI) * (p_lam * p_delta / sqr(p_r)) * dt /
         (cosh(dt) + 1.0);
}

auto Model::thermal_scattering_rate(double T) const -> double {
  return thermal_scattering_rate_small_q(T) + thermal_scattering_rate_mid_q(T) +
         thermal_scattering_rate_large_q(T) +
         thermal_scattering_rate_resonant(T);
}

auto thermal_scattering_rate(double T, double delta, double r, double g,
                             double lam) -> double {
  const Model model{delta, r, g, lam};
  return model.thermal_scattering_rate(T);
}
