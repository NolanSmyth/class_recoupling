#ifndef THERMAL_SCATTERING_RATES_H
#define THERMAL_SCATTERING_RATES_H

#include <cmath>

namespace thermal_scattering_rates {

class ModelMMS {

private:
  static double scalar_width(double g, double r) {
    return pow(g, 2) / (2 * M_PI * r);
  }
  static auto thermal_scattering_rate_prefactor(double delta, double r,
                                                double g, double lam)
      -> double {
    return g * g * lam / (96 * pow(M_PI * delta * r, 3));
  }

  [[nodiscard]] auto msqrd_ss_tavg(double) const -> double;
  [[nodiscard]] auto msqrd_uu_tavg(double) const -> double;
  [[nodiscard]] auto msqrd_su_tavg(double) const -> double;

  // Mass-splitting between DM and mediator: delta = mMed - m
  double p_delta;
  // Dark matter mass in units of the mass-splitting: r = m / dm
  double p_r;
  // Coupling constant between DM, DR and mediator
  double p_g;
  // g^2 * N
  double p_lam;
  // Width of the mediator in units of the mass-splitting: w = width / dm
  double p_w;
  // Rate prefactor
  double p_pre;

  // `switch_factor * r` is the value of q at which the integration switches to
  // an approximate result valid for large q
  double p_switch_factor = 5.0;

public:
  ModelMMS(double delta, double r, double g, double lam)
      : p_delta(delta), p_r(r), p_g(g), p_lam(lam), p_w(scalar_width(g, r)),
        p_pre(thermal_scattering_rate_prefactor(delta, r, g, lam)) {}

  [[nodiscard]] auto delta() const -> const double & { return p_delta; };
  [[nodiscard]] auto r() const -> const double & { return p_r; };
  [[nodiscard]] auto w() const -> const double & { return p_w; };
  [[nodiscard]] auto g() const -> const double & { return p_g; };
  [[nodiscard]] auto lam() const -> const double & { return p_lam; };
  [[nodiscard]] auto switch_factor() const -> const double & {
    return p_switch_factor;
  };

  auto delta(double val) -> void {
    p_delta = val;
    p_w = scalar_width(p_g, p_r);
    p_pre = thermal_scattering_rate_prefactor(p_delta, p_r, p_g, p_lam);
  };
  auto r(double val) -> void {
    p_r = val;
    p_w = scalar_width(p_g, p_r);
    p_pre = thermal_scattering_rate_prefactor(p_delta, p_r, p_g, p_lam);
  };
  auto g(double val) -> void {
    p_g = val;
    p_w = scalar_width(p_g, p_r);
    p_pre = thermal_scattering_rate_prefactor(p_delta, p_r, p_g, p_lam);
  };
  auto lam(double val) -> void {
    p_lam = val;
    p_w = scalar_width(p_g, p_r);
    p_pre = thermal_scattering_rate_prefactor(p_delta, p_r, p_g, p_lam);
  };
  auto switch_factor(double val) -> void { p_switch_factor = val; };

  /**
   * Compute the t-averaged squared matrix element for DM + DR -> DM + DR.
   * rate.
   *
   * @param q dark-radiation momentum scaled by dark-matter/mediator
   * mass-splitting: q = p / delta.
   */
  [[nodiscard]] auto msqrd_tavg(double) const -> double;

  /**
   * Compute the integrand of the thermal dark-matter/dark-radiation scattering
   * rate.
   *
   * @param q dark-radiation momentum scaled by dark-matter/mediator
   * mass-splitting: q = p / delta.
   * @param T temperature of the thermal bath (dark-radiation temperature)
   */
  [[nodiscard]] auto thermal_scattering_rate_integrand(double, double) const
      -> double;

  [[nodiscard]] auto thermal_scattering_rate_small_q(double) const -> double;
  [[nodiscard]] auto thermal_scattering_rate_mid_q(double) const -> double;
  [[nodiscard]] auto thermal_scattering_rate_large_q(double) const -> double;
  [[nodiscard]] auto thermal_scattering_rate_resonant(double) const -> double;

  /**
   * Compute the thermal dark-matter/dark-radiation scattering rate.
   *
   * @param T temperature of the thermal bath (dark-radiation temperature)
   */
  [[nodiscard]] auto thermal_scattering_rate(double) const -> double;
};

} // namespace thermal_scattering_rates

#endif // THERMAL_SCATTERING_RATES_H
