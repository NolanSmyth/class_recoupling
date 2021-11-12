#ifndef THERMAL_SCATTERING_CROSS_SECTIONS_H
#define THERMAL_SCATTERING_CROSS_SECTIONS_H

/**
 * Structure to hold the model parameters.
 */
struct Model {
  // Mass-splitting between DM and mediator: dm = mMed - m
  double dm;
  // Dark matter mass in units of the mass-splitting: r = m / dm
  double r;
  // Width of the mediator in units of the mass-splitting: w = width / dm
  double w;
  // Coupling constant between DM, DR and mediator
  double g;
  // g^2 * N
  double lam;
};

/**
 * Structure to hold the model parameters.
 */
struct IntegrationResults {
  double low;
  double mid;
  double high;
  double resonant;
  double error;
};

/**
 * Compute Gamma_heat / a
 *
 * @param[] T temperature of the thermal bath.
 * @param[in] model structure holding the model parameters
 * @param[out] error error estimate from GSL.
 */
void thermal_scattering_rate_components(double T, struct Model *model,
                                        struct IntegrationResults *results);

/**
 * Compute Gamma_heat / a
 *
 * @param[] T temperature of the thermal bath.
 * @param[in] model structure holding the model parameters
 * @param[out] error error estimate from GSL.
 */
double thermal_scattering_rate(double T, struct Model *model, double *error);

#endif // THERMAL_SCATTERING_CROSS_SECTIONS_H
