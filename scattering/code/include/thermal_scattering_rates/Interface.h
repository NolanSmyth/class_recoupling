#ifndef THERMAL_SCATTERING_RATE_INTERFACE_H
#define THERMAL_SCATTERING_RATE_INTERFACE_H

#ifdef __cplusplus
extern "C" {
#endif

auto thermal_scattering_rate(double T, double delta, double r, double g,
                             double lam) -> double;

#ifdef __cplusplus
}
#endif

#endif // THERMAL_SCATTERING_RATE_INTERFACE_H
