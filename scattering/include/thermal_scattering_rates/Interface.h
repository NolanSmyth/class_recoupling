#ifndef THERMAL_SCATTERING_RATE_INTERFACE_H
#define THERMAL_SCATTERING_RATE_INTERFACE_H

#ifdef __cplusplus
extern "C" {
#endif

double thermal_scattering_rate_mms(double T, double delta, double r, double g,
                                   double lam);

#ifdef __cplusplus
}
#endif

#endif // THERMAL_SCATTERING_RATE_INTERFACE_H
