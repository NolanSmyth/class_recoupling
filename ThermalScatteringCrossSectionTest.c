// This file contains code for computing the thermal scattering cross section
// for a model with a Majorana DM fermion, a Majorana DR fermion and a scalar
// mediator. The interaction term looks like: L > g_{ij} * \xi_i \xi_j S + h.c.

#include "ThermalScatteringCrossSections.h"
#include <float.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_dilog.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

int main(int argc, char *argv[])
{
  double g = 1e-2;
  double r = 25.0;
  double dm = 1.01;
  double lam = 1e4;

  if (argc > 1)
  {
    g = atof(argv[1]);
  }
  {
    if (argc > 2)
      lam = atof(argv[2]);
  }
  if (argc > 3)
  {
    r = atof(argv[3]);
  }
  if (argc > 4)
  {
    dm = atof(argv[4]);
  }

  printf("running with: %.16e, %.16e, %.16e, %.16e\n", g, lam, r, dm);

  double w = g * g / (2 * M_PI * r);
  struct Model model = {.dm = dm, .r = r, .g = g, .lam = lam};

  const int ws_size = 1000;

  int size = 500;
  double log_tmin = -4;
  double log_tmax = 8.0;
  double step = (log_tmax - log_tmin) / ((double)size);

  double T;
  double error = 0.0;
  double val1 = 0.0;
  double val2 = 0.0;

  FILE *file = fopen("data2.csv", "w+");

  struct IntegrationResults result = {};

  for (int i = 0; i < size; i++)
  {
    T = pow(10.0, log_tmin + i * step);
    thermal_scattering_rate_components(T, &model, &result);
    fprintf(file, "%f, %.16e \n", T, result.low + result.mid + result.high + result.resonant);
  }
  fclose(file);

  return 0;
}
