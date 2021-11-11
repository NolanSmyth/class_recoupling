// This file contains code for computing the thermal scattering cross section
// for a model with a Majorana DM fermion, a Majorana DR fermion and a scalar
// mediator. The interaction term looks like: L > g_{ij} * \xi_i \xi_j S + h.c.

#include <float.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_sf_dilog.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

//#define DEBUG
static const size_t WS_CQUAD_SIZE = 100;
// value of q at which we stop integrating and use approximate result valid for
// large q
static const double LARGE_Q_SWITCH = 10.0;

gsl_integration_cquad_workspace *get_cquad_integration_workspace() {
  return gsl_integration_cquad_workspace_alloc(WS_CQUAD_SIZE);
}

/**
 * Compute the derivative of the Fermi-Dirac distribution w.r.t. momentum.
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
 * @param p The momentum of the dark-radiation
 * @param T The temperature of the equilibrium bath
 */
double bose_einstein_deriv(double p, double T) {
  const double e = exp(-p / T);
  return (-1.0 / T) * e / pow(1.0 - e, 2);
}

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
 * Structure to hold parameters needed to perform integration of thermal
 * scattering rate.
 */
struct IntegrationParams {
  // Parameters of the model
  struct Model *model;
  // Dark-radiation temperature
  double T;
};

/**
 * Compute the real part of a scalar propagator.
 *
 * @param s invariant momentum: s = p^2.
 * @param m mass of the scalar.
 * @param w width of the scalar.
 */
double scalar_propagator_re(double s, double m, double w) {
  const double m2 = m * m;
  const double w2 = w * w;
  const double smm2 = (s - m2);
  return (s - m2) / (pow(s - m2, 2) + m2 * w2);
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

/**
 * Compute the squared matrix element integrated againts the first two Legendre
 * polynomials: A_0(p) - A_0(q) for a model with Majorana DM, Majorana DR, and
 * scalar mediator.
 * @param q The dark-radiation momentum scaled by mass-splitting: q = p / dm.
 * @param model Parameters of the model.
 */
double tavg_msqrd_ss(double q, struct Model *model) {
  const double r = model->r;
  const double w = model->w;
  const double d = model->dm;

  const double num =
      4 * pow(q, 6) * pow(r, 4) * pow(1 - 2 * (-1 + q) * r, 2) * pow(d, 4);
  const double den =
      pow(2 * q + r, 2) *
      pow(1 + pow(w, 2) + pow(r, 2) * (4 - 8 * q + 4 * pow(q, 2) + pow(w, 2)) +
              r * (4 - 4 * q + 2 * pow(w, 2)),
          2);

  return num / den;
}

double tavg_msqrd_uu(double q, struct Model *model) {
  const double r = model->r;
  const double w = model->w;
  const double d = model->dm;

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

double tavg_msqrd_su(double q, struct Model *model) {
  const double r = model->r;
  const double w = model->w;
  const double d = model->dm;

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

double tavg_msqrd(double q, struct Model *model) {
  return tavg_msqrd_ss(q, model) + tavg_msqrd_uu(q, model) +
         tavg_msqrd_su(q, model);
}

double rate_prefactor(struct Model *model) {
  const double r = model->r;
  const double d = model->dm;
  return 1.0 / (96 * pow(M_PI * d * r, 3));
}

/**
 * Intergrand of the non-resonant thermal-scattering rate.
 *
 * @param q scaled dark-radiation momentum: q = p / dm.
 * @param params structure holding the integration and model parameters.
 */
double thermal_scattering_rate_integrand(double q, void *params) {
  struct IntegrationParams *pars = (struct IntegrationParams *)params;
  const double jac = pars->model->dm;
  const double pre = pars->model->lam * pow(pars->model->g, 2);
  const double msqrd = tavg_msqrd(q, pars->model);
  const double ker = fermi_dirac_der(q * pars->model->dm, pars->T);
  return -jac * msqrd * ker;
}

/**
 * Intergrand of the non-resonant thermal-scattering rate.
 *
 * @param q scaled dark-radiation momentum: q = p / dm.
 * @param params structure holding the integration and model parameters.
 */
double thermal_scattering_rate_large_q(double T, double pre,
                                       struct Model *model) {
  const double r = model->r;
  const double g = model->g;
  const double dm = model->dm;
  const double lam = model->lam;
  const double d = dm / T;
  const double x = r * d;
  const double xp = LARGE_Q_SWITCH * x;
  const double ex = exp(-xp);
  return 3 * pre * lam * pow(g * x * dm * dm, 2) / (8 * T * pow(d, 5)) *
         (-4 * gsl_sf_dilog(-ex) +
          xp * (xp + 4 * log(1 + ex) - xp * tanh(xp / 2)));
}

/**
 * The resonant piece of the thermal scattering rate.
 *
 * @param T temperature.
 * @param params structure holding the model parameters.
 */
double resonant_thermal_scattering_rate(double T, struct Model *params) {
  const double d = params->dm;
  const double r2 = params->r * params->r;
  const double lam = params->lam;
  const double dt = d / T;
  return 1 / (48.0 * M_PI) * (lam * d / r2) * dt / (cosh(dt) + 1.0);
}

/**
 * Compute Gamma_heat / a
 * @param T Temperature of the dark-radiation
 * @param model Structure holding the model parameters
 * @param ws_size Size of the workspace
 * @param ws The workspace for the GSL integrator
 * @param error Pointer to where the error should be stored. If NULL, error is
 * not returned
 */
double thermal_scattering_rate(double T, struct Model *model, double *error) {
  static int got_ws = 0;
  static gsl_integration_cquad_workspace *ws = NULL;

  if (got_ws == 0) {
    ws = get_cquad_integration_workspace();
    got_ws = 1;
  }

  // GSL will abort on errors like roundoff which it shouldn't. Just turn the
  // errors off. We will handle them ourselves.
  gsl_set_error_handler_off();

  // Extract model parameters
  const double r = model->r;
  const double g = model->g;
  const double dm = model->dm;
  const double lam = model->lam;
  const double d = dm / T;
  const double x = r * d;

  // pre is the prefactor in font of Gamma_heat / a. This factor is
  // 1 /(48 * pi^3 * gx * mx^3). Since we are integrating over q = p/dm, we get
  // 4 extra factors of dm. Using mx = dm * r, we get
  // dm^2 / (48 * gx * pi^3 * r^3 ). See ArXiv 1706.07433 Eqn.(6).
  const double pre = rate_prefactor(model);

  // Integrator parameters
  // TODO This relative accuracy might be overkill. Could probly be more like
  // 1e-3.
  const double epsrel = 1e-7;
  const double epsabs = 0.0;
  // breakpt is the rough location of the peak (basically dm for small w)
  const double breakpt_peak = model->dm * (1.0 + pow(model->w / 2.0, 2));
  const double breakpt_large_r = LARGE_Q_SWITCH * r;
  struct IntegrationParams params = {.model = model, .T = T};

  gsl_function F;
  F.function = &thermal_scattering_rate_integrand;
  F.params = &params;

  // Split integral into two piece: first from 0 to peak of the Briet-Wigner and
  // the second from the peak to large r

  // Perform integration from 0 to peak
  double res1 = 0.0;
  double err1 = 0.0;
  size_t nevals1 = 0;
  gsl_integration_cquad(&F, 0.0, breakpt_peak, epsabs, epsrel, ws, &res1, &err1,
                        &nevals1);
  // Perform integral from peak to a cut-off for large r
  double res2 = 0.0;
  double err2 = 0.0;
  size_t nevals2 = 0;
  gsl_integration_cquad(&F, breakpt_peak, breakpt_large_r, epsabs, epsrel, ws,
                        &res2, &err2, &nevals2);

  // Scale the results by the prefactor
  res1 *= pre;
  err1 *= pre;
  res2 *= pre;
  err2 *= pre;

  // Make sure we aren't setting a null-pointer.
  if (error != NULL) {
    *error = sqrt(err1 * err1 + err2 * err2);
  }

  const double resonant = resonant_thermal_scattering_rate(T, model);

  // Appoximate result for integral from q >~r to infinity
  const double large = thermal_scattering_rate_large_q(T, pre, model);

  return resonant + res1 + res2 + large;
}

int main() {
  double g = 1e-2;
  double r = 25.0;
  double dm = 1.0;
  double lam = 1e4;
  double w = g * g / (2 * M_PI * r);
  struct Model model = {.dm = dm, .w = w, .r = r, .g = g, .lam = lam};

  const int ws_size = 1000;

  int size = 500;
  double log_tmin = -4.0;
  double log_tmax = 8.0;
  double step = (log_tmax - log_tmin) / ((double)size);

  double T;
  double error = 0.0;
  double val1 = 0.0;
  double val2 = 0.0;

  FILE *file = fopen("data3.csv", "w+");

  for (int i = 0; i < size; i++) {
    T = pow(10.0, log_tmin + i * step);
    val1 = thermal_scattering_rate(T, &model, &error);
    fprintf(file, "%f, %.16e\n", T, val1);
  }
  fclose(file);

  return 0;
}
