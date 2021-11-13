// This file contains code for computing the thermal scattering cross section
// for a model with a Majorana DM fermion, a Majorana DR fermion and a scalar
// mediator. The interaction term looks like: L > g_{ij} * \xi_i \xi_j S + h.c.

#include "thermal_scattering_rates/ThermalScatteringRates.h"
#include <pybind11/numpy.h>
#include <pybind11/pybind11.h>

namespace py = pybind11;
namespace tsr = thermal_scattering_rates;

auto add_ctor(py::class_<tsr::ModelMMS> *model) -> void {
  static std::string doc = R"Doc(
Create a new `tsr::ModelMMS` object representing a model with Majorana DM, Majorana DR
and N scalars mediating their interaction.

Parameters
----------
delta: float
  Mass splitting between DM and scalar mediator.
r: float
  Ratio of DM mass to mass splitting: r = (DM-mass) / (mass-splitting).
g: float
  Scalar-DR-DM coupling.
lam: float
  Number of scalars times the coupling squared: lam = g^2 * N. 
)Doc";

  model->def(py::init<double, double, double, double>(), doc.c_str(),
             py::arg("delta"), py::arg("r"), py::arg("g"), py::arg("lam"));
}

auto add_properties(py::class_<tsr::ModelMMS> *model) -> void {
  model->def_property(
      "delta", [](const tsr::ModelMMS &mod) { return mod.delta(); },
      [](tsr::ModelMMS &mod, double val) { mod.delta(val); },
      "Mass splitting between DM and scalar mediator.");
  model->def_property(
      "r", [](const tsr::ModelMMS &mod) { return mod.r(); },
      [](tsr::ModelMMS &mod, double val) { mod.r(val); },
      "Ratio of DM mass to mass splitting: r = (DM-mass) / (mass-splitting).");
  model->def_property_readonly(
      "w", [](const tsr::ModelMMS &mod) { return mod.w(); },
      "Scalar width scaled by mass-splitting: w = Gamma / (mass-splitting).");
  model->def_property(
      "g", [](const tsr::ModelMMS &mod) { return mod.g(); },
      [](tsr::ModelMMS &mod, double val) { mod.g(val); },
      "Scalar-DR-DM coupling.");
  model->def_property(
      "lam", [](const tsr::ModelMMS &mod) { return mod.lam(); },
      [](tsr::ModelMMS &mod, double val) { mod.lam(val); },
      "Number of scalars times the coupling squared: lam = g^2 * N. ");
}

auto add_tsr(py::class_<tsr::ModelMMS> *model) -> void {
  static std::string doc = R"Doc(
Compute the thermal dark-matter/dark-radiation scattering rate.

Parameters
----------
T: float
  Temperature of the thermal bath (dark-radiation temperature).

Returns
-------
tsr: float
  Thermal DM/DR scattering rate.
)Doc";
  static std::string doc_vec = R"Doc(
Compute the thermal dark-matter/dark-radiation scattering rate.

Parameters
----------
Ts: array-like
  Temperatures of the thermal bath (dark-radiation temperature).

Returns
-------
tsr: array-like
  Thermal DM/DR scattering rate.
)Doc";

  model->def("thermal_scattering_rate", &tsr::ModelMMS::thermal_scattering_rate,
             doc.c_str(), py::arg("T"));
  model->def(
      "thermal_scattering_rate",
      [](const tsr::ModelMMS &model, py::array_t<double> Ts) {
        py::buffer_info buf = Ts.request();
        auto result = py::array_t<double>(buf.size);
        py::buffer_info buf_res = result.request();

        double *ptrT = static_cast<double *>(buf.ptr);
        double *ptrR = static_cast<double *>(buf_res.ptr);

        for (size_t i = 0; i < buf.shape[0]; ++i) {
          ptrR[i] = model.thermal_scattering_rate(ptrT[i]);
        }
        return result;
      },
      doc_vec.c_str(), py::arg("Ts"));
}

auto add_tsr_integrand(py::class_<tsr::ModelMMS> *model) -> void {
  static std::string doc = R"Doc(
Compute the integrand of the thermal dark-matter/dark-radiation scattering
rate.

Parameters
----------
q: float
  Dark-radiation momentum scaled by dark-matter/mediator
  mass-splitting: q = p / delta.
T: flaot 
  Temperature of the thermal bath (dark-radiation temperature).

Returns
-------
tsr_integrand: float
  Integrand of the thermal DM/DR scattering rate.
)Doc";

  model->def("thermal_scattering_rate_integrand",
             &tsr::ModelMMS::thermal_scattering_rate_integrand, doc.c_str(),
             py::arg("q"), py::arg("T"));
}

auto add_msqrd(py::class_<tsr::ModelMMS> *model) -> void {
  static std::string doc = R"Doc(
Compute the t-averaged squared matrix element for DM + DR -> DM + DR.
rate.

Parameters
----------
q: float
  Dark-radiation momentum scaled by dark-matter/mediator
  mass-splitting: q = p / delta.

Returns
-------
msqrd: float
  t-averaged squared matrix element.
)Doc";

  model->def("msqrd_tavg", &tsr::ModelMMS::msqrd_tavg, doc.c_str(),
             py::arg("q"));
}

PYBIND11_MODULE(ts_rates, m) { // NOLINT

  auto mod = py::class_<tsr::ModelMMS>(m, "ModelMMS");
  add_ctor(&mod);
  add_properties(&mod);
  add_msqrd(&mod);
  add_tsr_integrand(&mod);
  add_tsr(&mod);
}
