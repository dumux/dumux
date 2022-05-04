# Richards test on annulus

This may corresponds to a rhizosphere model, where we have water uptake
by a root on the inner boundary. By transforming the Richards equation
using a Kirchhoff transformation, the flux term becomes linear and
we can compute an analytic solution for the stationary equation.
Under some assumptions on the storage term, we can also compute
approximative solutions to the instationary equation. Here, we
test convergence against the analytical solution by implementing
storage approximation as volumetric source terms (possibly nonlinear).

The benchmark case is described in https://doi.org/10.3389/fpls.2020.00316
See problem documentation for the analytic solution and how it's derived.

# Tests

There are two (four) tests implemented.

* `test_richards_annulus` is a convergence
test against the stationary analytical solution when the root-soil interface
pressure reached the wilting point pressure (-15000 cm pressure head).
The test runs for loam and clay and makes sure the convergence rate is better than 1.95.

* `test_richards_annulus_instationary` is a time-dependent test. At each time step,
the analytical solution corresponding to the steady-rate approximation with root-soil
pressure given by the numerical solution is computed. This solution corresponds to
a runtime (since the uptake rate is constant). The numerical runtime and the analytically
approximated runtime are compared. The analytic approximation is very accurate and the test
makes sure the difference is below 5%.
