# Part 1: Main program flow

The main program flow is implemented in file `main.cc`, described below. Unlike most
DuMu<sup>x</sup> examples, it does not assemble and solve a PDE -- `TwoPStatic` performs a
purely topological invasion-percolation sweep -- so `main()` is a straight-line driver: it
sets up a pore-network grid, precomputes the throat entry and snap-off capillary pressures
once, then steps a global capillary pressure up (drainage) and back down (imbibition),
calling `Dumux::PoreNetwork::TwoPStatic::updateInvasionState(...)` and
`updateTrappedState(...)` at every step while tracking per-pore saturations to build the
$`p_c`$-$`S_w`$ curve. For the underlying equations and the algorithm in pseudocode, see
[Mathematical and numerical model](../README.md#mathematical-and-numerical-model) in the
main description.

The code documentation is structured as follows:

[[_TOC_]]
