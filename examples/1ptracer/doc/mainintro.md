# Part 3: Main program flow

This file contains the main program flow. In this example, we solve a single-phase flow problem
to obtain a pressure distribution in the domain. Subsequently, the distribution of volume fluxes
is computed from that pressure distribution, which is then passed to a tracer problem to solve
the transport of an initial contamination through the model domain.

The subsequent code documentation is structured as follows:

[[_TOC_]]
