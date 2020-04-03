# Part 3: Main program flow

We want to solve a single-phase flow problem
to obtain a pressure distribution in the domain. Subsequently, we compute a volume fluxes field
based on the obtained pressure distribution and pass it to the tracer problem.
Finally, we solve a transient transport problem for a tracer using the computed volume fluxes.
This main program flow is implemented in the `main()` function
of the program which is defined in the file `main.cc` described below.

The code documentation is structured as follows:

[[_TOC_]]
