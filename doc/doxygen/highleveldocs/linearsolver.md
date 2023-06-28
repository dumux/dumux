@page linearsolver LinearSolver

The linear solvers solves a linear-system of equations in  the form Ax=b.
In the case of Dumux it is Jx=r, with the Jacobian J and the residual r.
The used linear solvers are implemented in DUNE.


### Key functionalities

- solve()
  - solves a system of equation of the Form Jx=r

### Overview

@mermaid{linearsolver}
