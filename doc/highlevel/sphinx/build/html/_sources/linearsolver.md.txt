## linearSolver

| The linear solvers basically solves a linear-system of equations in  the form Ax=b. </br> In the case of Dumux it is Jx=r, with the Jacobian J and the residual r. </br> The used linear solvers are implemented in DUNE. |
| :--- |

### Key functionalities

- solve()
  - Input: Reference of A,x and b
  - Output: x

### Overview

```{mermaid}
flowchart LR
    A[newtonSolver] -->|"solve()"| B[linearSolver]
    click A "./newtonsolver.html"
    click B "./linearsolver.html"
```
