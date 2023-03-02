## Assembler {#assembler}
<!-- @page assembler Assembler -->

The assembler is used to calculate the global Jacobian-matrix and the global residual vector.
In order to do so, it uses a discretization specific localAssembler.
The assembler creates a shared_ptr of the Jacobian-matrix and a shared_ptr of the residual.

### Key functionalites

<!-- 1. assembleJacobianAndResidual() -->
<!--    1. Calls the assembleJacobianAndResidual() function of the @ref localAssembler for each element. -->

### Overview

@mermaid{assembler}
