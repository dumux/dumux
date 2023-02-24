# Newton Loop

```{mermaid}
flowchart LR
    A(newtonSolver) -->|"assembleJacobianAndResidual()"| B(assembler)
    A -->|"solve()"| C(linearSolver)
    B -->|"assembleJacobianAndResidual()"| D(localAssembler)
    D -->|"evalFluxAndSourceResidual()"| E(localResidual)
    D -->|"evalStorageResidual()"| E
    click A "./newtonsolver.html"
    click B "./assembler.html"
    click C "./linearsolver.html"
    click D "./localassembler.html"
    click E "./localresidiual.html"
```