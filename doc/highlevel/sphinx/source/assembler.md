# Assembler
The assembler calculates the global Jacobian-matrix and the global residual vector.
In order to do so, it uses a discretization specific localAssembler.
Once the assembler is finished, a shared_ptr of the Jacobian-matrix and a shared_ptr of the residual is created.
A linear solver can then calculate a new solution.