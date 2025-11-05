# Printing system matrix

Debugging may require inspecting the system (Jacobian) matrix from time to time.

The following example shows a possible way to achieve this for a staggered-grid free-flow problem:

We add the output within the `solveLinearSystem_` overload for `Dune::MultiTypeBlockMatrix` Jacobians in `newtonsolver.hh`:

```cpp
    //  MultiTypeBlockMatrix is not supported for printing, needs to be converted first

    const auto& A = this->assembler().jacobian();
    const auto& b = this->assembler().residual();

    // create the bcrs matrix the IterativeSolver backend can handle
    const auto M = MatrixConverter<JacobianMatrix>::multiTypeToBCRSMatrix(A);

    // get the new matrix sizes
    const std::size_t numRows = M.N();
    assert(numRows == M.M());

    // create the vector the IterativeSolver backend can handle
    const auto bTmp = VectorConverter<SolutionVector>::multiTypeToBlockVector(b);
    assert(bTmp.size() == numRows);

    // create a blockvector to which the linear solver writes the solution
    using VectorBlock = typename Dune::FieldVector<Scalar, 1>;
    using BlockVector = typename Dune::BlockVector<VectorBlock>;
    BlockVector y(numRows);

    /////// BEGIN PRINTING /////////////////////////////////////////////////////
    static const bool printmatrix = getParam<bool>("Problem.PrintMatrix", false);

    if (printmatrix)
    {
        static int counter = 0;
        const auto rank = Dune::MPIHelper::getCommunication().rank();

        Dune::storeMatrixMarket(M, "matrix_" + std::to_string(rank) +  "_iter_" + std::to_string(counter) + ".log");
        Dune::storeMatrixMarket(bTmp, "rhs_" + std::to_string(rank) +  "_iter_" + std::to_string(counter) + ".log");

        ++counter;
    }
  /////// END PRINTING /////////////////////////////////////////////////////

// ....
```

This will write out files for the matrix and the residual vector for each Newton iteration and each process. `A` and `b` can also be accessed via `assembler.jacobian()` and `assembler.residual()` respectively.

You may then compare files using, e.g., `kompare` or `meld`.

Keep in mind that the `MatrixMarket`format uses indices starting with `1` instead of `0`.

For rather small matrices, you may also directly print the output to the terminal, giving you a graphical impression of the structure:

```cpp
Dune::printmatrix(std::cout, M, "", "", 10/*width*/, 2/*precision*/);
```
