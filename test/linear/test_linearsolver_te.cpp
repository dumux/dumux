#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>
#include <dune/istl/bvector.hh>
#include <dune/istl/multitypeblockvector.hh>
#include <dune/istl/multitypeblockmatrix.hh>
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/matrixmarket.hh>
#include <dune/common/indices.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/preconditioners.hh>
#include <dune/istl/umfpack.hh>

#include <dumux/linear/matrixconverter.hh>

using namespace Dune;
using namespace Indices;


int main() {

    typedef BCRSMatrix<FieldMatrix<double,1,1>> BCRS_MAT;
    // Create the single bcrs matrix
    BCRS_MAT A;
    BCRS_MAT B;
    BCRS_MAT C;
    BCRS_MAT D;

    //Use load matrix market to load the matrix
    loadMatrixMarket(A, "data/matrix_0_0.log");
    loadMatrixMarket(B, "data/matrix_0_1.log");
    loadMatrixMarket(C, "data/matrix_1_0.log");
    loadMatrixMarket(D, "data/matrix_1_1.log");

    // Create the multitype block matrix from the single bcrs matrices
    typedef Dune::MultiTypeBlockVector<BCRS_MAT,BCRS_MAT> BCRS_ROW;
    typedef Dune::MultiTypeBlockMatrix<BCRS_ROW,BCRS_ROW> MTBM;
    MTBM M;
    M[_0][_0] = A;
    M[_0][_1] = B;
    M[_1][_0] = C;
    M[_1][_1] = D;

    // Create the multitype block vector
    typedef Dune::BlockVector<FieldVector<double,1>> BV;
    typedef Dune::MultiTypeBlockVector<BV,BV> MTBV;
    BV R0;
    BV R1;
    loadMatrixMarket(R0, "data/rhs_0.log");
    loadMatrixMarket(R1, "data/rhs_1.log");


    MTBV R;
    R[_0] = R0;
    R[_1] = R1;

    auto r0n= R0.N();
    auto r1n=R1.N();
    MTBV X;
    X[_0].resize(r0n);
    X[_1].resize(r1n);

    InverseOperatorResult r;
    BCRS_MAT bcrsM = Dumux::MatrixConverter<decltype(M)>::multiTypeToBCRSMatrix(M);
    BV blockVectorR = Dumux::VectorConverter<decltype(R)>::multiTypeToBlockVector(R);
    BV blockVectorX = Dumux::VectorConverter<decltype(X)>::multiTypeToBlockVector(X);


    MatrixAdapter<BCRS_MAT ,BV,BV> op(bcrsM);
    Dune::UMFPack<BCRS_MAT> solver(bcrsM,1);
    solver.apply(blockVectorX,blockVectorR,r);








    return 0;
}
