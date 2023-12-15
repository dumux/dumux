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
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/common/parametertree.hh>

#include <dumux/linear/matrixconverter.hh>
#include <dumux/discretization/method.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/istlsolverfactorybackend.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearalgebratraits.hh>


namespace Dumux::Test{
struct MockGridGeometry
{
    using GridView = Dune::YaspGrid<2>::LeafGridView;
    using DofMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;
    using VertexMapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;
    using DiscretizationMethod = DiscretizationMethods::Box;
    static constexpr DiscretizationMethod discMethod{};
};

template<class M, class X, class V>
void solveWithFactory(M& A, X& x, V& b, const std::string& paramGroup)
{
    std::cout << std::endl;

    using LinearSolver = IstlSolverFactoryBackend<LinearSolverTraits<Test::MockGridGeometry>,
    LinearAlgebraTraits<M,V>>;
    LinearSolver solver(paramGroup);

    std::cout << "Solving Laplace problem with " << solver.name() << "\n";
    solver.solve(A, x, b);
    if (!solver.result().converged)
        DUNE_THROW(Dune::Exception, solver.name() << " did not converge!");
}

} // end namespace Dumux::Test

int main() {
    {
        using namespace Dune;
        {
            using namespace Indices;
            using BCRS_MAT= BCRSMatrix<FieldMatrix<double, 1, 1>>;
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
            using BCRS_ROW = Dune::MultiTypeBlockVector<BCRS_MAT, BCRS_MAT>;
            using MTBM = Dune::MultiTypeBlockMatrix<BCRS_ROW, BCRS_ROW>;
            MTBM M;
            M[_0][_0] = A;
            M[_0][_1] = B;
            M[_1][_0] = C;
            M[_1][_1] = D;

            // Create the multitype block vector
            using BV = Dune::BlockVector<FieldVector<double, 1>>;
            using MTBV = Dune::MultiTypeBlockVector<BV, BV>;
            BV R0;
            BV R1;
            loadMatrixMarket(R0, "data/rhs_0.log");
            loadMatrixMarket(R1, "data/rhs_1.log");


            MTBV R;
            R[_0] = R0;
            R[_1] = R1;

            auto r0n = R0.N();
            auto r1n = R1.N();
            MTBV X;
            X[_0].resize(r0n);
            X[_1].resize(r1n);

            BCRS_MAT bcrsM = Dumux::MatrixConverter<decltype(M)>::multiTypeToBCRSMatrix(M);
            BV blockVectorR = Dumux::VectorConverter<decltype(R)>::multiTypeToBlockVector(R);
            BV blockVectorX = Dumux::VectorConverter<decltype(X)>::multiTypeToBlockVector(X);
            {
                using namespace Dumux;
                Dune::ParameterTree params;
                params["maxit"] = "250";
                params["reduction"] = "1e-13";
                params["verbose"]="2";

                using LinearSolver = ILUBiCGSTABIstlSolver<LinearSolverTraits<Test::MockGridGeometry>,
                        LinearAlgebraTraits<BCRS_MAT ,BV>>;

                LinearSolver solver(params);
                std::cout<<"Solving Real Matrix Problem with:"<<"\n";
                auto r = solver.solve(bcrsM, blockVectorX, blockVectorR);
                std::cout<<"Converged: "<<r.converged<<"\n";


            }
        }
    }
    return 0;
}
