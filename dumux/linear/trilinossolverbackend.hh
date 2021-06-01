// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup Linear
 * \brief Trilinos solver backend
 */
#ifndef DUMUX_TRILINOS_SOLVER_BACKEND_HH
#define DUMUX_TRILINOS_SOLVER_BACKEND_HH

//#define TRILINOS_SOLVER_BACKEND_TEST_OUTPUT

#include <dune/istl/solver.hh>

#include <dumux/common/typetraits/matrix.hh>
#include <dumux/linear/solver.hh>

#if HAVE_TRILINOS
#include <ShyLU_DDFROSch_config.h>

#include <Xpetra_ConfigDefs.hpp>
#include <Xpetra_DefaultPlatform.hpp>

#include <Xpetra_MatrixFactory.hpp>
#include <Xpetra_CrsGraphFactory.hpp>
#include <Xpetra_MultiVectorFactory.hpp>
#include <Xpetra_CrsMatrixWrap.hpp>
#include "Xpetra_CrsMatrix.hpp"
#include <Xpetra_MapFactory.hpp>

#include <FROSch_Tools_def.hpp>

#include <Stratimikos_FROSch_def.hpp>

#include <Thyra_LinearOpWithSolveBase.hpp>
#include <Thyra_LinearOpWithSolveFactoryHelpers.hpp>
#endif

namespace Dumux {

#if HAVE_TRILINOS
/*!
 * \ingroup Linear
 * \brief Linear solver using Trilinos.
 */
class TrilinosSolverBackend : public LinearSolver
{
public:
    using LinearSolver::LinearSolver;

    template<class Matrix, class Vector>
    bool solve(const Matrix& A, Vector& x, const Vector& b)
    {
        using namespace FROSch;
        using namespace std;
        using namespace Teuchos;
        using namespace Thyra;
        using namespace Xpetra;
        using SC = double;
        using LO = int;
        using GO = DefaultGlobalOrdinal;
        // using GO = long long;
        using NO = KokkosClassic::DefaultNode::DefaultNodeType;
        using mapPtr = RCP<Map<LO,GO,NO> >;
        using mapFactory = MapFactory<LO,GO,NO>;
        using vectorPtr = RCP<MultiVector<SC,LO,GO,NO> >;
        using vectorFactory = MultiVectorFactory<SC,LO,GO,NO>;
        using matrixPtr = RCP<Xpetra::Matrix<SC,LO,GO,NO> >;
        using matrixFactory = MatrixFactory<SC,LO,GO,NO>;

        static_assert(isBCRSMatrix<Matrix>::value, "SuperLU only works with BCRS matrices!");
        using BlockType = typename Matrix::block_type;
        static_assert(BlockType::rows == BlockType::cols, "Matrix block must be quadratic!");
        constexpr auto blockSize = BlockType::rows;

        RCP<const Comm<int> > comm = DefaultPlatform::getDefaultPlatform().getComm();

        global_size_t numRows = x.size()*blockSize;
        size_t maxNumEntriesPerRow = 10; //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        mapPtr mapXpetra = mapFactory::Build(UseTpetra,numRows,0,comm);
        vectorPtr xXpetra = vectorFactory::Build(mapXpetra,1);
        vectorPtr bXpetra = vectorFactory::Build(mapXpetra,1);
        matrixPtr AXpetra = matrixFactory::Build(mapXpetra,maxNumEntriesPerRow);

        for (auto blockIdx = 0u; blockIdx < x.size(); ++blockIdx)
        {
#ifdef TRILINOS_SOLVER_BACKEND_TEST_OUTPUT
            cout << "block idx " << blockIdx << ": ";

            cout << "x = ";
#endif
            for (auto i = 0u; i < blockSize; ++i) {
                xXpetra->replaceGlobalValue(blockIdx*blockSize+i,0,x[blockIdx][i]);
#ifdef TRILINOS_SOLVER_BACKEND_TEST_OUTPUT
                cout << x[blockIdx][i] << ", ";
#endif
            }

#ifdef TRILINOS_SOLVER_BACKEND_TEST_OUTPUT
            cout << "b = ";
#endif
            for (auto i = 0u; i < blockSize; ++i) {
                bXpetra->replaceGlobalValue(blockIdx*blockSize+i,0,b[blockIdx][i]);
#ifdef TRILINOS_SOLVER_BACKEND_TEST_OUTPUT
                cout << b[blockIdx][i] << ", ";
#endif
            }

#ifdef TRILINOS_SOLVER_BACKEND_TEST_OUTPUT
            cout << endl;
#endif
        }

        // Tpetra::CrsMatrix At;

        for (auto rowIt = A.begin(); rowIt != A.end(); ++rowIt)
        {
#ifdef TRILINOS_SOLVER_BACKEND_TEST_OUTPUT
            cout << "block row idx " << rowIt.index() << ":" << endl;
#endif

            for (auto colIt = rowIt->begin(); colIt != rowIt->end(); ++colIt)
            {
#ifdef TRILINOS_SOLVER_BACKEND_TEST_OUTPUT
                cout << "  block col idx " << colIt.index() << ": block entry = { ";
#endif

                for (auto i = 0u; i < blockSize; ++i)
                {
#ifdef TRILINOS_SOLVER_BACKEND_TEST_OUTPUT
                    cout << "{ ";
#endif

                    Array<GO> cols(blockSize);
                    Array<SC> vals(blockSize);
                    for (auto j = 0u; j < blockSize; ++j)
                    {
                        cols[j] = colIt.index()*blockSize+j;
                        vals[j] = (*colIt)[i][j];
#ifdef TRILINOS_SOLVER_BACKEND_TEST_OUTPUT
                        cout << (*colIt)[i][j] << ", ";
#endif
                    }
                    AXpetra->insertGlobalValues(rowIt.index()*blockSize+i,cols(),vals());

#ifdef TRILINOS_SOLVER_BACKEND_TEST_OUTPUT
                    cout << " }, ";
#endif
                }

#ifdef TRILINOS_SOLVER_BACKEND_TEST_OUTPUT
                cout << " }" << endl;
#endif
            }
        }
        AXpetra->fillComplete();
#ifdef TRILINOS_SOLVER_BACKEND_TEST_OUTPUT
        RCP<FancyOStream> fancy = fancyOStream(rcpFromRef(cout));
        xXpetra->describe(*fancy,VERB_EXTREME);
        bXpetra->describe(*fancy,VERB_EXTREME);
        AXpetra->describe(*fancy,VERB_EXTREME);
#endif

        CrsMatrixWrap<SC,LO,GO,NO>& crsWrapA = dynamic_cast<CrsMatrixWrap<SC,LO,GO,NO>&>(*AXpetra);
        RCP<const LinearOpBase<SC> > AThyra = ThyraUtils<SC,LO,GO,NO>::toThyra(crsWrapA.getCrsMatrix());
        RCP<MultiVectorBase<SC> > xThyra = rcp_const_cast<MultiVectorBase<SC> >(ThyraUtils<SC,LO,GO,NO>::toThyraMultiVector(xXpetra));
        RCP<const MultiVectorBase<SC> > bThyra = ThyraUtils<SC,LO,GO,NO>::toThyraMultiVector(bXpetra);

        // Build the solver
        RCP<ParameterList> parameterList = getParametersFromXmlFile("ParameterList_Thyra.xml");
        DefaultLinearSolverBuilder linearSolverBuilder;
        enableFROSch<int,int,KokkosClassic::DefaultNode::DefaultNodeType>(linearSolverBuilder);
        linearSolverBuilder.setParameterList(parameterList);

        RCP<LinearOpWithSolveFactoryBase<double> > lowsFactory = linearSolverBuilder.createLinearSolveStrategy("");
        RCP<FancyOStream> out = VerboseObjectBase::getDefaultOStream();
        lowsFactory->setOStream(out);
        lowsFactory->setVerbLevel(VERB_HIGH);
        RCP<LinearOpWithSolveBase<double> > lows = linearOpWithSolve(*lowsFactory, AThyra);

        // Solve
        SolveStatus<double> status = Thyra::solve<double>(*lows, Thyra::NOTRANS, *bThyra, xThyra.ptr());
        xXpetra = ThyraUtils<SC,LO,GO,NO>::toXpetra(xThyra,comm);

#ifdef TRILINOS_SOLVER_BACKEND_TEST_OUTPUT
        xXpetra->describe(*fancy,VERB_EXTREME);
        bXpetra->describe(*fancy,VERB_EXTREME);
#endif

        ArrayRCP<const SC> vals = xXpetra->getData(0);
        for (auto blockIdx = 0u; blockIdx < x.size(); ++blockIdx)
        {
#ifdef TRILINOS_SOLVER_BACKEND_TEST_OUTPUT
            cout << "block idx " << blockIdx << ": ";

            cout << "x = ";
#endif
            for (auto i = 0u; i < blockSize; ++i) {
                x[blockIdx][i] = vals[blockIdx*blockSize+i];
#ifdef TRILINOS_SOLVER_BACKEND_TEST_OUTPUT
                cout << x[blockIdx][i] << ", ";
#endif
            }

#ifdef TRILINOS_SOLVER_BACKEND_TEST_OUTPUT
            cout << endl;
#endif
        }

        result_.converged = true;

        return result_.converged;
    }

    std::string name() const
    {
        return "Trilinos solver";
    }

    const Dune::InverseOperatorResult& result() const
    {
        return result_;
    }

private:
    Dune::InverseOperatorResult result_;
};
#endif // HAVE_TRILINOS

} // end namespace Dumux

#endif
