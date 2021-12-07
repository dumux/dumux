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

#define TRILINOS_SOLVER_BACKEND_TEST_OUTPUT

#include <dune/istl/solver.hh>

#include <dumux/common/typetraits/matrix.hh>
#include <dumux/linear/solver.hh>

#include <dumux/parallel/vectorcommdatahandle.hh>

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
template <class LinearSolverTraits>
class TrilinosSolverBackend : public LinearSolver
{
    using DuneComm = Dune::OwnerOverlapCopyCommunication<Dune::bigunsignedint<96>, int>;
    using GlobalIdSet = typename LinearSolverTraits::GridView::Grid::GlobalIdSet;
    using GridView = typename LinearSolverTraits::GridView;
    using Mapper = typename LinearSolverTraits::DofMapper;

public:
    TrilinosSolverBackend(const typename LinearSolverTraits::GridView& gridView,
                          const typename LinearSolverTraits::DofMapper& dofMapper,
                          const std::string& paramGroup = "")
    : LinearSolver(paramGroup)
#if HAVE_MPI
    , gridView_(gridView)
    , mapper_(dofMapper)
    , globalIdSet_(gridView.grid().globalIdSet())
    , isParallel_(Dune::MPIHelper::getCollectiveCommunication().size() > 1)
    , duneToTrilinosMap_(gridView.size(0),-1)
    , duneToTrilinosMapOverlap_(gridView.size(0),-1)
#endif
    {}

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

        // Construct the global indices
        size_t numInteriorRows = 0;
        for (const auto& element : elements(gridView_)) {
            if (element.partitionType() == Dune::InteriorEntity) numInteriorRows++;
        }
        mapPtr mapXpetra = mapFactory::Build(UseTpetra,-1,numInteriorRows,0,comm);

        // Build overlapping map and the Dune-to-Trilinos mappings
        std::vector<int> globalIDs(gridView_.size(0),-1);
        int indAll = 0, indInt = 0;
        for (const auto& element : elements(gridView_))
        {
            if (element.partitionType() == Dune::InteriorEntity)
            {
                globalIDs[indAll] = mapXpetra->getGlobalElement(indInt);
                duneToTrilinosMap_.at(mapper_.index(element)) = indInt;
                indInt++;
            }
            duneToTrilinosMapOverlap_.at(mapper_.index(element)) = indAll;
            indAll++;
        }
        VectorCommDataHandleMax<Mapper,std::vector<int>,0> maxHandle(mapper_,globalIDs);
        gridView_.communicate(maxHandle,Dune::InteriorBorder_All_Interface,Dune::ForwardCommunication);

        Teuchos::ArrayView<int> globalIDsOverlap(globalIDs);
        mapPtr mapXpetraOverlap = mapFactory::Build(UseTpetra,-1,globalIDsOverlap,0,comm);

        // Build matrix and vectors
        size_t maxNumEntriesPerRow = 10; // Can we estimate this?
        matrixPtr AXpetra = matrixFactory::Build(mapXpetra,mapXpetraOverlap,maxNumEntriesPerRow);

        vectorPtr xXpetra = vectorFactory::Build(mapXpetra,1);
        vectorPtr bXpetra = vectorFactory::Build(mapXpetra,1);

        // Fill vectors
        for (auto blockIdx = 0u; blockIdx < x.size(); ++blockIdx) {
            // x
            for (auto i = 0u; i < blockSize; ++i) {
                int rowidx = duneToTrilinosMap_.at(blockIdx*blockSize+i);
                if (rowidx!=-1) xXpetra->replaceLocalValue(rowidx,0,x[blockIdx][i]);
            }

            // b
            for (auto i = 0u; i < blockSize; ++i) {
                int rowidx = duneToTrilinosMap_.at(blockIdx*blockSize+i);
                if (rowidx!=-1) bXpetra->replaceLocalValue(rowidx,0,b[blockIdx][i]);
            }
        }

        // Fill matrix
        for (auto rowIt = A.begin(); rowIt != A.end(); ++rowIt) {
            int rowidx = duneToTrilinosMap_.at(rowIt.index()*blockSize); // only for blocksize = 1
            if (rowidx!=-1) {
                for (auto colIt = rowIt->begin(); colIt != rowIt->end(); ++colIt) {
                    int colidx = duneToTrilinosMapOverlap_.at(colIt.index()*blockSize); // only for blocksize = 1
                    if (colidx!=-1) {
                        for (auto i = 0u; i < blockSize; ++i) {
                            Array<GO> cols(blockSize);
                            Array<SC> vals(blockSize);
                            for (auto j = 0u; j < blockSize; ++j) {
                                cols[j] = colidx;
                                vals[j] = (*colIt)[i][j];
                            }
                            AXpetra->insertLocalValues(rowidx,cols(),vals());
                        }
                    }
                }
            }
        }
        AXpetra->fillComplete(mapXpetra,mapXpetra);

        // Thyra wrappers
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

        // write the solution back into x
        ArrayRCP<const SC> vals = xXpetra->getData(0);
        for (auto blockIdx = 0u; blockIdx < x.size(); ++blockIdx) {
            for (auto i = 0u; i < blockSize; ++i) {
                int rowidx = duneToTrilinosMap_.at(blockIdx*blockSize);
                if (rowidx!=-1) x[blockIdx][i] = vals[rowidx];
            }
        }

        result_.converged = true; // this should depend on status

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
    const GridView gridView_;
    const Mapper& mapper_;
    const GlobalIdSet& globalIdSet_;
    bool isParallel_;
    std::vector<int> duneToTrilinosMap_;
    std::vector<int> duneToTrilinosMapOverlap_;
};
#endif // HAVE_TRILINOS

} // end namespace Dumux

#endif
