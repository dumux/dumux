// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Linear
 * \brief Linear solvers from Trilinos
 * This solver backend requires MPI (even for sequential runs)
 */
#ifndef DUMUX_LINEAR_TRILINOS_SOLVERS_HH
#define DUMUX_LINEAR_TRILINOS_SOLVERS_HH

#include <string>
#include <dumux/assembly/jacobianpattern.hh>


#if DUMUX_HAVE_TRILINOS

#include <Teuchos_ScalarTraits.hpp>
#include <Teuchos_RCP.hpp>
#include <Teuchos_Tuple.hpp>
#include <Teuchos_VerboseObject.hpp>
#include <Teuchos_ParameterList.hpp>

#include <Tpetra_Core.hpp>
#include <Tpetra_Map.hpp>
#include <Tpetra_CrsMatrix.hpp>
#include <Tpetra_MultiVector.hpp>

#include <Amesos2.hpp>
#include <Amesos2_Version.hpp>
#include <Amesos2_MatrixAdapter_decl.hpp>

#define HAVE_MPI 1
#if HAVE_MPI
#include <mpi.h>
#endif

namespace Dumux::Detail {

struct TrilinosSolverResult : public Dune::InverseOperatorResult
{
    TrilinosSolverResult() = default;
    TrilinosSolverResult(const TrilinosSolverResult&) = default;
    TrilinosSolverResult(TrilinosSolverResult&&) = default;

    TrilinosSolverResult(const Dune::InverseOperatorResult& o) : InverseOperatorResult(o) {}
    TrilinosSolverResult(Dune::InverseOperatorResult&& o) : InverseOperatorResult(std::move(o)) {}

    operator bool() const { return this->converged; }
};

namespace TrilinosImpl {

template<class Graph, class GridGeometry, class CommPtr>
requires (GridGeometry::discMethod == DiscretizationMethods::box)
auto sparsityGraph(const GridGeometry& gridGeometry, CommPtr comm)
{
    using NodeType = typename Graph::node_type;
    using LocalIDType = typename Graph::local_ordinal_type;
    using GlobalIDType = typename Graph::global_ordinal_type;
    using Map = Tpetra::Map<LocalIDType, GlobalIDType, NodeType>;

    const auto rank = gridGeometry.gridView().comm().rank();
    const auto& mapper = gridGeometry.dofMapper();
    std::vector<std::size_t> isOwned(mapper.size(), rank); // row owned by this process



    const auto numRowsGlobal_ = gridGeometry.dofMapper().size();
    const auto rowMap_ = Teuchos::rcp(new Map (numRowsGlobal_, GlobalIDType{0}, comm));

    for (const auto& element : elements(gridGeometry.gridView(), Dune::Partitions::interior))
    {

    }

    Teuchos::ArrayRCP<std::size_t> numEntPerRow(numRows);
    // std::cout << "Dumux pattern: numRows: " << numRows << std::endl;
    for (std::size_t rowIdx = 0; rowIdx < numRows; ++rowIdx)
    {
        numEntPerRow[rowIdx] = pattern.rowsize(rowIdx);
        // std::cout << "Dumux pattern: rowIdx: " << rowIdx << ", numEntPerRow: " << numEntPerRow[rowIdx] << std::endl;
    }

    const auto colMap_ = rowMap_;
    auto graph = Teuchos::rcp(new Graph(rowMap_, colMap_, numEntPerRow()));
    for (std::size_t rowIdx = 0; rowIdx < numRows; ++rowIdx)
    {
        Teuchos::Array<LocalIDType> colIndices(
            std::visit([](const auto& indices) {
                return Teuchos::Array<LocalIDType>(indices.begin(), indices.end());
            }, pattern.columnIndices(rowIdx))
        );

        graph->insertLocalIndices(rowIdx, colIndices());
    }

    graph->fillComplete();

    //graph->describe(*fos_, Teuchos::VERB_EXTREME);
    return graph;
}

} // end namespace Tpetra

} // end namespace Dumux::Detail

namespace Dumux {

/*!
 * \ingroup Linear
 * \brief Direct linear solvers from Trilinos Amesos2
 *
 * Amesos2 is a package that provides direct solvers for large sparse linear systems.
 * and is part of the Trilinos project. It provides interface to a number of direct solvers.
 * It is based on Tpetra::CrsMatrix and Tpetra::MultiVector. As usual with Trilinos,
 * Teuchos is used as a utility library for parameter handling.
 */
template<class LSTraits, class LATraits>
class DirectSolverAmesos2 : public LinearSolver
{
    using Matrix = typename LATraits::Matrix;
    using XVector = typename LATraits::Vector;
    using BVector = typename LATraits::Vector;

    using Scalar = double;
    using Magnitude = Teuchos::ScalarTraits<Scalar>::magnitudeType;

    using NodeType = Tpetra::KokkosClassic::DefaultNode::DefaultNodeType;
    using LocalIDType = Tpetra::Map<>::local_ordinal_type;
    using GlobalIDType = Tpetra::Map<>::global_ordinal_type;
    using Map = Tpetra::Map<LocalIDType, GlobalIDType, NodeType>;

    using SolverMatrix = Tpetra::CrsMatrix<Scalar, LocalIDType, GlobalIDType, NodeType>;
    using SolverVector = Tpetra::MultiVector<Scalar, LocalIDType, GlobalIDType, NodeType>;
    using Solver = Amesos2::Solver<SolverMatrix, SolverVector>;
    using Graph = Tpetra::CrsGraph<LocalIDType, GlobalIDType, NodeType>;

    typedef Kokkos::DefaultHostExecutionSpace                     HostExecSpaceType;
    typedef typename HostExecSpaceType::memory_space              HostMemSpaceType;
    typedef Kokkos::View<LocalIDType*, HostExecSpaceType>  host_ordinal_type_view;
    using mumps_type = Amesos2::TypeMap<Amesos2::MUMPS, Scalar>::type;
    typedef Kokkos::View<mumps_type*, HostExecSpaceType>          host_value_type_view;

    using LinearSolverTraits = LSTraits;

public:
    using SolverResult = Detail::TrilinosSolverResult;

    template<class GridGeometry>
    DirectSolverAmesos2(const GridGeometry& gridGeometry,
                        const typename LinearSolverTraits::GridView& gridView,
                        const typename LinearSolverTraits::DofMapper& dofMapper,
                        const std::string& paramGroup = "")
    : LinearSolver(paramGroup)
#if HAVE_MPI
    , gridView_(gridView)
    , mapper_(dofMapper)
    , globalIdSet_(gridView.grid().globalIdSet())
    , isParallel_(gridView.comm().size() > 1)
#endif
    {
        params_ = Teuchos::ParameterList(
            "Amesos2"
            //getParamFromGroup<std::string>(paramGroup, "LinearSolver.ConfigFilename")
        );

        comm_ = Teuchos::rcp(new Teuchos::MpiComm<int>(static_cast<MPI_Comm>(gridView.comm())));
        fos_ = Teuchos::fancyOStream(Teuchos::rcpFromRef(std::cout));
        rank_ = comm_->getRank();
        if (rank_ == 0)
            *fos_ << "\n" << Amesos2::version() << std::endl;

        // build sequential communication infrastructure
        graph_ = Detail::TrilinosImpl::sparsityGraph<Graph>(gridGeometry, comm_);
    }

    /*!
     * \brief Solve the linear system Ax = b
     */
    SolverResult solve(const Matrix& A, XVector& x, const BVector& b)
    {
        return solve_(A, x, b);
    }

    Scalar norm(const XVector& x) const
    {
        return 0.0;
    }

    /*!
     * \brief name of the linear solver
     */
    std::string name() const
    {
        return "Trilinos Amesos2 Direct Solver";
    }

private:
    SolverResult solve_(const Matrix& A, XVector& x, const BVector& b)
    {
        SolverResult result;
        auto AA = convertMatrix_(A);
        auto BB = convertBVector_(b);
        auto XX = convertXVector_(x);

        mumps_par.comm_fortran = -987654;
        Teuchos::RCP<const Teuchos::Comm<int> > matComm = AA->getComm();

        if (matComm.is_null())
            throw std::logic_error("Amesos2::Comm");

        Teuchos::RCP<const Teuchos::MpiComm<int> > matMpiComm
            = Teuchos::rcp_dynamic_cast<const Teuchos::MpiComm<int> >(matComm);

        if (matMpiComm.is_null()) // required for Amesos2 backend
            throw std::logic_error("Amesos2::MPI The matrix's communicator is not an MpiComm");

        MPI_Comm rawMpiComm = (* (matMpiComm->getRawMpiComm()) )();
        mumps_par.comm_fortran = (int) MPI_Comm_c2f(rawMpiComm);

        std::cout << "MUMPS: par comm_fortran: " << mumps_par.comm_fortran << std::endl;

        Teuchos::RCP<Solver> solver = Amesos2::create<SolverMatrix, SolverVector>("MUMPS", AA, XX, BB);
        solver->setParameters( Teuchos::rcpFromRef(params_) );
        solver->symbolicFactorization().numericFactorization().solve();
        retrieveXVector_(XX, x);

        result.converged = true;
        return result;
    }

    auto convertXVector_(const XVector& x)
    {
        Teuchos::RCP<SolverVector> XX = Teuchos::rcp(new SolverVector(rowMap_, 1, true));
        return XX;
    }

    void retrieveXVector_(const Teuchos::RCP<SolverVector>& XX, XVector& x)
    {
        const auto view = XX->getLocalViewHost(Tpetra::Access::ReadOnly);
        const auto subview = Kokkos::subview(view, Kokkos::ALL(), 0);
        constexpr std::size_t blockSize = XVector::block_type::size();
        for (auto blockIdx = 0u; blockIdx < x.size(); ++blockIdx)
            for (auto j = 0u; j < blockSize; ++j)
                x[blockIdx][j] = subview(blockIdx*blockSize + j);
    }

    auto convertBVector_(const BVector& b)
    {
        Teuchos::RCP<SolverVector> BB = Teuchos::rcp(new SolverVector(rowMap_, 1, true));
        constexpr std::size_t blockSize = BVector::block_type::size();
        for (auto blockIdx = 0u; blockIdx < b.size(); ++blockIdx)
            for (auto j = 0u; j < blockSize; ++j)
                BB->replaceLocalValue(blockIdx*blockSize + j, 0, b[blockIdx][j]);
        return BB;
    }

    auto convertMatrix_(const Matrix& A)
    {
        // set sparsity pattern via the graph
        Teuchos::RCP<SolverMatrix> AA = Teuchos::rcp( new SolverMatrix(graph_) );

        using BlockType = typename Matrix::block_type;
        const auto blockSizeI = BlockType::rows;
        const auto blockSizeJ = BlockType::cols;

        // fill matrix entries
        for (auto rowIt = A.begin(); rowIt != A.end(); ++rowIt)
        {
            for (int i = 0; i < blockSizeI; ++i)
            {
                for (auto colIt = rowIt->begin(); colIt != rowIt->end(); ++colIt)
                {
                    Teuchos::Array<GlobalIDType> cols(blockSizeJ);
                    Teuchos::Array<Scalar> vals(blockSizeJ);
                    for (int j = 0; j < blockSizeJ; ++j)
                    {
                        cols[j] = colIt.index() * blockSizeJ + j;
                        vals[j] = (*colIt)[i][j];
                    }

                    AA->replaceLocalValues(rowIt.index() * blockSizeI + i, cols(), vals());
                }
            }
        }

        AA->fillComplete();
        return AA;
    }

#if HAVE_MPI
    const LinearSolverTraits::GridView gridView_;
    const LinearSolverTraits::DofMapper& mapper_;
    const LinearSolverTraits::GridView::Traits::Grid::Traits::GlobalIdSet& globalIdSet_;
    bool isParallel_;

    Teuchos::RCP<const Teuchos::Comm<int>> comm_;
    Teuchos::RCP<Teuchos::FancyOStream> fos_;
    Teuchos::RCP<const Map> rowMap_;
    Teuchos::RCP<const Map> colMap_;
    Teuchos::RCP<const Graph> graph_;
#endif

    mutable Amesos2::TypeMap<Amesos2::MUMPS, Scalar>::MUMPS_STRUC_C mumps_par;

    std::size_t rank_;
    std::vector<int> isOwned_;
    std::size_t numRowsGlobal_;

    host_ordinal_type_view host_rows_view_;
    host_ordinal_type_view host_col_ptr_view_;
    host_value_type_view host_nzvals_view_;

    Teuchos::ParameterList params_;
};

} // end namespace Dumux

#endif // DUMUX_HAVE_TRILINOS

#endif
