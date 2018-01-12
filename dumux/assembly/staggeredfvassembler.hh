// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 * \ingroup Assembly
 * \ingroup StaggeredDiscretization
 * \brief A linear system assembler (residual and Jacobian) for staggered finite volume schemes
 */
#ifndef DUMUX_STAGGERED_FV_ASSEMBLER_HH
#define DUMUX_STAGGERED_FV_ASSEMBLER_HH

#include <type_traits>

#include <dune/istl/matrixindexset.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/timeloop.hh>
#include <dumux/discretization/methods.hh>

#include "diffmethod.hh"
#include "staggeredlocalassembler.hh"
#include "staggeredlocalresidual.hh"

namespace Dumux {

/*!
 * \ingroup Assembly
 * \ingroup StaggeredDiscretization
 * \brief A linear system assembler (residual and Jacobian) for staggered finite volume schemes
 * \tparam TypeTag the TypeTag
 * \tparam diffMethod the differentiation method to residual compute derivatives
 * \tparam isImplicit if to use an implicit or explicit time discretization
 */
template<class TypeTag, DiffMethod diffMethod, bool isImplicit = true>
class StaggeredFVAssembler
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using LocalResidual = typename GET_PROP_TYPE(TypeTag, LocalResidual);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using GridVariables = typename GET_PROP_TYPE(TypeTag, GridVariables);
    using JacobianMatrix = typename GET_PROP_TYPE(TypeTag, JacobianMatrix);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using TimeLoop = TimeLoopBase<Scalar>;

    static constexpr int dim = GridView::dimension;

    using LocalAssembler =StaggeredLocalAssembler<TypeTag, diffMethod, isImplicit>;

    using DofTypeIndices = typename GET_PROP(TypeTag, DofTypeIndices);
    typename DofTypeIndices::CellCenterIdx cellCenterIdx;
    typename DofTypeIndices::FaceIdx faceIdx;

    using CCToCCMatrixBlock = typename GET_PROP(TypeTag, JacobianMatrix)::MatrixBlockCCToCC;
    using CCToFaceMatrixBlock = typename GET_PROP(TypeTag, JacobianMatrix)::MatrixBlockCCToFace;
    using FaceToFaceMatrixBlock = typename GET_PROP(TypeTag, JacobianMatrix)::MatrixBlockFaceToFace;
    using FaceToCCMatrixBlock = typename GET_PROP(TypeTag, JacobianMatrix)::MatrixBlockFaceToCC;

public:
    using ResidualType = SolutionVector;

    //! The constructor for stationary problems
    StaggeredFVAssembler(std::shared_ptr<const Problem> problem,
                std::shared_ptr<const FVGridGeometry> fvGridGeometry,
                std::shared_ptr<GridVariables> gridVariables)
    : problem_(problem)
    , fvGridGeometry_(fvGridGeometry)
    , gridVariables_(gridVariables)
    , stationary_(true)
    {
        static_assert(isImplicit, "Explicit assembler for stationary problem doesn't make sense!");
    }

    //! The constructor for instationary problems
    StaggeredFVAssembler(std::shared_ptr<const Problem> problem,
                std::shared_ptr<const FVGridGeometry> fvGridGeometry,
                std::shared_ptr<GridVariables> gridVariables,
                std::shared_ptr<TimeLoop> timeLoop)
    : problem_(problem)
    , fvGridGeometry_(fvGridGeometry)
    , gridVariables_(gridVariables)
    , localResidual_(timeLoop)
    , stationary_(false)
    {}

    /*!
     * \brief Assembles the global Jacobian of the residual
     *        and the residual for the current solution.
     */
    void assembleJacobianAndResidual(const SolutionVector& curSol)
    {
        if (!stationary_ && localResidual_.isStationary())
            DUNE_THROW(Dune::InvalidStateException, "Assembling instationary problem but previous solution was not set!");

        if(!jacobian_)
        {
            setLinearSystem();
        }

        if(!residual_)
        {
            residual_ = std::make_shared<SolutionVector>();
            setResidualSize();
        }

        resetJacobian_();
        resetResidual_();

        bool succeeded;
        // try assembling the global linear system
        try
        {
            // let the local assembler add the element contributions
            for (const auto& element : elements(gridView()))
                LocalAssembler::assembleJacobianAndResidual(*this, *jacobian_, *residual_, element, curSol);

            // if we get here, everything worked well
            succeeded = true;
            if (gridView().comm().size() > 1)
                succeeded = gridView().comm().min(succeeded);
        }
        // throw exception if a problem ocurred
        catch (NumericalProblem &e)
        {
            std::cout << "rank " << gridView().comm().rank()
                      << " caught an exception while assembling:" << e.what()
                      << "\n";
            succeeded = false;
            if (gridView().comm().size() > 1)
                succeeded = gridView().comm().min(succeeded);
        }
        if (!succeeded)
            DUNE_THROW(NumericalProblem, "A process did not succeed in linearizing the system");
    }


    //! compute the residuals
    void assembleResidual(const SolutionVector& curSol)
    {
        if(!residual_)
        {
            residual_ = std::make_shared<SolutionVector>();
            setResidualSize();
        }

        assembleResidual(*residual_, curSol);
    }

    //! compute the residuals
    void assembleResidual(ResidualType& r, const SolutionVector& curSol) const
    {
        if (!stationary_ && localResidual_.isStationary())
            DUNE_THROW(Dune::InvalidStateException, "Assembling instationary problem but previous solution was not set!");

        // let the local assembler add the element contributions
        for (const auto& element : elements(gridView()))
            LocalAssembler::assembleResidual(*this, r, element, curSol);
    }

    //! computes the residual norm
    Scalar residualNorm(const SolutionVector& curSol) const
    {
        ResidualType residual;
        residual[cellCenterIdx].resize(fvGridGeometry().numCellCenterDofs());
        residual[faceIdx].resize(fvGridGeometry().numFaceDofs());
        assembleResidual(residual, curSol);

        // calculate the square norm of the residual
        Scalar result2 = residual.two_norm2();
        if (gridView().comm().size() > 1)
            result2 = gridView().comm().sum(result2);

        using std::sqrt;
        return sqrt(result2);
    }

    /*!
     * \brief Tells the assembler which jacobian and residual to use.
     *        This also resizes the containers to the required sizes and sets the
     *        sparsity pattern of the jacobian matrix.
     */
    void setLinearSystem(std::shared_ptr<JacobianMatrix> A,
                         std::shared_ptr<SolutionVector> r)
    {
        jacobian_ = A;
        residual_ = r;

        // set the BCRS matrix's build mode
        // convenience references
        CCToCCMatrixBlock& A11 = (*jacobian_)[cellCenterIdx][cellCenterIdx];
        CCToFaceMatrixBlock& A12 = (*jacobian_)[cellCenterIdx][faceIdx];
        FaceToCCMatrixBlock& A21 = (*jacobian_)[faceIdx][cellCenterIdx];
        FaceToFaceMatrixBlock& A22 = (*jacobian_)[faceIdx][faceIdx];

        A11.setBuildMode(CCToCCMatrixBlock::random);
        A12.setBuildMode(CCToFaceMatrixBlock::random);
        A21.setBuildMode(FaceToCCMatrixBlock::random);
        A22.setBuildMode(FaceToFaceMatrixBlock::random);

        setJacobianPattern();
        setResidualSize();
    }

    /*!
     * \brief The version without arguments uses the default constructor to create
     *        the jacobian and residual objects in this assembler.
     */
    void setLinearSystem()
    {
        jacobian_ = std::make_shared<JacobianMatrix>();

        // convenience references
        CCToCCMatrixBlock& A11 = (*jacobian_)[cellCenterIdx][cellCenterIdx];
        CCToFaceMatrixBlock& A12 = (*jacobian_)[cellCenterIdx][faceIdx];
        FaceToCCMatrixBlock& A21 = (*jacobian_)[faceIdx][cellCenterIdx];
        FaceToFaceMatrixBlock& A22 = (*jacobian_)[faceIdx][faceIdx];

        A11.setBuildMode(CCToCCMatrixBlock::random);
        A12.setBuildMode(CCToFaceMatrixBlock::random);
        A21.setBuildMode(FaceToCCMatrixBlock::random);
        A22.setBuildMode(FaceToFaceMatrixBlock::random);

        residual_ = std::make_shared<SolutionVector>();

        setJacobianPattern();
        setResidualSize();
    }

    /*!
     * \brief Sets the solution from which to start the time integration. Has to be
     *        called prior to assembly for time-dependent problems.
     */
    void setPreviousSolution(const SolutionVector& u)
    { localResidual_.setPreviousSolution(u); }

    /*!
     * \brief Return the solution that has been set as the previous one.
     */
    const SolutionVector& prevSol() const
    { return localResidual_.prevSol(); }

    /*!
     * \brief Resizes the jacobian and sets the jacobian' sparsity pattern.
     */
    void setJacobianPattern()
    {
        // resize the jacobian and the residual
        const auto numDofsCC = fvGridGeometry().numCellCenterDofs();
        const auto numDofsFace = fvGridGeometry().numFaceDofs();

        // convenience references
        CCToCCMatrixBlock& A11 = (*jacobian_)[cellCenterIdx][cellCenterIdx];
        CCToFaceMatrixBlock& A12 = (*jacobian_)[cellCenterIdx][faceIdx];
        FaceToCCMatrixBlock& A21 = (*jacobian_)[faceIdx][cellCenterIdx];
        FaceToFaceMatrixBlock& A22 = (*jacobian_)[faceIdx][faceIdx];

        // set the size of the sub-matrizes
        A11.setSize(numDofsCC, numDofsCC);
        A12.setSize(numDofsCC, numDofsFace);
        A21.setSize(numDofsFace, numDofsCC);
        A22.setSize(numDofsFace, numDofsFace);


        // set occupation pattern of the jacobian
        Dune::MatrixIndexSet occupationPatternA11;
        Dune::MatrixIndexSet occupationPatternA12;
        Dune::MatrixIndexSet occupationPatternA21;
        Dune::MatrixIndexSet occupationPatternA22;
        occupationPatternA11.resize(numDofsCC, numDofsCC);
        occupationPatternA12.resize(numDofsCC, numDofsFace);
        occupationPatternA21.resize(numDofsFace, numDofsCC);
        occupationPatternA22.resize(numDofsFace, numDofsFace);

        const auto& connectivityMap = fvGridGeometry().connectivityMap();

        // evaluate the acutal pattern
        for (const auto& element : elements(fvGridGeometry().gridView()))
        {
            // the global index of the element at hand
            const auto ccGlobalI = fvGridGeometry().elementMapper().index(element);

            for (auto&& ccGlobalJ : connectivityMap(cellCenterIdx, cellCenterIdx, ccGlobalI))
                occupationPatternA11.add(ccGlobalI, ccGlobalJ);
            for (auto&& faceGlobalJ : connectivityMap(cellCenterIdx, faceIdx, ccGlobalI))
                occupationPatternA12.add(ccGlobalI, faceGlobalJ);

            auto fvGeometry = localView(fvGridGeometry());
            fvGeometry.bindElement(element);

            // loop over sub control faces
            for (auto&& scvf : scvfs(fvGeometry))
            {
                const auto faceGlobalI = scvf.dofIndex();
                for (auto&& ccGlobalJ : connectivityMap(faceIdx, cellCenterIdx, scvf.index()))
                    occupationPatternA21.add(faceGlobalI, ccGlobalJ);
                for (auto&& faceGlobalJ : connectivityMap(faceIdx, faceIdx, scvf.index()))
                    occupationPatternA22.add(faceGlobalI, faceGlobalJ);
            }
        }

        occupationPatternA11.exportIdx(A11);
        occupationPatternA12.exportIdx(A12);
        occupationPatternA21.exportIdx(A21);
        occupationPatternA22.exportIdx(A22);
    }

    /*!
     * \brief Resizes the residual
     */
    void setResidualSize()
    {
        (*residual_)[cellCenterIdx].resize(fvGridGeometry().numCellCenterDofs());
        (*residual_)[faceIdx].resize(fvGridGeometry().numFaceDofs());
    }

    const Problem& problem() const
    { return *problem_; }

    const FVGridGeometry& fvGridGeometry() const
    { return *fvGridGeometry_; }

    const GridView& gridView() const
    { return fvGridGeometry().gridView(); }

    GridVariables& gridVariables()
    { return *gridVariables_; }

    const GridVariables& gridVariables() const
    { return *gridVariables_; }

    JacobianMatrix& jacobian()
    {
       if (!residual_)
            DUNE_THROW(Dune::InvalidStateException, "No jacobian was set.");
        return *jacobian_;
    }

    SolutionVector& residual()
    {
        if (!residual_)
            DUNE_THROW(Dune::InvalidStateException, "No residual was set.");
        return *residual_;
    }

    const LocalResidual& localResidual() const
    { return localResidual_; }

    bool isStationaryProblem() const
    { return stationary_; }

private:
    //! reset the residual to 0.0
    void resetResidual_()
    {
        (*residual_)[cellCenterIdx] = 0.0;
        (*residual_)[faceIdx] = 0.0;
    }

    //! reset the jacobian to 0.0
    void resetJacobian_()
    {
        (*jacobian_)[cellCenterIdx][cellCenterIdx] = 0.0;
        (*jacobian_)[cellCenterIdx][faceIdx] = 0.0;
        (*jacobian_)[faceIdx][cellCenterIdx] = 0.0;
        (*jacobian_)[faceIdx][faceIdx] = 0.0;
    }

    //! pointer to the problem to be solved
    std::shared_ptr<const Problem> problem_;

    //! the finite volume geometry of the grid
    std::shared_ptr<const FVGridGeometry> fvGridGeometry_;

    //! the variables container for the grid
    std::shared_ptr<GridVariables> gridVariables_;

    //! shared pointers to the jacobian matrix and residual
    std::shared_ptr<JacobianMatrix> jacobian_;
    std::shared_ptr<SolutionVector> residual_;

    //! class computing the residual of an element
    LocalResidual localResidual_;

    //! if this assembler is assembling a time dependent problem
    bool stationary_;
};

} // namespace Dumux

#endif
