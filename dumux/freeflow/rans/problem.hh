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
 * \ingroup RANSModel
 * \copydoc Dumux::RANSProblem
 */
#ifndef DUMUX_RANS_PROBLEM_HH
#define DUMUX_RANS_PROBLEM_HH

#include <dumux/common/properties.hh>
#include <dumux/common/staggeredfvproblem.hh>
#include <dumux/discretization/localview.hh>
#include <dumux/discretization/methods.hh>
#include <dumux/freeflow/navierstokes/problem.hh>

#include "model.hh"

namespace Dumux
{

/*!
 * \ingroup RANSModel
 * \brief Reynolds-Averaged Navier-Stokes problem base class.
 *
 * This implements gravity (if desired) and a function returning the temperature.
 * Includes a specialized method used only by the staggered grid discretization.
 *
 * \todo inherit all functions (especially gravity and temperature from Navier-Stokes)
 */
template<class TypeTag>
class RANSProblem : public NavierStokesParentProblem<TypeTag>
{
    using ParentType = NavierStokesParentProblem<TypeTag>;
    using Implementation = typename GET_PROP_TYPE(TypeTag, Problem);

    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Grid = typename GridView::Grid;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);

    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);

    enum {
        dim = Grid::dimension,
        dimWorld = Grid::dimensionworld
      };
    // TODO: dim or dimWorld appropriate here?
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
    using DimVector = Dune::FieldVector<Scalar, dimWorld>;
    using DimMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;

    enum {
        massBalanceIdx = Indices::massBalanceIdx,
        momentumBalanceIdx = Indices::momentumBalanceIdx/*,
        momentumXBalanceIdx = Indices::momentumXBalanceIdx,
        momentumYBalanceIdx = Indices::momentumYBalanceIdx,
        pressureIdx = Indices::pressureIdx,
        velocityXIdx = Indices::velocityXIdx,
        velocityYIdx = Indices::velocityYIdx*/
    };

//     using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
//
//     using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
//
//     using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
//
//     using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
//     using NumEqVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);
//
//     using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
//
    using DofTypeIndices = typename GET_PROP(TypeTag, DofTypeIndices);
    typename DofTypeIndices::CellCenterIdx cellCenterIdx;
    typename DofTypeIndices::FaceIdx faceIdx;

public:
    //! The constructor sets the gravity, if desired by the user.
    RANSProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
    {
        if (getParamFromGroup<bool>(GET_PROP_VALUE(TypeTag, ModelParameterGroup), "Problem.EnableGravity"))
            gravity_[dim-1]  = -9.81;

        asImp_().updateStaticWallProperties();
    }

    /*!
     * \brief \todo please doc me
     */
    void updateStaticWallProperties() const
    {
        std::cout << "Update static wall properties. ";

        // update size and initial values of the global vectors
        wallElementIDs_.resize(this->fvGridGeometry().elementMapper().size());
        wallDistances_.resize(this->fvGridGeometry().elementMapper().size());
        velocity_.resize(this->fvGridGeometry().elementMapper().size());
        velocityGradients_.resize(this->fvGridGeometry().elementMapper().size());
        kinematicViscosity_.resize(this->fvGridGeometry().elementMapper().size());
        for (unsigned int i = 0; i < wallElementIDs_.size(); ++i)
        {
            wallDistances_[i] = std::numeric_limits<Scalar>::max();
            velocity_[i] = DimVector(0.0);
            velocityGradients_[i] = DimMatrix(0.0);
            kinematicViscosity_[i] = 0.0;
        }

        // retrieve all wall intersections and corresponding elements
        std::vector<unsigned int> wallElements;
        std::vector<GlobalPosition> wallPositions;
        auto& gridView(this->fvGridGeometry().gridView());
        for (const auto& element : elements(gridView))
        {
            for (const auto& intersection : intersections(gridView, element))
            {
                GlobalPosition global = intersection.geometry().center();
                if (asImp_().isOnWall(global))
                {
                    wallElements.push_back(this->fvGridGeometry().elementMapper().index(element));
                    wallPositions.push_back(global);
                }
            }
        }
        std::cout << "NumWallIntersections=" << wallPositions.size() << std::endl;

        // search for shortest distance to wall for each element
        for (const auto& element : elements(gridView))
        {
            for (unsigned int i = 0; i < wallPositions.size(); ++i)
            {
                GlobalPosition global = element.geometry().center();
                global -= wallPositions[i];
                unsigned int elementID = this->fvGridGeometry().elementMapper().index(element);
                if (global.two_norm() < wallDistances_[elementID])
                {
                    wallDistances_[elementID] = global.two_norm();
                    wallElementIDs_[elementID] = wallElements[i];
                }
            }
        }
    }

    /*!
     * \brief \todo please doc me
     */
    void updateDynamicWallProperties(const SolutionVector& curSol) const
    {
        std::cout << "Update dynamic wall properties." << std::endl;

        for (const auto& element : elements(this->fvGridGeometry().gridView()))
        {
            auto fvGeometry = localView(this->fvGridGeometry());
            fvGeometry.bindElement(element);
            unsigned int elementID = this->fvGridGeometry().elementMapper().index(element);

            // calculate velocities
            DimVector velocityTemp(0.0);
            for (auto&& scvf : scvfs(fvGeometry))
            {
                const int dofIdxFace = scvf.dofIndex();
                const int dirIdx = scvf.directionIndex();
                const auto numericalSolutionFace = curSol[faceIdx][dofIdxFace][momentumBalanceIdx];
                velocityTemp[scvf.directionIndex()] += numericalSolutionFace;
            }
            for (unsigned int dimIdx = 0; dimIdx < dim; ++dimIdx)
                velocity_[elementID][dimIdx] = velocityTemp[dimIdx] * 0.5; // faces are equidistant to cell center

//             // TODO: calculate velocity gradients
//             for (auto&& scv : scvs(fvGeometry))
//             {
//                 // treat cell-center dofs
//                 const auto dofIdxCellCenter = scv.dofIndex();
//                 const auto& posCellCenter = scv.dofPosition();
//                 const auto analyticalSolutionCellCenter = analyticalSolution(posCellCenter)[pressureIdx];
//                 const auto numericalSolutionCellCenter = curSol[cellCenterIdx][dofIdxCellCenter][pressureIdx];
//                 sumError[pressureIdx] += squaredDiff_(analyticalSolutionCellCenter, numericalSolutionCellCenter) * scv.volume();
//                 sumReference[pressureIdx] += analyticalSolutionCellCenter * analyticalSolutionCellCenter * scv.volume();
//                 totalVolume += scv.volume();
//
//
//                 // treat face dofs
//                 for (auto&& scvf : scvfs(fvGeometry))
//                 {
//                     const int dofIdxFace = scvf.dofIndex();
//                     const int dirIdx = scvf.directionIndex();
//                     const auto analyticalSolutionFace = analyticalSolution(scvf.center())[Indices::velocity(dirIdx)];
//                     const auto numericalSolutionFace = curSol[faceIdx][dofIdxFace][momentumBalanceIdx];
//                     directionIndex[dofIdxFace] = dirIdx;
//                     errorVelocity[dofIdxFace] = squaredDiff_(analyticalSolutionFace, numericalSolutionFace);
//                     velocityReference[dofIdxFace] = squaredDiff_(analyticalSolutionFace, 0.0);
//                     const Scalar staggeredHalfVolume = 0.5 * scv.volume();
//                     staggeredVolume[dofIdxFace] = staggeredVolume[dofIdxFace] + staggeredHalfVolume;
//                     const int dofIdxFace = scvf.dofIndex();
//                     const int dirIdx = scvf.directionIndex();
//                     const auto numericalSolutionFace = curSol[faceIdx][dofIdxFace][momentumBalanceIdx];
//                     velocityTemp[scvf.directionIndex()] += numericalSolutionFace;
//                 }
//             }
            for (unsigned int dimIdx = 0; dimIdx < dim; ++dimIdx)
                for (unsigned int velIdx = 0; velIdx < dim; ++velIdx)
                    velocityGradients_[elementID][velIdx][dimIdx] = 0.12345;

//             // TODO calculate or call all secondary variables
//             // TODO call kinematic viscosity value from vol vars
            kinematicViscosity_[elementID] = 15e-6;
        }
    }

    /*!
     * \brief Returns whether a given point is on a wall
     */
    const bool isOnWall(GlobalPosition &globalPos) const
    {
        // Throw an exception if no walls are implemented
        DUNE_THROW(Dune::InvalidStateException,
                   "The problem does not provide an isOnWall() method.");
    }

    /*!
     * \brief Returns the temperature \f$\mathrm{[K]}\f$ at a given global position.
     *
     * This is not specific to the discretization. By default it just
     * calls temperature().
     *
     * \param globalPos The position in global coordinates where the temperature should be specified.
     */
    Scalar temperatureAtPos(const GlobalPosition &globalPos) const
    { return asImp_().temperature(); }

    /*!
     * \brief Returns the temperature within the domain.
     *
     * This method MUST be overwritten by the actual problem.
     */
    Scalar temperature() const
    { DUNE_THROW(Dune::NotImplemented, "temperature() method not implemented by the actual problem"); }

    /*!
     * \brief Returns the acceleration due to gravity.
     *
     * If the <tt>Problem.EnableGravity</tt> parameter is true, this means
     * \f$\boldsymbol{g} = ( 0,\dots,\ -9.81)^T \f$, else \f$\boldsymbol{g} = ( 0,\dots, 0)^T \f$
     */
    const GlobalPosition &gravity() const
    { return gravity_; }

    //! Applys the initial face solution (velocities on the faces). Specialization for staggered grid discretization.
    template <class T = TypeTag>
    typename std::enable_if<GET_PROP_TYPE(T, FVGridGeometry)::discMethod == DiscretizationMethod::staggered, void>::type
    applyInititalFaceSolution(SolutionVector& sol,
                              const SubControlVolumeFace& scvf,
                              const PrimaryVariables& initSol) const
    {
        typename GET_PROP(TypeTag, DofTypeIndices)::FaceIdx faceIdx;
        const auto numEqCellCenter = GET_PROP_VALUE(TypeTag, NumEqCellCenter);
        sol[faceIdx][scvf.dofIndex()][numEqCellCenter] = initSol[Indices::velocity(scvf.directionIndex())];
    }

public:
    mutable std::vector<unsigned int> wallElementIDs_;
    mutable std::vector<Scalar> wallDistances_;
    mutable std::vector<DimVector> velocity_;
    mutable std::vector<DimMatrix> velocityGradients_;
    mutable std::vector<Scalar> kinematicViscosity_;

private:
    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }

    GlobalPosition gravity_;
};

}

#endif
