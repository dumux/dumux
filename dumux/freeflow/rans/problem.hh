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
#include <dumux/discretization/staggered/elementsolution.hh>
#include <dumux/discretization/methods.hh>
#include <dumux/freeflow/navierstokes/problem.hh>

#include "model.hh"

namespace Dumux
{

/*!
 * \ingroup RANSModel
 * \brief Reynolds-Averaged Navier-Stokes problem base class.
 *
 * This implements some base functionality for RANS models.
 * Especially vectors containing all wall-relevant properties, which are accessed
 * by the volumevariables.
 * \todo inherit all functions (especially gravity and temperature from Navier-Stokes)
 * This implements gravity (if desired) and a function returning the temperature.
 * Includes a specialized method used only by the staggered grid discretization.
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
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using CellCenterPrimaryVariables = typename GET_PROP_TYPE(TypeTag, CellCenterPrimaryVariables);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);

    enum {
        dim = Grid::dimension,
      };
    using GlobalPosition = Dune::FieldVector<Scalar, dim>;
    using DimVector = Dune::FieldVector<Scalar, dim>;
    using DimMatrix = Dune::FieldMatrix<Scalar, dim, dim>;

    enum {
        massBalanceIdx = Indices::massBalanceIdx,
        momentumBalanceIdx = Indices::momentumBalanceIdx
    };

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
     * \brief Update the static (solution independent) relations to the walls
     *
     * This function determines all element with a wall intersection,
     * the wall distances and the relation to the neighboring elements.
     */
    void updateStaticWallProperties() const
    {
        using std::abs;
        std::cout << "Update static wall properties. ";

        // update size and initial values of the global vectors
        wallElementIDs_.resize(this->fvGridGeometry().elementMapper().size());
        wallDistances_.resize(this->fvGridGeometry().elementMapper().size());
        neighborIDs_.resize(this->fvGridGeometry().elementMapper().size());
        cellCenters_.resize(this->fvGridGeometry().elementMapper().size());
        velocity_.resize(this->fvGridGeometry().elementMapper().size());
        velocityGradients_.resize(this->fvGridGeometry().elementMapper().size());
        kinematicViscosity_.resize(this->fvGridGeometry().elementMapper().size());
        for (unsigned int i = 0; i < wallElementIDs_.size(); ++i)
        {
            wallDistances_[i] = std::numeric_limits<Scalar>::max();
            cellCenters_[i] = GlobalPosition(0.0);
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
            unsigned int elementID = this->fvGridGeometry().elementMapper().index(element);
            cellCenters_[elementID] = element.geometry().center();
            for (unsigned int i = 0; i < wallPositions.size(); ++i)
            {
                GlobalPosition global = element.geometry().center();
                global -= wallPositions[i];
                if (global.two_norm() < wallDistances_[elementID])
                {
                    wallDistances_[elementID] = global.two_norm();
                    wallElementIDs_[elementID] = wallElements[i];
                }
            }
        }

        // search for neighbor IDs
        for (const auto& element : elements(gridView))
        {
            unsigned int elementID = this->fvGridGeometry().elementMapper().index(element);
            for (unsigned int dimIdx = 0; dimIdx < dim; ++dimIdx)
            {
                neighborIDs_[elementID][dimIdx][0] = elementID;
                neighborIDs_[elementID][dimIdx][1] = elementID;
            }

            for (const auto& intersection : intersections(gridView, element))
            {
                if (intersection.boundary())
                    continue;

                unsigned int neighborID = this->fvGridGeometry().elementMapper().index(intersection.outside());
                for (unsigned int dimIdx = 0; dimIdx < dim; ++dimIdx)
                {
                    if (abs(cellCenters_[elementID][dimIdx] - cellCenters_[neighborID][dimIdx]) > 1e-8)
                    {
                        if (cellCenters_[elementID][dimIdx] > cellCenters_[neighborID][dimIdx])
                        {
                            neighborIDs_[elementID][dimIdx][0] = neighborID;
                        }
                        if (cellCenters_[elementID][dimIdx] < cellCenters_[neighborID][dimIdx])
                        {
                            neighborIDs_[elementID][dimIdx][1] = neighborID;
                        }
                    }
                }
            }
        }
    }



    /*!
     * \brief Update the dynamic (solution dependent) relations to the walls
     *
     * The basic function calcuates the cell-centered velocities and
     * the respective gradients.
     * Further, the kinematic viscosity at the wall is stored.
     *
     * \param curSol The solution vector.
     */
    void updateDynamicWallProperties(const SolutionVector& curSol) const
    {
        std::cout << "Update dynamic wall properties." << std::endl;

        // calculate cell-center-averaged velocities
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
                const auto numericalSolutionFace = curSol[faceIdx][dofIdxFace][momentumBalanceIdx];
                velocityTemp[scvf.directionIndex()] += numericalSolutionFace;
            }
            for (unsigned int dimIdx = 0; dimIdx < dim; ++dimIdx)
                velocity_[elementID][dimIdx] = velocityTemp[dimIdx] * 0.5; // faces are equidistant to cell center
        }

        // calculate cell-center-averaged velocity gradients
        for (const auto& element : elements(this->fvGridGeometry().gridView()))
        {
            unsigned int elementID = this->fvGridGeometry().elementMapper().index(element);
            for (unsigned int dimIdx = 0; dimIdx < dim; ++dimIdx)
            {
                for (unsigned int velIdx = 0; velIdx < dim; ++velIdx)
                {
                    velocityGradients_[elementID][velIdx][dimIdx]
                        = (velocity_[neighborIDs_[elementID][dimIdx][1]][velIdx]
                              - velocity_[neighborIDs_[elementID][dimIdx][0]][velIdx])
                          / (cellCenters_[neighborIDs_[elementID][dimIdx][1]][dimIdx]
                              - cellCenters_[neighborIDs_[elementID][dimIdx][0]][dimIdx]);
                }
            }
        }

        // calculate or call all secondary variables
        for (const auto& element : elements(this->fvGridGeometry().gridView()))
        {
            unsigned int elementID = this->fvGridGeometry().elementMapper().index(element);

            auto fvGeometry = localView(this->fvGridGeometry());
            fvGeometry.bindElement(element);
            for (auto&& scv : scvs(fvGeometry))
            {
                const int dofIdx = scv.dofIndex();
                CellCenterPrimaryVariables priVars(curSol[cellCenterIdx][dofIdx]);
                auto elemSol = elementSolution<SolutionVector, FVGridGeometry>(std::move(priVars));
                VolumeVariables volVars;
                volVars.update(elemSol, asImp_(), element,scv);
                kinematicViscosity_[elementID] = volVars.viscosity() / volVars.density();
            }
        }
    }

    /*!
     * \brief Returns whether a given point is on a wall
     *
     * \param globalPos The position in global coordinates where the temperature should be specified.
     */
    bool isOnWall(const GlobalPosition &globalPos) const
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
    mutable std::vector<std::array<std::array<unsigned int, 2>, dim>> neighborIDs_;
    mutable std::vector<GlobalPosition> cellCenters_;
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
