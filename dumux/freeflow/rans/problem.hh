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

namespace Dumux {

/*!
 * \ingroup RANSModel
 * \brief Reynolds-Averaged Navier-Stokes problem base class.
 *
 * This implements some base functionality for RANS models.
 * Especially vectors containing all wall-relevant properties, which are accessed
 * by the volumevariables.
 */
template<class TypeTag>
class RANSProblem : public NavierStokesProblem<TypeTag>
{
    using ParentType = NavierStokesProblem<TypeTag>;
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
    using Indices = typename GET_PROP_TYPE(TypeTag, ModelTraits)::Indices;

    enum {
        dim = Grid::dimension,
      };
    using GlobalPosition = typename SubControlVolumeFace::GlobalPosition;
    using DimVector = Dune::FieldVector<Scalar, dim>;
    using DimMatrix = Dune::FieldMatrix<Scalar, dim, dim>;

public:
    /*!
     * \brief The constructor
     * \param fvGridGeometry The finite volume grid geometry
     * \param paramGroup The parameter group in which to look for runtime parameters first (default is "")
     */
    RANSProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry, const std::string& paramGroup = "")
    : ParentType(fvGridGeometry, paramGroup)
    { }

    /*!
     * \brief Update the static (solution independent) relations to the walls
     *
     * This function determines all element with a wall intersection,
     * the wall distances and the relation to the neighboring elements.
     */
    void updateStaticWallProperties()
    {
        using std::abs;
        std::cout << "Update static wall properties. ";

        // update size and initial values of the global vectors
        wallElementID_.resize(this->fvGridGeometry().elementMapper().size());
        wallDistance_.resize(this->fvGridGeometry().elementMapper().size(), std::numeric_limits<Scalar>::max());
        neighborID_.resize(this->fvGridGeometry().elementMapper().size());
        cellCenter_.resize(this->fvGridGeometry().elementMapper().size(), GlobalPosition(0.0));
        velocity_.resize(this->fvGridGeometry().elementMapper().size(), DimVector(0.0));
        velocityMaximum_.resize(this->fvGridGeometry().elementMapper().size(), DimVector(0.0));
        velocityGradients_.resize(this->fvGridGeometry().elementMapper().size(), DimMatrix(0.0));
        stressTensorScalarProduct_.resize(this->fvGridGeometry().elementMapper().size(), 0.0);
        flowNormalAxis_.resize(this->fvGridGeometry().elementMapper().size(), 0);
        wallNormalAxis_.resize(this->fvGridGeometry().elementMapper().size(), 1);
        kinematicViscosity_.resize(this->fvGridGeometry().elementMapper().size(), 0.0);
        sandGrainRoughness_.resize(this->fvGridGeometry().elementMapper().size(), 0.0);

        // retrieve all wall intersections and corresponding elements
        std::vector<unsigned int> wallElements;
        std::vector<GlobalPosition> wallPositions;
        std::vector<unsigned int> wallNormalAxisTemp;
        auto& gridView(this->fvGridGeometry().gridView());
        for (const auto& element : elements(gridView))
        {
            for (const auto& intersection : intersections(gridView, element))
            {
                // only search for walls at a global boundary
                if (!intersection.boundary())
                    continue;

                GlobalPosition global = intersection.geometry().center();
                if (asImp_().isOnWall(global))
                {
                    wallElements.push_back(this->fvGridGeometry().elementMapper().index(element));
                    wallPositions.push_back(global);
                    for (unsigned int dimIdx = 0; dimIdx < dim; ++dimIdx)
                    {
                        if (abs(intersection.centerUnitOuterNormal()[dimIdx]) > 1e-8)
                            wallNormalAxisTemp.push_back(dimIdx);
                    }
                }
            }
        }
        std::cout << "NumWallIntersections=" << wallPositions.size() << std::endl;

        // search for shortest distance to wall for each element
        for (const auto& element : elements(gridView))
        {
            unsigned int elementID = this->fvGridGeometry().elementMapper().index(element);
            cellCenter_[elementID] = element.geometry().center();
            for (unsigned int i = 0; i < wallPositions.size(); ++i)
            {
                static const int problemWallNormalAxis
                    = getParamFromGroup<int>(this->paramGroup(), "RANS.WallNormalAxis", -1);
                int searchAxis = problemWallNormalAxis;

                // search along wall normal axis of the intersection
                if (problemWallNormalAxis < 0 || problemWallNormalAxis >= dim)
                {
                    searchAxis = wallNormalAxisTemp[i];
                }

                GlobalPosition global = element.geometry().center();
                global -= wallPositions[i];
                // second and argument ensures to use only aligned elements
                if (abs(global[searchAxis]) < wallDistance_[elementID]
                    && abs(global[searchAxis]) < global.two_norm() + 1e-8
                    && abs(global[searchAxis]) > global.two_norm() - 1e-8)
                {
                    wallDistance_[elementID] = abs(global[searchAxis]);
                    wallElementID_[elementID] = wallElements[i];
                    wallNormalAxis_[elementID] = searchAxis;
                    sandGrainRoughness_[elementID] = asImp_().sandGrainRoughnessAtPos(wallPositions[i]);
                }
            }
        }

        // search for neighbor IDs
        for (const auto& element : elements(gridView))
        {
            unsigned int elementID = this->fvGridGeometry().elementMapper().index(element);
            for (unsigned int dimIdx = 0; dimIdx < dim; ++dimIdx)
            {
                neighborID_[elementID][dimIdx][0] = elementID;
                neighborID_[elementID][dimIdx][1] = elementID;
            }

            for (const auto& intersection : intersections(gridView, element))
            {
                if (intersection.boundary())
                    continue;

                unsigned int neighborID = this->fvGridGeometry().elementMapper().index(intersection.outside());
                for (unsigned int dimIdx = 0; dimIdx < dim; ++dimIdx)
                {
                    if (abs(cellCenter_[elementID][dimIdx] - cellCenter_[neighborID][dimIdx]) > 1e-8)
                    {
                        if (cellCenter_[elementID][dimIdx] > cellCenter_[neighborID][dimIdx])
                        {
                            neighborID_[elementID][dimIdx][0] = neighborID;
                        }
                        if (cellCenter_[elementID][dimIdx] < cellCenter_[neighborID][dimIdx])
                        {
                            neighborID_[elementID][dimIdx][1] = neighborID;
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
    void updateDynamicWallProperties(const SolutionVector& curSol)
    {
        using std::abs;
        using std::max;
        using std::min;
        std::cout << "Update dynamic wall properties." << std::endl;

        static const int flowNormalAxis
            = getParamFromGroup<int>(this->paramGroup(), "RANS.FlowNormalAxis", -1);

        // re-initialize min and max values
        velocityMaximum_.assign(this->fvGridGeometry().elementMapper().size(), DimVector(1e-16));
        velocityMinimum_.assign(this->fvGridGeometry().elementMapper().size(), DimVector(std::numeric_limits<Scalar>::max()));

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
                const auto numericalSolutionFace = curSol[FVGridGeometry::faceIdx()][dofIdxFace][Indices::velocity(scvf.directionIndex())];
                velocityTemp[scvf.directionIndex()] += numericalSolutionFace;
            }
            for (unsigned int dimIdx = 0; dimIdx < dim; ++dimIdx)
                velocity_[elementID][dimIdx] = velocityTemp[dimIdx] * 0.5; // faces are equidistant to cell center
        }

        // calculate cell-center-averaged velocity gradients, maximum, and minimum values
        for (const auto& element : elements(this->fvGridGeometry().gridView()))
        {
            unsigned int elementID = this->fvGridGeometry().elementMapper().index(element);
            unsigned int wallElementID = wallElementID_[elementID];

            Scalar maxVelocity = 0.0;
            for (unsigned int dimIdx = 0; dimIdx < dim; ++dimIdx)
            {
                for (unsigned int velIdx = 0; velIdx < dim; ++velIdx)
                {
                    velocityGradients_[elementID][velIdx][dimIdx]
                        = (velocity_[neighborID_[elementID][dimIdx][1]][velIdx]
                              - velocity_[neighborID_[elementID][dimIdx][0]][velIdx])
                          / (cellCenter_[neighborID_[elementID][dimIdx][1]][dimIdx]
                              - cellCenter_[neighborID_[elementID][dimIdx][0]][dimIdx]);
                }

                if (abs(velocity_[elementID][dimIdx]) > abs(velocityMaximum_[wallElementID][dimIdx]))
                {
                    velocityMaximum_[wallElementID][dimIdx] = velocity_[elementID][dimIdx];
                }
                if (abs(velocity_[elementID][dimIdx]) < abs(velocityMinimum_[wallElementID][dimIdx]))
                {
                    velocityMinimum_[wallElementID][dimIdx] = velocity_[elementID][dimIdx];
                }

                if (0 <= flowNormalAxis && flowNormalAxis < dim)
                {
                    flowNormalAxis_[elementID] = flowNormalAxis;
                }
                else if (abs(maxVelocity) < abs(velocity_[elementID][dimIdx]))
                {
                    maxVelocity = abs(velocity_[elementID][dimIdx]);
                    flowNormalAxis_[elementID] = dimIdx;
                }
            }

            auto fvGeometry = localView(this->fvGridGeometry());
            fvGeometry.bindElement(element);
            for (auto&& scvf : scvfs(fvGeometry))
            {
                unsigned int normDim = scvf.directionIndex();
                if (scvf.boundary() && asImp_().boundaryTypes(element, scvf).isDirichlet(Indices::velocity(normDim)))
                {
                    for (unsigned int velIdx = 0; velIdx < dim; ++velIdx)
                    {
                        // face Value
                        Scalar dirichletVelocity = asImp_().dirichlet(element, scvf)[Indices::velocity(velIdx)];

                        unsigned int neighborID = neighborID_[elementID][normDim][0];
                        if (scvf.center()[normDim] < cellCenter_[elementID][normDim])
                            neighborID = neighborID_[elementID][normDim][1];

                        velocityGradients_[elementID][velIdx][normDim]
                            = (velocity_[neighborID][velIdx] - dirichletVelocity)
                              / (cellCenter_[neighborID][normDim] - scvf.center()[normDim]);
                    }
                }
            }
        }

        // calculate or call all secondary variables
        for (const auto& element : elements(this->fvGridGeometry().gridView()))
        {
            unsigned int elementID = this->fvGridGeometry().elementMapper().index(element);

            Dune::FieldMatrix<Scalar, GridView::dimension, GridView::dimension> stressTensor(0.0);
            for (unsigned int dimIdx = 0; dimIdx < dim; ++dimIdx)
            {
                for (unsigned int velIdx = 0; velIdx < dim; ++velIdx)
                {
                    stressTensor[dimIdx][velIdx] = 0.5 * velocityGradients_[elementID][dimIdx][velIdx]
                                                   + 0.5 * velocityGradients_[elementID][velIdx][dimIdx];
              }
            }
            stressTensorScalarProduct_[elementID] = 0.0;
            for (unsigned int dimIdx = 0; dimIdx < dim; ++dimIdx)
            {
                for (unsigned int velIdx = 0; velIdx < dim; ++velIdx)
                {
                    stressTensorScalarProduct_[elementID] += stressTensor[dimIdx][velIdx] * stressTensor[dimIdx][velIdx];
                }
            }

            auto fvGeometry = localView(this->fvGridGeometry());
            fvGeometry.bindElement(element);
            for (auto&& scv : scvs(fvGeometry))
            {
                const int dofIdx = scv.dofIndex();

                // construct a privars object from the cell center solution vector
                const auto& cellCenterPriVars = curSol[FVGridGeometry::cellCenterIdx()][dofIdx];
                PrimaryVariables priVars = makePriVarsFromCellCenterPriVars<PrimaryVariables>(cellCenterPriVars);
                auto elemSol = elementSolution<typename FVGridGeometry::LocalView>(std::move(priVars));

                VolumeVariables volVars;
                volVars.update(elemSol, asImp_(), element, scv);
                kinematicViscosity_[elementID] = volVars.viscosity() / volVars.density();
            }
        }
    }

    /*!
     * \brief Returns whether a given point is on a wall
     *
     * \param globalPos The position in global coordinates.
     */
    bool isOnWall(const GlobalPosition &globalPos) const
    {
        // Throw an exception if no walls are implemented
        DUNE_THROW(Dune::InvalidStateException,
                   "The problem does not provide an isOnWall() method.");
    }

    /*!
     * \brief Returns the sand-grain roughness \f$\mathrm{[m]}\f$ at a given position
     *
     * \param globalPos The position in global coordinates.
     */
    Scalar sandGrainRoughnessAtPos(const GlobalPosition &globalPos) const
    {
        return 0.0;
    }

    /*!
     * \brief Returns the Karman constant
     */
    const Scalar karmanConstant() const
    { return 0.41; }

public:
    std::vector<unsigned int> wallElementID_;
    std::vector<Scalar> wallDistance_;
    std::vector<std::array<std::array<unsigned int, 2>, dim>> neighborID_;
    std::vector<GlobalPosition> cellCenter_;
    std::vector<DimVector> velocity_;
    std::vector<DimVector> velocityMaximum_;
    std::vector<DimVector> velocityMinimum_;
    std::vector<DimMatrix> velocityGradients_;
    std::vector<Scalar> stressTensorScalarProduct_;
    std::vector<unsigned int> wallNormalAxis_;
    std::vector<unsigned int> flowNormalAxis_;
    std::vector<Scalar> kinematicViscosity_;
    std::vector<Scalar> sandGrainRoughness_;

private:
    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }
};

} // end namespace Dumux

#endif
