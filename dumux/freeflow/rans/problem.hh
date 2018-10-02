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

#include <dune/common/fmatrix.hh>
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

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);

    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using GridView = typename FVGridGeometry::GridView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using PrimaryVariables = typename VolumeVariables::PrimaryVariables;
    using CellCenterPrimaryVariables = typename GET_PROP_TYPE(TypeTag, CellCenterPrimaryVariables);
    using Indices = typename GET_PROP_TYPE(TypeTag, ModelTraits)::Indices;

    using GlobalPosition = typename SubControlVolumeFace::GlobalPosition;

    static constexpr auto dim = GridView::dimension;
    using DimVector = GlobalPosition;
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
        calledUpdateStaticWallProperties = true;

        // update size and initial values of the global vectors
        wallElementIdx_.resize(this->fvGridGeometry().elementMapper().size());
        wallDistance_.resize(this->fvGridGeometry().elementMapper().size(), std::numeric_limits<Scalar>::max());
        neighborIdx_.resize(this->fvGridGeometry().elementMapper().size());
        cellCenter_.resize(this->fvGridGeometry().elementMapper().size(), GlobalPosition(0.0));
        velocity_.resize(this->fvGridGeometry().elementMapper().size(), DimVector(0.0));
        velocityMaximum_.resize(this->fvGridGeometry().elementMapper().size(), DimVector(0.0));
        velocityGradients_.resize(this->fvGridGeometry().elementMapper().size(), DimMatrix(0.0));
        stressTensorScalarProduct_.resize(this->fvGridGeometry().elementMapper().size(), 0.0);
        vorticityTensorScalarProduct_.resize(this->fvGridGeometry().elementMapper().size(), 0.0);
        flowNormalAxis_.resize(this->fvGridGeometry().elementMapper().size(), 0);
        wallNormalAxis_.resize(this->fvGridGeometry().elementMapper().size(), 1);
        kinematicViscosity_.resize(this->fvGridGeometry().elementMapper().size(), 0.0);
        sandGrainRoughness_.resize(this->fvGridGeometry().elementMapper().size(), 0.0);

        // retrieve all wall intersections and corresponding elements
        std::vector<unsigned int> wallElements;
        std::vector<GlobalPosition> wallPositions;
        std::vector<unsigned int> wallNormalAxisTemp;

        const auto gridView = this->fvGridGeometry().gridView();
        auto fvGeometry = localView(this->fvGridGeometry());

        for (const auto& element : elements(gridView))
        {
            fvGeometry.bindElement(element);
            for (const auto& scvf : scvfs(fvGeometry))
            {
                // only search for walls at a global boundary
                if (!scvf.boundary())
                    continue;

                if (asImp_().isOnWall(scvf))
                {
                    wallElements.push_back(this->fvGridGeometry().elementMapper().index(element));
                    wallPositions.push_back(scvf.center());
                    wallNormalAxisTemp.push_back(scvf.directionIndex());
                }
            }
        }
        std::cout << "NumWallIntersections=" << wallPositions.size() << std::endl;
        if (wallPositions.size() == 0)
            DUNE_THROW(Dune::InvalidStateException,
                       "No wall intersections have been found. Make sure that the isOnWall(globalPos) is working properly.");


        // search for shortest distance to wall for each element
        for (const auto& element : elements(gridView))
        {
            unsigned int elementIdx = this->fvGridGeometry().elementMapper().index(element);
            cellCenter_[elementIdx] = element.geometry().center();
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
                if (abs(global[searchAxis]) < wallDistance_[elementIdx]
                    && abs(global[searchAxis]) < global.two_norm() + 1e-8
                    && abs(global[searchAxis]) > global.two_norm() - 1e-8)
                {
                    wallDistance_[elementIdx] = abs(global[searchAxis]);
                    wallElementIdx_[elementIdx] = wallElements[i];
                    wallNormalAxis_[elementIdx] = searchAxis;
                    sandGrainRoughness_[elementIdx] = asImp_().sandGrainRoughnessAtPos(wallPositions[i]);
                }
            }
        }

        // search for neighbor Idxs
        for (const auto& element : elements(gridView))
        {
            unsigned int elementIdx = this->fvGridGeometry().elementMapper().index(element);
            for (unsigned int dimIdx = 0; dimIdx < dim; ++dimIdx)
            {
                neighborIdx_[elementIdx][dimIdx][0] = elementIdx;
                neighborIdx_[elementIdx][dimIdx][1] = elementIdx;
            }

            for (const auto& intersection : intersections(gridView, element))
            {
                if (intersection.boundary())
                    continue;

                unsigned int neighborIdx = this->fvGridGeometry().elementMapper().index(intersection.outside());
                for (unsigned int dimIdx = 0; dimIdx < dim; ++dimIdx)
                {
                    if (abs(cellCenter_[elementIdx][dimIdx] - cellCenter_[neighborIdx][dimIdx]) > 1e-8)
                    {
                        if (cellCenter_[elementIdx][dimIdx] > cellCenter_[neighborIdx][dimIdx])
                        {
                            neighborIdx_[elementIdx][dimIdx][0] = neighborIdx;
                        }
                        if (cellCenter_[elementIdx][dimIdx] < cellCenter_[neighborIdx][dimIdx])
                        {
                            neighborIdx_[elementIdx][dimIdx][1] = neighborIdx;
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
        if (!calledUpdateStaticWallProperties)
            DUNE_THROW(Dune::InvalidStateException,
                       "You have to call updateStaticWallProperties() once before you call updateDynamicWallProperties().");

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
            unsigned int elementIdx = this->fvGridGeometry().elementMapper().index(element);

            // calculate velocities
            DimVector velocityTemp(0.0);
            for (auto&& scvf : scvfs(fvGeometry))
            {
                const int dofIdxFace = scvf.dofIndex();
                const auto numericalSolutionFace = curSol[FVGridGeometry::faceIdx()][dofIdxFace][Indices::velocity(scvf.directionIndex())];
                velocityTemp[scvf.directionIndex()] += numericalSolutionFace;
            }
            for (unsigned int dimIdx = 0; dimIdx < dim; ++dimIdx)
                velocity_[elementIdx][dimIdx] = velocityTemp[dimIdx] * 0.5; // faces are equidistant to cell center
        }

        // calculate cell-center-averaged velocity gradients, maximum, and minimum values
        for (const auto& element : elements(this->fvGridGeometry().gridView()))
        {
            unsigned int elementIdx = this->fvGridGeometry().elementMapper().index(element);
            unsigned int wallElementIdx = wallElementIdx_[elementIdx];

            Scalar maxVelocity = 0.0;
            for (unsigned int dimIdx = 0; dimIdx < dim; ++dimIdx)
            {
                for (unsigned int velIdx = 0; velIdx < dim; ++velIdx)
                {
                    velocityGradients_[elementIdx][velIdx][dimIdx]
                        = (velocity_[neighborIdx_[elementIdx][dimIdx][1]][velIdx]
                              - velocity_[neighborIdx_[elementIdx][dimIdx][0]][velIdx])
                          / (cellCenter_[neighborIdx_[elementIdx][dimIdx][1]][dimIdx]
                              - cellCenter_[neighborIdx_[elementIdx][dimIdx][0]][dimIdx]);
                    if (abs(cellCenter_[neighborIdx_[elementIdx][dimIdx][1]][dimIdx]
                            - cellCenter_[neighborIdx_[elementIdx][dimIdx][0]][dimIdx]) < 1e-8)
                        velocityGradients_[elementIdx][velIdx][dimIdx] = 0.0;
                }

                if (abs(velocity_[elementIdx][dimIdx]) > abs(velocityMaximum_[wallElementIdx][dimIdx]))
                {
                    velocityMaximum_[wallElementIdx][dimIdx] = velocity_[elementIdx][dimIdx];
                }
                if (abs(velocity_[elementIdx][dimIdx]) < abs(velocityMinimum_[wallElementIdx][dimIdx]))
                {
                    velocityMinimum_[wallElementIdx][dimIdx] = velocity_[elementIdx][dimIdx];
                }

                if (0 <= flowNormalAxis && flowNormalAxis < dim)
                {
                    flowNormalAxis_[elementIdx] = flowNormalAxis;
                }
                else if (abs(maxVelocity) < abs(velocity_[elementIdx][dimIdx]))
                {
                    maxVelocity = abs(velocity_[elementIdx][dimIdx]);
                    flowNormalAxis_[elementIdx] = dimIdx;
                }
            }

            auto fvGeometry = localView(this->fvGridGeometry());
            fvGeometry.bindElement(element);
            for (auto&& scvf : scvfs(fvGeometry))
            {
                // adapt calculations for Dirichlet condition
                unsigned int scvfNormDim = scvf.directionIndex();
                if (scvf.boundary())
                {
                    for (unsigned int velIdx = 0; velIdx < dim; ++velIdx)
                    {
                        if (!asImp_().boundaryTypes(element, scvf).isDirichlet(Indices::velocity(velIdx)))
                            continue;

                        Scalar dirichletVelocity = asImp_().dirichlet(element, scvf)[Indices::velocity(velIdx)];

                        unsigned int neighborIdx = neighborIdx_[elementIdx][scvfNormDim][0];
                        if (scvf.center()[scvfNormDim] < cellCenter_[elementIdx][scvfNormDim])
                            neighborIdx = neighborIdx_[elementIdx][scvfNormDim][1];

                        velocityGradients_[elementIdx][velIdx][scvfNormDim]
                            = (velocity_[neighborIdx][velIdx] - dirichletVelocity)
                              / (cellCenter_[neighborIdx][scvfNormDim] - scvf.center()[scvfNormDim]);
                    }
                }

                // Calculate the BJS-velocity by accounting for all sub faces.
                std::vector<int> bjsNumFaces(dim, 0);
                std::vector<unsigned int> bjsNeighbor(dim, 0);
                DimVector bjsVelocityAverage(0.0);
                DimVector normalNormCoordinate(0.0);
                unsigned int velIdx = Indices::velocity(scvfNormDim);
                const int numSubFaces = scvf.pairData().size();
                for(int localSubFaceIdx = 0; localSubFaceIdx < numSubFaces; ++localSubFaceIdx)
                {
                    const auto& normalFace = fvGeometry.scvf(scvf.insideScvIdx(), scvf.pairData()[localSubFaceIdx].localNormalFaceIdx);

                    // adapt calculations for Beavers-Joseph-Saffman condition
                    unsigned int normalNormDim = normalFace.directionIndex();
                    if (normalFace.boundary() && (asImp_().boundaryTypes(element, normalFace).isBJS(Indices::velocity(velIdx))))
                    {
                        unsigned int neighborIdx = neighborIdx_[elementIdx][normalNormDim][0];
                        if (normalFace.center()[normalNormDim] < cellCenter_[elementIdx][normalNormDim])
                            neighborIdx = neighborIdx_[elementIdx][normalNormDim][1];

                        bjsVelocityAverage[normalNormDim] += ParentType::bjsVelocity(scvf, normalFace, localSubFaceIdx, velocity_[elementIdx][velIdx]);
                        if (bjsNumFaces[normalNormDim] > 0 && neighborIdx != bjsNeighbor[normalNormDim])
                            DUNE_THROW(Dune::InvalidStateException, "Two different neighborIdx should not occur");
                        bjsNeighbor[normalNormDim] = neighborIdx;
                        normalNormCoordinate[normalNormDim] = normalFace.center()[normalNormDim];
                        bjsNumFaces[normalNormDim]++;
                    }
                }
                for (unsigned dirIdx = 0; dirIdx < dim; ++dirIdx)
                {
                    if (bjsNumFaces[dirIdx] == 0)
                        continue;

                    unsigned int neighborIdx = bjsNeighbor[dirIdx];
                    bjsVelocityAverage[dirIdx] /= bjsNumFaces[dirIdx];

                    velocityGradients_[elementIdx][velIdx][dirIdx]
                        = (velocity_[neighborIdx][velIdx] - bjsVelocityAverage[dirIdx])
                          / (cellCenter_[neighborIdx][dirIdx] - normalNormCoordinate[dirIdx]);

                }
            }
        }

        // calculate or call all secondary variables
        for (const auto& element : elements(this->fvGridGeometry().gridView()))
        {
            unsigned int elementIdx = this->fvGridGeometry().elementMapper().index(element);

            Dune::FieldMatrix<Scalar, GridView::dimension, GridView::dimension> stressTensor(0.0);
            for (unsigned int dimIdx = 0; dimIdx < dim; ++dimIdx)
            {
                for (unsigned int velIdx = 0; velIdx < dim; ++velIdx)
                {
                    stressTensor[dimIdx][velIdx] = 0.5 * velocityGradients_[elementIdx][dimIdx][velIdx]
                                                   + 0.5 * velocityGradients_[elementIdx][velIdx][dimIdx];
              }
            }
            stressTensorScalarProduct_[elementIdx] = 0.0;
            for (unsigned int dimIdx = 0; dimIdx < dim; ++dimIdx)
            {
                for (unsigned int velIdx = 0; velIdx < dim; ++velIdx)
                {
                    stressTensorScalarProduct_[elementIdx] += stressTensor[dimIdx][velIdx] * stressTensor[dimIdx][velIdx];
                }
            }

            Dune::FieldMatrix<Scalar, GridView::dimension, GridView::dimension> vorticityTensor(0.0);
            for (unsigned int dimIdx = 0; dimIdx < dim; ++dimIdx)
            {
                for (unsigned int velIdx = 0; velIdx < dim; ++velIdx)
                {
                    vorticityTensor[dimIdx][velIdx] = 0.5 * velocityGradients_[elementIdx][dimIdx][velIdx]
                                                      - 0.5 * velocityGradients_[elementIdx][velIdx][dimIdx];
              }
            }
            vorticityTensorScalarProduct_[elementIdx] = 0.0;
            for (unsigned int dimIdx = 0; dimIdx < dim; ++dimIdx)
            {
                for (unsigned int velIdx = 0; velIdx < dim; ++velIdx)
                {
                    vorticityTensorScalarProduct_[elementIdx] += vorticityTensor[dimIdx][velIdx] * vorticityTensor[dimIdx][velIdx];
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
                kinematicViscosity_[elementIdx] = volVars.viscosity() / volVars.density();
            }
        }
    }

    /*!
     * \brief Returns whether a given sub control volume face is on a wall
     *
     * \param scvf The sub control volume face.
     */
    bool isOnWall(const SubControlVolumeFace& scvf) const
    {
        return asImp_().isOnWallAtPos(scvf.center());
    }

    /*!
     * \brief Returns whether a given point is on a wall
     *
     * \param globalPos The position in global coordinates.
     */
    bool isOnWallAtPos(const GlobalPosition &globalPos) const
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

    /*!
     * \brief Return the turbulent Prandtl number \f$ [-] \f$ which is used to convert
     *        the eddy viscosity to an eddy thermal conductivity
     */
    Scalar turbulentPrandtlNumber() const
    {
        static const Scalar turbulentPrandtlNumber
            = getParamFromGroup<Scalar>(this->paramGroup(), "RANS.TurbulentPrandtlNumber", 1.0);
        return turbulentPrandtlNumber;
    }

    /*!
     * \brief Return the turbulent Schmidt number \f$ [-] \f$ which is used to convert
     *        the eddy viscosity to an eddy diffusivity
     */
    Scalar turbulentSchmidtNumber() const
    {
        static const Scalar turbulentSchmidtNumber
            = getParamFromGroup<Scalar>(this->paramGroup(), "RANS.TurbulentSchmidtNumber", 1.0);
        return turbulentSchmidtNumber;
    }

public:
    bool calledUpdateStaticWallProperties = false;
    std::vector<unsigned int> wallElementIdx_;
    std::vector<Scalar> wallDistance_;
    std::vector<std::array<std::array<unsigned int, 2>, dim>> neighborIdx_;
    std::vector<GlobalPosition> cellCenter_;
    std::vector<DimVector> velocity_;
    std::vector<DimVector> velocityMaximum_;
    std::vector<DimVector> velocityMinimum_;
    std::vector<DimMatrix> velocityGradients_;
    std::vector<Scalar> stressTensorScalarProduct_;
    std::vector<Scalar> vorticityTensorScalarProduct_;
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
