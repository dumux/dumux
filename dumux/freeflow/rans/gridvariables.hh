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
 * \ingroup StaggeredDiscretization
 * \copydoc Dumux::StaggeredGridVariables
 */
#ifndef DUMUX_RANS_GRID_VARIABLES_HH
#define DUMUX_RANS_GRID_VARIABLES_HH

#include <dumux/discretization/staggered/gridvariables.hh>
#include <dune/common/fmatrix.hh>
#include <dune/common/float_cmp.hh>

namespace Dumux {

/*!
 * \ingroup StaggeredDiscretization
 * \brief Class storing data associated to scvs and scvfs
 * \tparam GG the type of the grid geometry
 * \tparam GVV the type of the grid volume variables
 * \tparam GFVC the type of the grid flux variables cache
 * \tparam GFV the type of the grid face variables
 */
template<class NavierStokesGridVariables, class Implementation>
class RANSGridVariables : public NavierStokesGridVariables
{
    using ParentType = NavierStokesGridVariables;
    using FVGridGeometry = typename ParentType::GridGeometry;
    using Element = typename FVGridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    static constexpr auto dim = FVGridGeometry::GridView::dimension;

    using DimVector = GlobalPosition;
    using DimMatrix = Dune::FieldMatrix<typename ParentType::Scalar, dim, dim>;

    using Indices = typename ParentType::VolumeVariables::Indices;

    // using ThisType = StaggeredGridVariables<GG, GVV, GFVC, GFV>;
    // friend class StaggeredGridVariablesView<ThisType>;

    // static constexpr auto cellCenterIdx = FVGridGeometry::cellCenterIdx();
    // static constexpr auto faceIdx = FVGridGeometry::faceIdx();

public:
    //! Export the Scalar type
    using Scalar = typename ParentType::Scalar;

    // using ParentType::ParentType;

    //! Constructor
    template<class Problem>
    RANSGridVariables(std::shared_ptr<Problem> problem,
                      std::shared_ptr<FVGridGeometry> fvGridGeometry)
    : ParentType(problem, fvGridGeometry)
    {
        asImp_().updateStaticWallProperties();
    }

    //! update all variables
    template<class SolutionVector>
    void update(const SolutionVector& curSol)
    {
        ParentType::update(curSol);
        asImp_().updateDynamicWallProperties(curSol);
    }

    //! initialize all variables (stationary case)
    template<class SolutionVector>
    void init(const SolutionVector& curSol)
    {
        ParentType::init(curSol);
        asImp_().updateDynamicWallProperties(curSol);
    }

    //! initialize all variables (instationary case)
    template<class SolutionVector>
    void init(const SolutionVector& curSol, const SolutionVector& initSol)
    {
        ParentType::init(curSol, initSol);
        asImp_().updateDynamicWallProperties(curSol);
    }

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

        const auto& problem = this->curGridVolVars().problem();

        const auto numElements =  this->fvGridGeometry().gridView().size(0);

        // update size and initial values of the global vectors
        correspondingWallElementIdx_.resize(numElements);
        wallDistance_.resize(numElements, std::numeric_limits<Scalar>::max());
        neighborElementIdx_.resize(numElements);
        cellCenter_.resize(numElements, GlobalPosition(0.0));
        velocity_.resize(numElements, DimVector(0.0));
        velocityMaximum_.resize(numElements, DimVector(0.0));
        velocityGradients_.resize(numElements, DimMatrix(0.0));
        stressTensorScalarProduct_.resize(numElements, 0.0);
        vorticityTensorScalarProduct_.resize(numElements, 0.0);
        flowNormalAxis_.resize(numElements, 0);
        wallNormalAxis_.resize(numElements, 1);
        kinematicViscosity_.resize(numElements, 0.0);
        sandGrainRoughness_.resize(numElements, 0.0);

        // retrieve all wall intersections and corresponding elements
        std::vector<std::size_t> wallElementIdx;
        std::vector<GlobalPosition> wallPositions;
        std::vector<std::size_t> wallNormalAxisTemp;

        // reserve some memory for the vectors for increased performance
        const auto numBoundaryScvf = this->fvGridGeometry().numBoundaryScvf();
        wallElementIdx.reserve(numBoundaryScvf);
        wallPositions.reserve(numBoundaryScvf);
        wallNormalAxisTemp.reserve(numBoundaryScvf);

        const auto gridView = this->fvGridGeometry().gridView();
        auto fvGeometry = localView(this->fvGridGeometry());

        // move along the domains boundary and find elements and scvfs adjacent to the wall
        for (const auto& element : elements(gridView))
        {
            fvGeometry.bindElement(element);
            for (const auto& scvf : scvfs(fvGeometry))
            {
                // only search for walls at a global boundary
                if (!scvf.boundary())
                    continue;

                if (problem.isOnWall(scvf))
                {
                    wallElementIdx.push_back(this->fvGridGeometry().elementMapper().index(element));
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
            std::size_t elementIdx = this->fvGridGeometry().elementMapper().index(element);
            cellCenter_[elementIdx] = element.geometry().center();
            for (std::size_t i = 0; i < wallPositions.size(); ++i)
            {
                static const int problemWallNormalAxis
                    = getParamFromGroup<int>(problem.paramGroup(), "RANS.WallNormalAxis", -1);
                int searchAxis = problemWallNormalAxis;

                // search along wall normal axis of the intersection
                if (problemWallNormalAxis < 0 || problemWallNormalAxis >= dim)
                    searchAxis = wallNormalAxisTemp[i];

                const GlobalPosition distanceToWall = cellCenter_[elementIdx] - wallPositions[i];
                // find minimum wall distance and make sure the search axis is honored
                if (abs(distanceToWall[searchAxis]) < wallDistance_[elementIdx]
                    && Dune::FloatCmp::eq(distanceToWall[searchAxis], distanceToWall.two_norm())) // TODO: check if this works
                {
                    wallDistance_[elementIdx] = abs(distanceToWall[searchAxis]);
                    correspondingWallElementIdx_[elementIdx] = wallElementIdx[i];
                    wallNormalAxis_[elementIdx] = searchAxis;
                    sandGrainRoughness_[elementIdx] = problem.sandGrainRoughnessAtPos(wallPositions[i]);
                }
            }
        }

        // search for neighbor Idxs
        for (const auto& element : elements(gridView))
        {
            std::size_t elementIdx = this->fvGridGeometry().elementMapper().index(element);
            for (std::size_t dimIdx = 0; dimIdx < dim; ++dimIdx)
            {
                neighborElementIdx_[elementIdx][dimIdx][0] = elementIdx;
                neighborElementIdx_[elementIdx][dimIdx][1] = elementIdx;
            }

            for (const auto& intersection : intersections(gridView, element))
            {
                if (intersection.boundary())
                    continue;

                std::size_t neighborElementIdx = this->fvGridGeometry().elementMapper().index(intersection.outside());
                for (std::size_t dimIdx = 0; dimIdx < dim; ++dimIdx)
                {
                    if (abs(cellCenter_[elementIdx][dimIdx] - cellCenter_[neighborElementIdx][dimIdx]) > 1e-8) // TODO: required?
                    {
                        if (cellCenter_[elementIdx][dimIdx] > cellCenter_[neighborElementIdx][dimIdx])
                            neighborElementIdx_[elementIdx][dimIdx][0] = neighborElementIdx;
                        if (cellCenter_[elementIdx][dimIdx] < cellCenter_[neighborElementIdx][dimIdx])
                            neighborElementIdx_[elementIdx][dimIdx][1] = neighborElementIdx;
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
    template<class SolutionVector>
    void updateDynamicWallProperties(const SolutionVector& curSol)
    {
        using std::abs;
        using std::max;
        using std::min;
        std::cout << "Update dynamic wall properties." << std::endl;
        if (!calledUpdateStaticWallProperties)
            DUNE_THROW(Dune::InvalidStateException,
                       "You have to call updateStaticWallProperties() once before you call updateDynamicWallProperties().");

        const auto& problem = this->curGridVolVars().problem();
        auto fvGeometry = localView(this->fvGridGeometry());

        static const int flowNormalAxis
            = getParamFromGroup<int>(problem.paramGroup(), "RANS.FlowNormalAxis", -1);

        // re-initialize min and max values
        std::fill(velocityMaximum_.begin(), velocityMaximum_.end(), DimVector(std::numeric_limits<Scalar>::min()));
        std::fill(velocityMinimum_.begin(), velocityMinimum_.end(), DimVector(std::numeric_limits<Scalar>::max()));

        // calculate cell-center-averaged velocities
        for (const auto& element : elements(this->fvGridGeometry().gridView()))
        {
            fvGeometry.bindElement(element);
            const auto elementIdx = this->fvGridGeometry().elementMapper().index(element);

            // calculate velocities
            DimVector velocityTemp(0.0);
            for (auto&& scvf : scvfs(fvGeometry))
            {
                const auto numericalSolutionFace = curSol[FVGridGeometry::faceIdx()][scvf.dofIndex()][Indices::velocity(scvf.directionIndex())];
                velocityTemp[scvf.directionIndex()] += numericalSolutionFace;
            }
            for (std::size_t dimIdx = 0; dimIdx < dim; ++dimIdx)
                velocity_[elementIdx][dimIdx] = velocityTemp[dimIdx] * 0.5; // faces are equidistant to cell center
        }

        // calculate cell-center-averaged velocity gradients, maximum, and minimum values
        for (const auto& element : elements(this->fvGridGeometry().gridView()))
        {
            const auto elementIdx = this->fvGridGeometry().elementMapper().index(element);
            const auto correspondingWallElementIdx = correspondingWallElementIdx_[elementIdx];

            Scalar maxVelocity = 0.0;
            for (std::size_t dimIdx = 0; dimIdx < dim; ++dimIdx)
            {
                for (std::size_t velIdx = 0; velIdx < dim; ++velIdx)
                {
                    const auto forwardElementIdx = neighborElementIdx_[elementIdx][dimIdx][1];
                    const auto backWardElementIdx = neighborElementIdx_[elementIdx][dimIdx][0];

                    velocityGradients_[elementIdx][velIdx][dimIdx]
                        = (velocity_[forwardElementIdx][velIdx]
                              - velocity_[backWardElementIdx][velIdx])
                          / (cellCenter_[forwardElementIdx][dimIdx]
                              - cellCenter_[backWardElementIdx][dimIdx]);
                    if (abs(cellCenter_[forwardElementIdx][dimIdx]
                            - cellCenter_[backWardElementIdx][dimIdx]) < 1e-8)
                        velocityGradients_[elementIdx][velIdx][dimIdx] = 0.0; // TODO check via indices
                }

                if (abs(velocity_[elementIdx][dimIdx]) > abs(velocityMaximum_[correspondingWallElementIdx][dimIdx]))
                {
                    velocityMaximum_[correspondingWallElementIdx][dimIdx] = velocity_[elementIdx][dimIdx];
                }
                if (abs(velocity_[elementIdx][dimIdx]) < abs(velocityMinimum_[correspondingWallElementIdx][dimIdx]))
                {
                    velocityMinimum_[correspondingWallElementIdx][dimIdx] = velocity_[elementIdx][dimIdx];
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
    //
    //         auto fvGeometry = localView(this->fvGridGeometry());
    //         fvGeometry.bindElement(element);
    //         for (auto&& scvf : scvfs(fvGeometry))
    //         {
    //             // adapt calculations for Dirichlet condition
    //             std::size_t scvfNormDim = scvf.directionIndex();
    //             if (scvf.boundary())
    //             {
    //                 for (std::size_t velIdx = 0; velIdx < dim; ++velIdx)
    //                 {
    //                     if (!asImp_().boundaryTypes(element, scvf).isDirichlet(Indices::velocity(velIdx)))
    //                         continue;
    //
    //                     Scalar dirichletVelocity = asImp_().dirichlet(element, scvf)[Indices::velocity(velIdx)];
    //
    //                     std::size_t neighborIdx = neighborIdx_[elementIdx][scvfNormDim][0];
    //                     if (scvf.center()[scvfNormDim] < cellCenter_[elementIdx][scvfNormDim])
    //                         neighborIdx = neighborIdx_[elementIdx][scvfNormDim][1];
    //
    //                     velocityGradients_[elementIdx][velIdx][scvfNormDim]
    //                         = (velocity_[neighborIdx][velIdx] - dirichletVelocity)
    //                           / (cellCenter_[neighborIdx][scvfNormDim] - scvf.center()[scvfNormDim]);
    //                 }
    //             }
    //
    //             // Calculate the BJS-velocity by accounting for all sub faces.
    //             std::vector<int> bjsNumFaces(dim, 0);
    //             std::vector<std::size_t> bjsNeighbor(dim, 0);
    //             DimVector bjsVelocityAverage(0.0);
    //             DimVector normalNormCoordinate(0.0);
    //             std::size_t velIdx = Indices::velocity(scvfNormDim);
    //             const int numSubFaces = scvf.pairData().size();
    //             for(int localSubFaceIdx = 0; localSubFaceIdx < numSubFaces; ++localSubFaceIdx)
    //             {
    //                 const auto& normalFace = fvGeometry.scvf(scvf.insideScvIdx(), scvf.pairData()[localSubFaceIdx].localNormalFaceIdx);
    //
    //                 // adapt calculations for Beavers-Joseph-Saffman condition
    //                 std::size_t normalNormDim = normalFace.directionIndex();
    //                 if (normalFace.boundary() && (asImp_().boundaryTypes(element, normalFace).isBJS(Indices::velocity(velIdx))))
    //                 {
    //                     std::size_t neighborIdx = neighborIdx_[elementIdx][normalNormDim][0];
    //                     if (normalFace.center()[normalNormDim] < cellCenter_[elementIdx][normalNormDim])
    //                         neighborIdx = neighborIdx_[elementIdx][normalNormDim][1];
    //
    //                     bjsVelocityAverage[normalNormDim] += ParentType::bjsVelocity(scvf, normalFace, localSubFaceIdx, velocity_[elementIdx][velIdx]);
    //                     if (bjsNumFaces[normalNormDim] > 0 && neighborIdx != bjsNeighbor[normalNormDim])
    //                         DUNE_THROW(Dune::InvalidStateException, "Two different neighborIdx should not occur");
    //                     bjsNeighbor[normalNormDim] = neighborIdx;
    //                     normalNormCoordinate[normalNormDim] = normalFace.center()[normalNormDim];
    //                     bjsNumFaces[normalNormDim]++;
    //                 }
    //             }
    //             for (unsigned dirIdx = 0; dirIdx < dim; ++dirIdx)
    //             {
    //                 if (bjsNumFaces[dirIdx] == 0)
    //                     continue;
    //
    //                 std::size_t neighborIdx = bjsNeighbor[dirIdx];
    //                 bjsVelocityAverage[dirIdx] /= bjsNumFaces[dirIdx];
    //
    //                 velocityGradients_[elementIdx][velIdx][dirIdx]
    //                     = (velocity_[neighborIdx][velIdx] - bjsVelocityAverage[dirIdx])
    //                       / (cellCenter_[neighborIdx][dirIdx] - normalNormCoordinate[dirIdx]);
    //
    //             }
    //         }
        }
    //
    //     // calculate or call all secondary variables
    //     for (const auto& element : elements(this->fvGridGeometry().gridView()))
    //     {
    //         std::size_t elementIdx = this->fvGridGeometry().elementMapper().index(element);
    //
    //         Dune::FieldMatrix<Scalar, GridView::dimension, GridView::dimension> stressTensor(0.0);
    //         for (std::size_t dimIdx = 0; dimIdx < dim; ++dimIdx)
    //         {
    //             for (std::size_t velIdx = 0; velIdx < dim; ++velIdx)
    //             {
    //                 stressTensor[dimIdx][velIdx] = 0.5 * velocityGradients_[elementIdx][dimIdx][velIdx]
    //                                                + 0.5 * velocityGradients_[elementIdx][velIdx][dimIdx];
    //           }
    //         }
    //         stressTensorScalarProduct_[elementIdx] = 0.0;
    //         for (std::size_t dimIdx = 0; dimIdx < dim; ++dimIdx)
    //         {
    //             for (std::size_t velIdx = 0; velIdx < dim; ++velIdx)
    //             {
    //                 stressTensorScalarProduct_[elementIdx] += stressTensor[dimIdx][velIdx] * stressTensor[dimIdx][velIdx];
    //             }
    //         }
    //
    //         Dune::FieldMatrix<Scalar, GridView::dimension, GridView::dimension> vorticityTensor(0.0);
    //         for (std::size_t dimIdx = 0; dimIdx < dim; ++dimIdx)
    //         {
    //             for (std::size_t velIdx = 0; velIdx < dim; ++velIdx)
    //             {
    //                 vorticityTensor[dimIdx][velIdx] = 0.5 * velocityGradients_[elementIdx][dimIdx][velIdx]
    //                                                   - 0.5 * velocityGradients_[elementIdx][velIdx][dimIdx];
    //           }
    //         }
    //         vorticityTensorScalarProduct_[elementIdx] = 0.0;
    //         for (std::size_t dimIdx = 0; dimIdx < dim; ++dimIdx)
    //         {
    //             for (std::size_t velIdx = 0; velIdx < dim; ++velIdx)
    //             {
    //                 vorticityTensorScalarProduct_[elementIdx] += vorticityTensor[dimIdx][velIdx] * vorticityTensor[dimIdx][velIdx];
    //             }
    //         }
    //
    //         auto fvGeometry = localView(this->fvGridGeometry());
    //         fvGeometry.bindElement(element);
    //         for (auto&& scv : scvs(fvGeometry))
    //         {
    //             const int dofIdx = scv.dofIndex();
    //
    //             // construct a privars object from the cell center solution vector
    //             const auto& cellCenterPriVars = curSol[FVGridGeometry::cellCenterIdx()][dofIdx];
    //             PrimaryVariables priVars = makePriVarsFromCellCenterPriVars<PrimaryVariables>(cellCenterPriVars);
    //             auto elemSol = elementSolution<typename FVGridGeometry::LocalView>(std::move(priVars));
    //
    //             VolumeVariables volVars;
    //             volVars.update(elemSol, asImp_(), element, scv);
    //             kinematicViscosity_[elementIdx] = volVars.viscosity() / volVars.density();
    //         }
        // }
    }
protected:
    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }

private:
    bool calledUpdateStaticWallProperties = false;
    std::vector<std::size_t> correspondingWallElementIdx_;
    std::vector<Scalar> wallDistance_;
    std::vector<std::array<std::array<std::size_t, 2>, dim>> neighborElementIdx_;
    std::vector<GlobalPosition> cellCenter_;
    std::vector<DimVector> velocity_;
    std::vector<DimVector> velocityMaximum_;
    std::vector<DimVector> velocityMinimum_;
    std::vector<DimMatrix> velocityGradients_;
    std::vector<Scalar> stressTensorScalarProduct_;
    std::vector<Scalar> vorticityTensorScalarProduct_;
    std::vector<std::size_t> wallNormalAxis_;
    std::vector<std::size_t> flowNormalAxis_;
    std::vector<Scalar> kinematicViscosity_;
    std::vector<Scalar> sandGrainRoughness_;

};

} // end namespace Dumux

#endif
