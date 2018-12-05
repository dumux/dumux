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
#include <dumux/discretization/method.hh>
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
    using Implementation = GetPropType<TypeTag, Properties::Problem>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    using FVGridGeometry = GetPropType<TypeTag, Properties::FVGridGeometry>;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using GridView = typename FVGridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using PrimaryVariables = typename VolumeVariables::PrimaryVariables;
    using CellCenterPrimaryVariables = GetPropType<TypeTag, Properties::CellCenterPrimaryVariables>;
    using FacePrimaryVariables = GetPropType<TypeTag, Properties::FacePrimaryVariables>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

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
    {
        // update size and initial values of the global vectors
        int numElements = this->fvGridGeometry().elementMapper().size();
        cellCenter_.resize(numElements, GlobalPosition(0.0));
        wallElementIdx_.resize(numElements);
        wallDistance_.resize(numElements, std::numeric_limits<Scalar>::max());
        wallProfileIdx_.resize(numElements, 0);
        neighborIdx_.resize(numElements);
        velocity_.resize(numElements, DimVector(0.0));
        velocityGradients_.resize(numElements, DimMatrix(0.0));
        stressTensorScalarProduct_.resize(numElements, 0.0);
        vorticityTensorScalarProduct_.resize(numElements, 0.0);
        flowNormalAxis_.resize(numElements, 0);
        wallNormalAxis_.resize(numElements, 1);
        kinematicViscosity_.resize(numElements, 0.0);
        sandGrainRoughness_.resize(numElements, 0.0);
    }

    /*!
     * \brief Update the static (solution independent) relations to the walls
     *
     * This function first stores the location of each element (storeCCPositions),
     * then stores a list of all elements with a wall intersection (findAllWallElements),
     * then stores the distance from the wall to the element center for each element (storeWallInformation),
     * and finally stores the indexes of the neighboring elements for each element (storeNeighborCellInformation).
     */
    void updateStaticWallProperties()
    {
        std::cout << "Update static wall properties. ";
        calledUpdateStaticWallProperties = true;

        // update the the cell centered locations
        storeCCPositions(this->fvGridGeometry().gridView());

        // fill the vector with elements on the wall
        findAllWallElements(this->fvGridGeometry().gridView());

        std::cout << "NumWallIntersections=" << wallPositions_.size() << std::endl;
        if (wallPositions_.size() == 0)
            DUNE_THROW(Dune::InvalidStateException,
                       "No wall intersections have been found. Make sure that the isOnWall(globalPos) is working properly.");

        velocityMaximum_.resize(wallPositions_.size(), DimVector(1e-16));
        velocityMinimum_.resize(wallPositions_.size(), DimVector(std::numeric_limits<Scalar>::max()));

        // store the distance from the wall to the element
        storeWallInformation(this->fvGridGeometry().gridView());

        // stores the indexes of the neighboring elements for each element
        storeNeighborCellInformation(this->fvGridGeometry().gridView());

        std::cout << "\n";
        std::cout << "wallElements_ size: "<< wallElements_.size() <<"\n";
        for (int i = 0; i < wallElements_.size(); i++)
            std::cout << wallElements_[i] << "\n";
        std::cout << "\n";

        std::cout << "\n";
        std::cout << "wallPositions_ size: "<< wallPositions_.size() <<"\n";
        for (int i = 0; i < wallPositions_.size(); i++)
            std::cout << wallPositions_[i] << "\n";
        std::cout << "\n";

        std::cout << "\n";
        std::cout << "wallElementIdx_ size: "<< wallElementIdx_.size() <<"\n";
        for (int i = 0; i < wallElementIdx_.size(); i++)
            std::cout << wallElementIdx_[i] << "\n";
        std::cout << "\n";

        std::cout << "\n";
        std::cout << "wallProfileIdx_ size: "<< wallProfileIdx_.size() <<"\n";
        for (int i = 0; i < wallProfileIdx_.size(); i++)
            std::cout << wallProfileIdx_[i] << "\n";
        std::cout << "\n";

    }

    /*!
     * \brief Update the dynamic (solution dependent) relations to the walls
     *
     * This function first stores the cell centered velocities (storeCCVelocites),
     * then stores the maximum and minimum velocities located along the axis perpendicular to the wall (storeCCmaxMinVelocities),
     * then stores the velocity gradients (storeCCVelocityGradients),
     * and then calculates the stress tensor and vorticity tensors (calculateStressTensorInformation, calculateVorticityTensorInformation).
     * Further, the kinematic viscosity at the wall is stored.
     *
     * \param curSol The solution vector.
     */
    void updateDynamicWallProperties(const SolutionVector& curSol)
    {
        std::cout << "Update dynamic wall properties." << std::endl;
        if (!calledUpdateStaticWallProperties)
            DUNE_THROW(Dune::InvalidStateException,
                       "You have to call updateStaticWallProperties() once before you call updateDynamicWallProperties().");

        // store the cell centered velocities
        storeCCVelocites(this->fvGridGeometry().gridView(), curSol);

        // store the maximum and minimum velocities located along the axis perpendicular to the wall
        storeCCMaxMinVelocities(this->fvGridGeometry().gridView());

        // calculate and store the velocity gradients
        storeCCVelocityGradients(this->fvGridGeometry().gridView());

        // calculate and store the stress tensor and tensor product
        calculateStressTensorInformation(this->fvGridGeometry().gridView());

        // calculate and store the vortuosity tensor and tensor product
        calculateVorticityTensorInformation(this->fvGridGeometry().gridView());

        // store the kinematic viscosity at the wall
        updateKinematicViscosity(this->fvGridGeometry().gridView(), curSol);
    }

    /*!
     * \brief Stores all cell centered global wallPositions_
     *
     *  Fills a vector of cell centered locations indexed with the elementIdx
     *
     * \param gridView the grid view we are solving on
     */
    void storeCCPositions(const GridView& gridView)
    {
        // search for shortest distance to wall for each element
        for (const auto& element : elements(gridView))
        {
            cellCenter_[this->fvGridGeometry().elementMapper().index(element)] = element.geometry().center();
        }
    }

    /*!
     * \brief Find all of the elements will a wall intersection
     *
     *  Stores all of the elements with a wall intersection (wallElements_),
     *  their position (wallPositions_),
     *  and the direction of the normal vector (wallNormalAxisTemp_).
     *
     * \param gridView the grid view we are solving on
     */
    void findAllWallElements(const GridView& gridView)
    {
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
                    wallElements_.push_back(this->fvGridGeometry().elementMapper().index(element));
                    wallPositions_.push_back(scvf.center());
                    wallNormalAxisTemp_.push_back(scvf.directionIndex());
                }
            }
        }
    }

    /*!
     * \brief Store the distance from the wall to the element center and further wall information
     *
     *  A vector wallDistance_ is filled with the distance from the cell center to the wall for each element
     *  A vector wallElementIdx_ is filled with the index of the element at the closest wallClockTime
     *  A vector wallNormalAxis_ is filled with the
     *
     * \param gridView the grid view we are solving on
     */
    void storeWallInformation(const GridView& gridView)
    {
        using std::abs;
        // search for shortest distance to wall for each element
        for (const auto& element : elements(gridView))
        {
            unsigned int elementIdx = this->fvGridGeometry().elementMapper().index(element);
            for (unsigned int i = 0; i < wallPositions_.size(); ++i)
            {
                int searchAxis = wallNormalAxisTemp_[i];

                GlobalPosition global = element.geometry().center();
                global -= wallPositions_[i];
                // second and argument ensures to use only aligned elements
                if (abs(global[searchAxis]) < wallDistance_[elementIdx]
                    && abs(global[searchAxis]) < global.two_norm() + 1e-8
                    && abs(global[searchAxis]) > global.two_norm() - 1e-8)
                {
                    wallDistance_[elementIdx] = abs(global[searchAxis]);
                    wallElementIdx_[elementIdx] = wallElements_[i];
                    wallProfileIdx_[elementIdx]= i;
                    wallNormalAxis_[elementIdx] = searchAxis;
                    sandGrainRoughness_[elementIdx] = asImp_().sandGrainRoughnessAtPos(wallPositions_[i]);
                }
            }
        }
    }

    /*!
     * \brief For each element, the indicies of the neighboring elements are stored
     *
     *  A vector of matrices(dimx2) called neighborIdx is filled.
     *  Each entry in this vector will represent an element, indexed with elementIdx.
     *  Each entry will have two neighboring indexes per dimension,
     *  one for the element further from the origin, one for the element closer to the origin
     *
     * \param gridView the grid view we are solving on
     */
    void storeNeighborCellInformation(const GridView& gridView)
    {
        using std::abs;
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
     * \brief Stores the cell centered velocity
     *
     *  Builds a vector velocity_ where for each element,
     *  the average of the velocities on the faces is stored at the cell center
     *
     * \param gridView the grid view we are solving on
     * \param curSol The solution vector.
     */
    void storeCCVelocites(const GridView& gridView, const SolutionVector& curSol)
    {
        // calculate cell-center-averaged velocities
        for (const auto& element : elements(gridView))
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
    }

    /*!
     * \brief Stores the maxiumum and minimum velocities
     *
     *  For each profile perpendicular to a wall position, the maxmimum and minimum velocities are logged.
     *
     * \param gridView the grid view we are solving on
     */
    void storeCCMaxMinVelocities(const GridView& gridView)
    {
        using std::abs;
        using std::max;
        using std::min;

        // re-initialize min and max values
        std::fill(velocityMaximum_.begin(), velocityMaximum_.end(), DimVector(1e-16));
        std::fill(velocityMinimum_.begin(), velocityMinimum_.end(), DimVector(std::numeric_limits<Scalar>::max()));

        // calculate the maximum and minumum velocities in each profile
        for (const auto& element : elements(gridView))
        {
            unsigned int elementIdx = this->fvGridGeometry().elementMapper().index(element);
            unsigned int profileIdx = this->wallProfileIdx_[elementIdx];
            for (unsigned int dimIdx = 0; dimIdx < dim; ++dimIdx)
            {
                if (abs(this->velocity_[elementIdx][dimIdx]) > abs(velocityMaximum_[profileIdx][dimIdx]))
                {
                    velocityMaximum_[profileIdx][dimIdx] = this->velocity_[elementIdx][dimIdx];
                }
                if (abs(this->velocity_[elementIdx][dimIdx]) < abs(velocityMinimum_[profileIdx][dimIdx]))
                {
                    velocityMinimum_[profileIdx][dimIdx] = this->velocity_[elementIdx][dimIdx];
                }
            }
        }
    }

    /*!
     * \brief Calculates the velocity gradients in each element
     *
     *  For each cell, a velocity gradient matrix is filled.
     *  In the case of a dirichlet boundary, the velocity gradient matrix is adapted.
     *  In the case of a Beavers Joeseph conditon at the boundary, the matrix is adapted.
     *
     * \param gridView the grid view we are solving on
     */
    void storeCCVelocityGradients(const GridView& gridView)
    {
        using std::abs;
        // calculate cell-center-averaged velocity gradients
        for (const auto& element : elements(gridView))
        {
            unsigned int elementIdx = this->fvGridGeometry().elementMapper().index(element);
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

                if (abs(maxVelocity) < abs(velocity_[elementIdx][dimIdx]))
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

                        bjsVelocityAverage[normalNormDim] += ParentType::bjsVelocity(element, scvf, normalFace, localSubFaceIdx, velocity_[elementIdx][velIdx]);
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
    }

    /*!
     * \brief Calculates the Stress Tensor Information in each element
     *
     *  First the Stress Tensor matrix is filled,
     *  then the scalar product of this matrix is calculated.
     *
     * \param gridView the grid view we are solving on
     */
    void calculateStressTensorInformation(const GridView& gridView)
    {
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
        }
    }

    /*!
     * \brief Calculates the Vorticity Tensor Information in each element
     *
     *  First the Vorticity Tensor matrix is filled,
     *  then the scalar product of this matrix is calculated.
     *
     * \param gridView the grid view we are solving on
     */
    void calculateVorticityTensorInformation(const GridView& gridView)
    {
        // calculate or call all secondary variables
        for (const auto& element : elements(this->fvGridGeometry().gridView()))
        {
            unsigned int elementIdx = this->fvGridGeometry().elementMapper().index(element);
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
        }
    }

    /*!
     * \brief Store the kinematic viscosity in each element
     *
     *  Using the previous solution, the kinematic viscosity for each element is stored.
     *
     * \param gridView the grid view we are solving on
     */
    void updateKinematicViscosity(const GridView& gridView, const SolutionVector& curSol)
    {
        // calculate or call all secondary variables
        for (const auto& element : elements(this->fvGridGeometry().gridView()))
        {
            unsigned int elementIdx = this->fvGridGeometry().elementMapper().index(element);
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
     * \brief Returns whether a wall function should be used at a given face
     *
     * \param element The element.
     * \param scvf The sub control volume face.
     * \param eqIdx The equation index.
     */
    bool useWallFunction(const Element& element,
                         const SubControlVolumeFace& scvf,
                         const int& eqIdx) const
    { return false; }

    /*!
     * \brief Returns an additional wall function momentum flux
     */
    template<class ElementVolumeVariables, class ElementFaceVariables>
    FacePrimaryVariables wallFunction(const Element& element,
                                      const FVElementGeometry& fvGeometry,
                                      const ElementVolumeVariables& elemVolVars,
                                      const ElementFaceVariables& elemFaceVars,
                                      const SubControlVolumeFace& scvf,
                                      const SubControlVolumeFace& localSubFace) const
    { return FacePrimaryVariables(0.0); }

    /*!
     * \brief  Returns an additional wall function flux for cell-centered quantities
     */
    template<class ElementVolumeVariables, class ElementFaceVariables>
    CellCenterPrimaryVariables wallFunction(const Element& element,
                                            const FVElementGeometry& fvGeometry,
                                            const ElementVolumeVariables& elemVolVars,
                                            const ElementFaceVariables& elemFaceVars,
                                            const SubControlVolumeFace& scvf) const
    { return CellCenterPrimaryVariables(0.0); }

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
    std::vector<DimMatrix> velocityGradients_;
    std::vector<DimVector> velocityMinimum_;
    std::vector<DimVector> velocityMaximum_;
    std::vector<Scalar> stressTensorScalarProduct_;
    std::vector<Scalar> vorticityTensorScalarProduct_;
    std::vector<unsigned int> wallNormalAxis_;
    std::vector<unsigned int> flowNormalAxis_;
    std::vector<Scalar> kinematicViscosity_;
    std::vector<Scalar> sandGrainRoughness_;
    std::vector<unsigned int> wallElements_;
    std::vector<GlobalPosition> wallPositions_;
    std::vector<unsigned int> wallProfileIdx_;
    std::vector<unsigned int> wallNormalAxisTemp_;

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
