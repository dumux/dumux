// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup RANSModel
 * \copydoc Dumux::RANSProblem
 */
#ifndef DUMUX_RANS_PROBLEM_HH
#define DUMUX_RANS_PROBLEM_HH

#include <algorithm>

#include <dune/common/fmatrix.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/staggeredfvproblem.hh>
#include <dumux/discretization/localview.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/walldistance.hh>
#include <dumux/discretization/staggered/elementsolution.hh>
#include <dumux/freeflow/navierstokes/staggered/problem.hh>
#include "model.hh"

namespace Dumux {

//! forward declare
template<class TypeTag, TurbulenceModel turbulenceModel>
class RANSProblemImpl;

//! the turbulence-model-specfic RANS problem
template<class TypeTag>
using RANSProblem = RANSProblemImpl<TypeTag, GetPropType<TypeTag, Properties::ModelTraits>::turbulenceModel()>;

/*!
 * \ingroup RANSModel
 * \brief Reynolds-Averaged Navier-Stokes problem base class.
 *
 * This implements some base functionality for RANS models.
 * Especially vectors containing all wall-relevant properties, which are accessed
 * by the volumevariables.
 */
template<class TypeTag>
class RANSProblemBase : public NavierStokesStaggeredProblem<TypeTag>
{
    using ParentType = NavierStokesStaggeredProblem<TypeTag>;
    using Implementation = GetPropType<TypeTag, Properties::Problem>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using PrimaryVariables = typename VolumeVariables::PrimaryVariables;
    using CellCenterPrimaryVariables = GetPropType<TypeTag, Properties::CellCenterPrimaryVariables>;
    using FacePrimaryVariables = GetPropType<TypeTag, Properties::FacePrimaryVariables>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    using GlobalPosition = typename SubControlVolumeFace::GlobalPosition;

    static constexpr auto dim = GridView::dimension;
    static constexpr int numCorners = SubControlVolumeFace::numCornersPerFace;
    using DimVector = GlobalPosition;
    using DimMatrix = Dune::FieldMatrix<Scalar, dim, dim>;

    struct WallElementInformation
    {
        // store the element indices for all elements with an intersection on the wall
        unsigned int wallElementIdx;
        // for each wall element, store the faces normal axis
        unsigned int wallFaceNormalAxis;
        // for each wall element, store the location of the face center and each corner.
        GlobalPosition wallFaceCenter;
        std::array<GlobalPosition, numCorners> wallFaceCorners;
    };

public:

    /*!
     * \brief The constructor
     * \param gridGeometry The finite volume grid geometry
     * \param paramGroup The parameter group in which to look for runtime parameters first (default is "")
     */
    RANSProblemBase(std::shared_ptr<const GridGeometry> gridGeometry, const std::string& paramGroup = "")
    : ParentType(gridGeometry, paramGroup)
    {
        if ( !(hasParamInGroup(this->paramGroup(), "RANS.IsFlatWallBounded")))
        {
            std::cout << "The parameter \"Rans.IsFlatWallBounded\" is not specified. \n"
                    << " -- Based on the grid and the boundary conditions specified by the user,"
                    << " this parameter is set to be "<< std::boolalpha << isFlatWallBounded() << "\n";
        }

        // update size and initial values of the global vectors
        wallDistance_.resize(this->gridGeometry().elementMapper().size(), std::numeric_limits<Scalar>::max());
        neighborIdx_.resize(this->gridGeometry().elementMapper().size());
        velocity_.resize(this->gridGeometry().elementMapper().size(), DimVector(0.0));
        velocityGradients_.resize(this->gridGeometry().elementMapper().size(), DimMatrix(0.0));
        stressTensorScalarProduct_.resize(this->gridGeometry().elementMapper().size(), 0.0);
        vorticityTensorScalarProduct_.resize(this->gridGeometry().elementMapper().size(), 0.0);
        flowDirectionAxis_.resize(this->gridGeometry().elementMapper().size(), fixedFlowDirectionAxis_);
        storedViscosity_.resize(this->gridGeometry().elementMapper().size(), 0.0);
        storedDensity_.resize(this->gridGeometry().elementMapper().size(), 0.0);
    }

    /*!
     * \brief Update the static (solution independent) relations to the walls and neighbors
     */
    void updateStaticWallProperties()
    {
        std::cout << "Update static wall properties. ";
        calledUpdateStaticWallProperties = true;

        checkForWalls_();
        findWallDistances_();
        findNeighborIndices_();
    }

    /*!
     * \brief Update the dynamic (solution dependent) turbulence parameters
     *
     * \param curSol The solution vector.
     */
    template<class SolutionVector>
    void updateDynamicWallProperties(const SolutionVector& curSol)
    {
        std::cout << "Update dynamic wall properties." << std::endl;
        if (!calledUpdateStaticWallProperties)
            DUNE_THROW(Dune::InvalidStateException,
                       "You have to call updateStaticWallProperties() once before you call updateDynamicWallProperties().");

        calculateCCVelocities_(curSol);
        calculateCCVelocityGradients_();
        calculateMaxMinVelocities_();
        calculateStressTensor_();
        calculateVorticityTensor_();
        storeViscosities_(curSol);
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
                                      const SubControlVolumeFace& lateralBoundaryFace) const
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
     */
    bool isFlatWallBounded() const
    {
        static const bool hasAlignedWalls = hasAlignedWalls_();
        return hasAlignedWalls;
    }

    /*!
     * \brief Returns the Karman constant
     */
    const Scalar karmanConstant() const
    { return 0.41; }

    //! \brief Returns the \f$ \beta_{\omega} \f$ constant
    const Scalar betaOmega() const
    { return 0.0708; }

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

    int wallNormalAxis(const int elementIdx) const
    {
        if (!isFlatWallBounded())
            DUNE_THROW(Dune::NotImplemented, "\n Due to grid/geometric concerns, models requiring a wallNormalAxis "
                                          << "can only be used for flat wall bounded flows. "
                                          << "\n If your geometry is a flat channel, "
                                          << "please set the runtime parameter RANS.IsFlatWallBounded to true. \n");
        return wallNormalAxis_[elementIdx];
    }

    int flowDirectionAxis(const int elementIdx) const
    {
        if (!isFlatWallBounded())
            DUNE_THROW(Dune::NotImplemented, "\n Due to grid/geometric concerns, models requiring a flowDirectionAxis "
                                          << "can only be used for flat wall bounded flows. "
                                          << "\n If your geometry is a flat channel, "
                                          << "please set the runtime parameter RANS.IsFlatWallBounded to true. \n");
        return flowDirectionAxis_[elementIdx];
    }

    unsigned int wallElementIndex(const int elementIdx) const
    {
        if (!isFlatWallBounded())
            DUNE_THROW(Dune::NotImplemented, "\n Due to grid/geometric concerns, models requiring a wallElementIndex "
                                          << "can only be used for flat wall bounded flows. "
                                          << "\n If your geometry is a flat channel, "
                                          << "please set the runtime parameter RANS.IsFlatWallBounded to true. \n");
        return wallElementIdx_[elementIdx];

    }

    Scalar wallDistance(const int elementIdx) const
    { return wallDistance_[elementIdx]; }

    GlobalPosition cellCenter(const int elementIdx) const
    {
        const auto& element = this->gridGeometry().element(elementIdx);
        return element.geometry().center();
    }

    unsigned int neighborIndex(const int elementIdx, const int axisIdx, const int sideIdx) const
    { return neighborIdx_[elementIdx][axisIdx][sideIdx];}

    DimVector ccVelocityVector(const int elementIdx) const
    { return velocity_[elementIdx]; }

    Scalar ccVelocity(const int elementIdx, const int axisIdx) const
    { return velocity_[elementIdx][axisIdx]; }

    DimVector velocityMaximum(const int elementIdx) const
    { return velocityMaximum_[elementIdx]; }

    DimVector velocityMinimum(const int elementIdx) const
    { return velocityMinimum_[elementIdx]; }

    DimMatrix velocityGradientTensor(const int elementIdx) const
    { return velocityGradients_[elementIdx]; }

    Scalar velocityGradient(const int elementIdx, const int i, const int j) const
    { return velocityGradients_[elementIdx][i][j]; }

    Scalar stressTensorScalarProduct(const int elementIdx) const
    { return stressTensorScalarProduct_[elementIdx]; }

    Scalar vorticityTensorScalarProduct(const int elementIdx) const
    { return vorticityTensorScalarProduct_[elementIdx]; }

    Scalar storedViscosity(const int elementIdx) const
    { return storedViscosity_[elementIdx]; }

    Scalar storedDensity(const int elementIdx) const
    { return storedDensity_[elementIdx]; }

    Scalar kinematicViscosity(const int elementIdx) const
    { return storedViscosity(elementIdx) / storedDensity(elementIdx); }

    bool calledUpdateStaticWallProperties = false;

private:

    bool hasAlignedWalls_() const
    {
        if ( hasParamInGroup(this->paramGroup(), "RANS.IsFlatWallBounded"))
        {
            static const bool isFlatWallBounded = getParamFromGroup<bool>(this->paramGroup(), "RANS.IsFlatWallBounded");
            return isFlatWallBounded;
        }

        std::vector<int> wallFaceAxis;
        wallFaceAxis.reserve(this->gridGeometry().numBoundaryScvf());

        const auto gridView = this->gridGeometry().gridView();
        auto fvGeometry = localView(this->gridGeometry());
        for (const auto& element : elements(gridView))
        {
            fvGeometry.bindElement(element);
            for (const auto& scvf : scvfs(fvGeometry))
                if (!scvf.boundary() && asImp_().boundaryTypes(element, scvf).hasWall())  // only search for walls at a global boundary
                    wallFaceAxis.push_back(scvf.directionIndex());
        }

        // Returns if all wall directions are the same
        return std::all_of(wallFaceAxis.begin(), wallFaceAxis.end(), [firstDir=wallFaceAxis[0]](auto dir){ return (dir == firstDir);} ) ;
    }

    void checkForWalls_()
    {
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            auto fvGeometry = localView(this->gridGeometry());
            fvGeometry.bindElement(element);
            for (auto&& scvf : scvfs(fvGeometry))
                if (asImp_().boundaryTypes(element, scvf).hasWall())
                    return;
        }
        // If reached, no walls were found using the boundary types has wall function.
        DUNE_THROW(Dune::InvalidStateException, "No walls are are specified with the setWall() function");
    }

    /*!
     * \brief Use the boundary search algorithm to find the shortest distance to a wall for each element
     *
     *  Also store the wall element's index, and its direction in the case of flat wall bounded problems
     */
    void findWallDistances_()
    {
        WallDistance wallInformation(this->gridGeometry(), WallDistance<GridGeometry>::atElementCenters,
            [this] (const FVElementGeometry& fvGeometry, const SubControlVolumeFace& scvf)
            { return asImp_().boundaryTypes(fvGeometry.element(), scvf).hasWall(); });
        wallDistance_ = wallInformation.wallDistance();
        storeWallElementAndDirectionIndex_(wallInformation.wallData());
    }

    template <class WallData>
    void storeWallElementAndDirectionIndex_(const WallData& wallData)
    {
        // The wall Direction Index is used for flat quadrilateral channel problems only
        if (!(GridGeometry::discMethod == DiscretizationMethods::staggered))
            DUNE_THROW(Dune::NotImplemented, "The wall direction Index can only be calculated for quadrilateral structured grids");

        // If isFlatWallBounded, the corresponding wall element is stored for each element
        if (isFlatWallBounded())
        {
            wallNormalAxis_.resize(wallData.size());
            wallElementIdx_.resize(wallData.size());

            for (const auto& element : elements(this->gridGeometry().gridView()))
            {
                unsigned int elementIdx = this->gridGeometry().elementMapper().index(element);
                wallElementIdx_[elementIdx] = wallData[elementIdx].eIdx;
                if ( ! (hasParam("RANS.WallNormalAxis")) )
                {
                    GlobalPosition wallOuterNormal = wallData[elementIdx].scvfOuterNormal;
                    if constexpr (dim == 2) // 2D
                        wallNormalAxis_[elementIdx] = (wallOuterNormal[0] == 1) ? 0 : 1;
                    else // 3D
                        wallNormalAxis_[elementIdx] = (wallOuterNormal[0] == 1) ? 0 : ((wallOuterNormal[1] == 1) ? 1 : 2);
                }
                else
                    wallNormalAxis_[elementIdx] = fixedWallNormalAxis_;
            }
        }
    }

    /*!
     * \brief Store all direct neighbor indices for each element
     */
    void findNeighborIndices_()
    {
        // search for neighbor Idxs
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            unsigned int elementIdx = this->gridGeometry().elementMapper().index(element);
            for (unsigned int axisIdx = 0; axisIdx < dim; ++axisIdx)
            {
                neighborIdx_[elementIdx][axisIdx][0] = elementIdx;
                neighborIdx_[elementIdx][axisIdx][1] = elementIdx;
            }

            for (const auto& intersection : intersections(this->gridGeometry().gridView(), element))
            {
                if (intersection.boundary())
                    continue;

                unsigned int neighborIdx = this->gridGeometry().elementMapper().index(intersection.outside());
                for (unsigned int axisIdx = 0; axisIdx < dim; ++axisIdx)
                {
                    if (abs(cellCenter(elementIdx)[axisIdx] - cellCenter(neighborIdx)[axisIdx]) > 1e-8)
                    {
                        if (cellCenter(elementIdx)[axisIdx] > cellCenter(neighborIdx)[axisIdx])
                            neighborIdx_[elementIdx][axisIdx][0] = neighborIdx;

                        if (cellCenter(elementIdx)[axisIdx] < cellCenter(neighborIdx)[axisIdx])
                            neighborIdx_[elementIdx][axisIdx][1] = neighborIdx;
                    }
                }
            }
        }
    }

    template<class SolutionVector>
    void calculateCCVelocities_(const SolutionVector& curSol)
    {
        auto fvGeometry = localView(this->gridGeometry());
        // calculate cell-center-averaged velocities
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            fvGeometry.bindElement(element);
            unsigned int elementIdx = this->gridGeometry().elementMapper().index(element);

            // calculate velocities
            DimVector velocityTemp(0.0);
            for (auto&& scvf : scvfs(fvGeometry))
            {
                const int dofIdxFace = scvf.dofIndex();
                const auto numericalSolutionFace = curSol[GridGeometry::faceIdx()][dofIdxFace][Indices::velocity(scvf.directionIndex())];
                velocityTemp[scvf.directionIndex()] += numericalSolutionFace;
            }
            for (unsigned int axisIdx = 0; axisIdx < dim; ++axisIdx)
                velocity_[elementIdx][axisIdx] = velocityTemp[axisIdx] * 0.5; // faces are equidistant to cell center
        }
    }


    void calculateCCVelocityGradients_()
    {
        using std::abs;

        // calculate cell-center-averaged velocity gradients, maximum, and minimum values
        auto fvGeometry = localView(this->gridGeometry());
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            const unsigned int elementIdx = this->gridGeometry().elementMapper().index(element);

            for (unsigned int j = 0; j < dim; ++j)
            {
                for (unsigned int i = 0; i < dim; ++i)
                {
                    const unsigned int neighborIndex0 = neighborIndex(elementIdx, j, 0);
                    const unsigned int neighborIndex1 = neighborIndex(elementIdx, j, 1);

                    velocityGradients_[elementIdx][i][j]
                        = (ccVelocity(neighborIndex1, i) - ccVelocity(neighborIndex0, i))
                        / (cellCenter(neighborIndex1)[j] - cellCenter(neighborIndex0)[j]);

                    if (abs(cellCenter(neighborIndex1)[j] - cellCenter(neighborIndex0)[j]) < 1e-8)
                        velocityGradients_[elementIdx][i][j] = 0.0;
                }
            }

            fvGeometry.bindElement(element);
            for (auto&& scvf : scvfs(fvGeometry))
            {
                // adapt calculations for Dirichlet condition
                unsigned int axisIdx = scvf.directionIndex();
                if (scvf.boundary())
                {
                    for (unsigned int velIdx = 0; velIdx < dim; ++velIdx)
                    {
                        if (!asImp_().boundaryTypes(element, scvf).isDirichlet(Indices::velocity(velIdx)))
                            continue;

                        Scalar dirichletVelocity = asImp_().dirichlet(element, scvf)[Indices::velocity(velIdx)];

                        unsigned int neighborIdx = neighborIndex(elementIdx, axisIdx, 0);
                        if (scvf.center()[axisIdx] < cellCenter(elementIdx)[axisIdx])
                            neighborIdx = neighborIndex(elementIdx, axisIdx, 1);

                        velocityGradients_[elementIdx][velIdx][axisIdx]
                            = (ccVelocity(neighborIdx, velIdx) - dirichletVelocity)
                              / (cellCenter(neighborIdx)[axisIdx] - scvf.center()[axisIdx]);
                    }
                }

                // Calculate the BJS-velocity by accounting for all sub faces.
                std::vector<int> bjsNumFaces(dim, 0);
                std::vector<unsigned int> bjsNeighbor(dim, 0);
                DimVector bjsVelocityAverage(0.0);
                DimVector normalNormCoordinate(0.0);
                unsigned int velCompIdx = Indices::velocity(scvf.directionIndex());
                const int numSubFaces = scvf.pairData().size();
                for(int localSubFaceIdx = 0; localSubFaceIdx < numSubFaces; ++localSubFaceIdx)
                {
                    const auto& lateralFace = fvGeometry.scvf(scvf.insideScvIdx(), scvf.pairData()[localSubFaceIdx].localLateralFaceIdx);

                    // adapt calculations for Beavers-Joseph-Saffman condition
                    unsigned int lateralAxisIdx = lateralFace.directionIndex();
                    if (lateralFace.boundary() && (asImp_().boundaryTypes(element, lateralFace).isBeaversJoseph(Indices::velocity(velCompIdx))))
                    {
                        unsigned int neighborIdx = neighborIndex(elementIdx, lateralAxisIdx, 0);
                        if (lateralFace.center()[lateralAxisIdx] < cellCenter(elementIdx)[lateralAxisIdx])
                            neighborIdx = neighborIndex(elementIdx, lateralAxisIdx, 1);

                        const SubControlVolume& scv = fvGeometry.scv(scvf.insideScvIdx());
                        bjsVelocityAverage[lateralAxisIdx] += ParentType::beaversJosephVelocity(element, scv, scvf, lateralFace, ccVelocity(elementIdx, velCompIdx), 0.0);
                        if (bjsNumFaces[lateralAxisIdx] > 0 && neighborIdx != bjsNeighbor[lateralAxisIdx])
                            DUNE_THROW(Dune::InvalidStateException, "Two different neighborIdx should not occur");
                        bjsNeighbor[lateralAxisIdx] = neighborIdx;
                        normalNormCoordinate[lateralAxisIdx] = lateralFace.center()[lateralAxisIdx];
                        bjsNumFaces[lateralAxisIdx]++;
                    }
                }
                for (unsigned axisIdx = 0; axisIdx < dim; ++axisIdx)
                {
                    if (bjsNumFaces[axisIdx] == 0)
                        continue;

                    unsigned int neighborIdx = bjsNeighbor[axisIdx];
                    bjsVelocityAverage[axisIdx] /= bjsNumFaces[axisIdx];

                    velocityGradients_[elementIdx][velCompIdx][axisIdx]
                        = (ccVelocity(neighborIdx, velCompIdx) - bjsVelocityAverage[axisIdx])
                        / (cellCenter(neighborIdx)[axisIdx] - normalNormCoordinate[axisIdx]);
                }
            }
        }
    }

    void calculateMaxMinVelocities_()
    {
        using std::abs;
        if (isFlatWallBounded())
        {
            // If the parameter isFlatWallBounded is set to true,
            // the maximum/minimum velocities are calculated along a profile perpendicular to the corresponding wall face.

            // re-initialize min and max values
            velocityMaximum_.assign(this->gridGeometry().elementMapper().size(), DimVector(1e-16));
            velocityMinimum_.assign(this->gridGeometry().elementMapper().size(), DimVector(std::numeric_limits<Scalar>::max()));

            // For each profile perpendicular to the channel wall, find the max and minimum velocities
            for (const auto& element : elements(this->gridGeometry().gridView()))
            {
                const unsigned int elementIdx = this->gridGeometry().elementMapper().index(element);
                Scalar maxVelocity = 0.0;
                const unsigned int wallElementIdx = wallElementIndex(elementIdx);

                for (unsigned int axisIdx = 0; axisIdx < dim; ++axisIdx)
                {
                    if (abs(ccVelocity(elementIdx, axisIdx)) > abs(velocityMaximum_[wallElementIdx][axisIdx]))
                        velocityMaximum_[wallElementIdx][axisIdx] = ccVelocity(elementIdx, axisIdx);

                    if (abs(ccVelocity(elementIdx, axisIdx)) < abs(velocityMinimum_[wallElementIdx][axisIdx]))
                        velocityMinimum_[wallElementIdx][axisIdx] = ccVelocity(elementIdx, axisIdx);

                    // Set the flow direction axis as the direction of the max velocity
                    if ((hasParam("RANS.FlowDirectionAxis") != 1) && (maxVelocity) < abs(ccVelocity(elementIdx, axisIdx)))
                    {
                        maxVelocity = abs(ccVelocity(elementIdx, axisIdx));
                        flowDirectionAxis_[elementIdx] = axisIdx;
                    }
                }
            }
        }
        else
        {
            // If the parameter isFlatWallBounded is set to false, or not set,
            // the maximum/minimum velocities are calculated as a global max/min throughout the domain.

            DimVector maxVelocity(0.0);
            DimVector minVelocity(std::numeric_limits<Scalar>::max());
            // Find the max and minimum velocities in the full domain
            for (const auto& element : elements(this->gridGeometry().gridView()))
            {
                const unsigned int elementIdx = this->gridGeometry().elementMapper().index(element);

                for (unsigned int axisIdx = 0; axisIdx < dim; ++axisIdx)
                {
                    if (abs(ccVelocity(elementIdx, axisIdx)) > abs(maxVelocity[axisIdx]))
                        maxVelocity[axisIdx] = ccVelocity(elementIdx, axisIdx);

                    if (abs(ccVelocity(elementIdx, axisIdx)) < abs(minVelocity[axisIdx]))
                        minVelocity[axisIdx] = ccVelocity(elementIdx, axisIdx);
                }
            }
            velocityMaximum_.assign(this->gridGeometry().elementMapper().size(), maxVelocity);
            velocityMinimum_.assign(this->gridGeometry().elementMapper().size(), minVelocity);
        }
    }

    void calculateStressTensor_()
    {
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            unsigned int elementIdx = this->gridGeometry().elementMapper().index(element);
            Dune::FieldMatrix<Scalar, GridView::dimension, GridView::dimension> stressTensor(0.0);
            for (unsigned int j = 0; j < dim; ++j)
            {
                for (unsigned int i = 0; i < dim; ++i)
                {
                    stressTensor[j][i] = 0.5 * velocityGradient(elementIdx, j, i)
                                                 + 0.5 * velocityGradient(elementIdx, i, j);
              }
            }
            stressTensorScalarProduct_[elementIdx] = 0.0;
            for (unsigned int j = 0; j < dim; ++j)
            {
                for (unsigned int i = 0; i < dim; ++i)
                {
                    stressTensorScalarProduct_[elementIdx] += stressTensor[j][i] * stressTensor[j][i];
                }
            }
        }
    }

    void calculateVorticityTensor_()
    {
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            unsigned int elementIdx = this->gridGeometry().elementMapper().index(element);
            Dune::FieldMatrix<Scalar, GridView::dimension, GridView::dimension> vorticityTensor(0.0);
            for (unsigned int j = 0; j < dim; ++j)
            {
                for (unsigned int i = 0; i < dim; ++i)
                {
                    vorticityTensor[j][i] = 0.5 * velocityGradient(elementIdx, j, i)
                                                    - 0.5 * velocityGradient(elementIdx, i, j);
              }
            }
            vorticityTensorScalarProduct_[elementIdx] = 0.0;
            for (unsigned int j = 0; j < dim; ++j)
            {
                for (unsigned int i = 0; i < dim; ++i)
                {
                    vorticityTensorScalarProduct_[elementIdx] += vorticityTensor[j][i] * vorticityTensor[j][i];
                }
            }
        }
    }

    template<class SolutionVector>
    void storeViscosities_(const SolutionVector& curSol)
    {
        // calculate or call all secondary variables
        auto fvGeometry = localView(this->gridGeometry());
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            unsigned int elementIdx = this->gridGeometry().elementMapper().index(element);
            fvGeometry.bindElement(element);
            for (auto&& scv : scvs(fvGeometry))
            {
                const int dofIdx = scv.dofIndex();
                // construct a privars object from the cell center solution vector
                const auto& cellCenterPriVars = curSol[GridGeometry::cellCenterIdx()][dofIdx];
                PrimaryVariables priVars = makePriVarsFromCellCenterPriVars<PrimaryVariables>(cellCenterPriVars);
                auto elemSol = elementSolution<typename GridGeometry::LocalView>(std::move(priVars));

                VolumeVariables volVars;
                volVars.update(elemSol, asImp_(), element, scv);
                storedDensity_[elementIdx] = volVars.density();
                storedViscosity_[elementIdx] = volVars.viscosity();
            }
        }
    }

    const int fixedFlowDirectionAxis_ = getParam<int>("RANS.FlowDirectionAxis", 0);
    const int fixedWallNormalAxis_ = getParam<int>("RANS.WallNormalAxis", 1);

    std::vector<unsigned int> wallNormalAxis_;
    std::vector<unsigned int> flowDirectionAxis_;
    std::vector<Scalar> wallDistance_;
    std::vector<unsigned int> wallElementIdx_;
    std::vector<std::array<std::array<unsigned int, 2>, dim>> neighborIdx_;

    std::vector<DimVector> velocity_;
    std::vector<DimVector> velocityMaximum_;
    std::vector<DimVector> velocityMinimum_;
    std::vector<DimMatrix> velocityGradients_;

    std::vector<Scalar> stressTensorScalarProduct_;
    std::vector<Scalar> vorticityTensorScalarProduct_;

    std::vector<Scalar> storedDensity_;
    std::vector<Scalar> storedViscosity_;

    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }
};

} // end namespace Dumux

#endif
