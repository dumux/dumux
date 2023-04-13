// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//

#ifndef DUMUX_POISEUILLE_FLOW_VERTICAL_TEST_PROBLEM_HH
#define DUMUX_POISEUILLE_FLOW_VERTICAL_TEST_PROBLEM_HH

#include <algorithm>
#include <cctype>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/numeqvector.hh>

#include <dumux/freeflow/shallowwater/problem.hh>
#include <dumux/freeflow/shallowwater/boundaryfluxes.hh>

namespace Dumux {

/*!
 * \ingroup ShallowWaterTests
 * \brief A simple test for the 2D flow in a channel with viscous bottom friction (plane Poiseuille flow).
 * In comparison to the other Poiseuille flow test, here we assume a rectangular channel with zero friction walls
 * but viscous friction due to bottom shear stress.
 * The result is a parabolic velocity profile \f$ U(z) \f$ over the height:
 * \f[
 * U(z) = 3u(z/h - 0.5(z/h)^2)
 * \f]
 * bottom shear stress
 * \f[
 * \tau = -\mu 3u/h
 * \f]
 * and mean velocity
 * \f[
 * u = \frac{1}{3} \frac{h^2 \rho g S}{\mu}
 * \f]
 * where \f$ S \f$ denotes the bed slope in m/m.
 * Therefore, in difference to the other test, the heigh-averaged velocity is constant over the channel width
 * and to model the bottom surface, we need to include the effect of the bottom shear stress.
 * This problem uses the \ref ShallowWaterModel
 */
template<class TypeTag>
class PoiseuilleFlowProblem
: public ShallowWaterProblem<TypeTag>
{
    using ParentType = ShallowWaterProblem<TypeTag>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using ElementFluxVariablesCache = typename GridVariables::GridFluxVariablesCache::LocalView;
    using VolumeVariables = typename ElementVolumeVariables::VolumeVariables;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using NeumannFluxes = NumEqVector;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using FluidSystem = typename VolumeVariables::FluidSystem;

public:
    PoiseuilleFlowProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        // assume box grid with flow in x-direction driven by gravity
        channelWidth_ = this->gridGeometry().bBoxMax()[1] - this->gridGeometry().bBoxMin()[1];
        name_ = getParam<std::string>("Problem.Name");
        bedSlope_ = getParam<Scalar>("Problem.BedSlope");
        waterDepthBoundary_ = getParam<Scalar>("Problem.WaterDepth");
        dynamicViscosity_ = FluidSystem::viscosity(293.15, 1e5);

        exactWaterDepth_.resize(gridGeometry->numDofs(), 0.0);
        exactVelocityX_.resize(gridGeometry->numDofs(), 0.0);
        exactVelocityY_.resize(gridGeometry->numDofs(), 0.0);
    }

    const std::string& name() const
    { return name_; }

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param globalPos The position for which the boundary type is set
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes bcTypes;
        bcTypes.setAllNeumann();
        return bcTypes;
    }

    /*!
     * \brief Specifies the neumann boundary
     *
     *  We need the Riemann invariants to compute the values depending of the boundary type.
     *  Since we use a weak imposition we do not have a dirichlet value. We impose fluxes
     *  based on q, h, etc. computed with the Riemann invariants
     */
    NeumannFluxes neumann(const Element& element,
                          const FVElementGeometry& fvGeometry,
                          const ElementVolumeVariables& elemVolVars,
                          const ElementFluxVariablesCache& elemFluxVarsCache,
                          const SubControlVolumeFace& scvf) const
    {
        NeumannFluxes values(0.0);

        const auto& globalPos = scvf.ipGlobal();
        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& insideVolVars = elemVolVars[insideScv];
        const auto& unitNormal = scvf.unitOuterNormal();
        const auto gravity = this->spatialParams().gravity(globalPos);
        std::array<Scalar, 3> boundaryStateVariables;

        if (globalPos[1] > this->gridGeometry().bBoxMax()[1] - eps_
            || globalPos[1] < this->gridGeometry().bBoxMin()[1] + eps_)
        {
            // full slip for tangential part and reflection (no-flow) for normal part
            const Scalar insideVelocityNormalWall =  insideVolVars.velocity(0)*unitNormal[0] + insideVolVars.velocity(1)*unitNormal[1];
            const Scalar insideVelocityTangentWall = -insideVolVars.velocity(0)*unitNormal[1] + insideVolVars.velocity(1)*unitNormal[0];
            const Scalar outsideVelocityNormalWall = -insideVelocityNormalWall;
            const Scalar outsideVelocityTangentWall = insideVelocityTangentWall;
            const Scalar outsideVelocityXWall = outsideVelocityNormalWall*unitNormal[0] - outsideVelocityTangentWall*unitNormal[1];
            const Scalar outsideVelocityYWall = outsideVelocityNormalWall*unitNormal[1] + outsideVelocityTangentWall*unitNormal[0];
            boundaryStateVariables = { insideVolVars.waterDepth(), outsideVelocityXWall, outsideVelocityYWall };
        }
        else
        {
            // for inlet and outlet, the same water depth is prescribed
            boundaryStateVariables = ShallowWater::fixedWaterDepthBoundary(
                waterDepthBoundary_,
                insideVolVars.waterDepth(), insideVolVars.velocity(0), insideVolVars.velocity(1),
                gravity, unitNormal
            );
        }

        auto riemannFlux = ShallowWater::riemannProblem(
            insideVolVars.waterDepth(), boundaryStateVariables[0],
            insideVolVars.velocity(0), boundaryStateVariables[1],
            insideVolVars.velocity(1), boundaryStateVariables[2],
            insideVolVars.bedSurface(), insideVolVars.bedSurface(),
            gravity, unitNormal
        );

        values[Indices::massBalanceIdx] = riemannFlux[0];
        values[Indices::velocityXIdx] = riemannFlux[1];
        values[Indices::velocityYIdx] = riemannFlux[2];

        // zero viscous part of the flux since we assume full slip on the walls

        return values;
    }

    // bottom friction source
    NumEqVector source(const Element &element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolume& scv) const
    {
        NumEqVector bottomFrictionSource(0.0);

        const auto& volVars = elemVolVars[scv];
        Dune::FieldVector<Scalar, 2> bottomShearStress =
            this->spatialParams().frictionLaw(element, scv).bottomShearStress(volVars);

        bottomFrictionSource[0] = 0.0;
        bottomFrictionSource[1] = -bottomShearStress[0] / volVars.density();
        bottomFrictionSource[2] = -bottomShearStress[1] / volVars.density();

        return bottomFrictionSource;
    }

    // Set initial sol to the exact solution
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    {
       return exactSol_(globalPos);
    };

    //! Update the analytical solution
    void updateAnalyticalSolution()
    {
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            const auto eIdx = this->gridGeometry().elementMapper().index(element);
            const auto& globalPos = element.geometry().center();
            const auto sol = exactSol_(globalPos);
            exactWaterDepth_[eIdx] = sol[Indices::waterdepthIdx];
            exactVelocityX_[eIdx] = sol[Indices::velocityXIdx];
            exactVelocityY_[eIdx] = sol[Indices::velocityYIdx];
        }
    }

    //! Get the analytical water depth
    const std::vector<Scalar>& getExactWaterDepth() const
    { return exactWaterDepth_; }

    //! Get the analytical U-velocity
    const std::vector<Scalar>& getExactVelocityX() const
    { return exactVelocityX_; }

    //! Get the analytical V-velocity
    const std::vector<Scalar>& getExactVelocityY() const
    { return exactVelocityY_; }

private:
    PrimaryVariables exactSol_(const GlobalPosition& globalPos) const
    {
        PrimaryVariables values(0.0);
        const auto gravity = this->spatialParams().gravity(globalPos);
        const auto hSquared = waterDepthBoundary_*waterDepthBoundary_;
        const auto pressureGradient = 1000.0*gravity*bedSlope_;
        values[0] = waterDepthBoundary_;
        values[1] = pressureGradient*hSquared/(3.0*dynamicViscosity_);
        values[2] = 0.0;
        return values;
    }

    std::vector<Scalar> exactWaterDepth_;
    std::vector<Scalar> exactVelocityX_;
    std::vector<Scalar> exactVelocityY_;

    Scalar channelWidth_;
    Scalar waterDepthBoundary_; // water level
    Scalar bedSlope_; // bed slope (positive downwards)
    Scalar dynamicViscosity_;

    static constexpr Scalar eps_ = 1.0e-6;
    std::string name_;
};

} //end namespace Dumux

#endif
