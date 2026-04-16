// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_TEST_FREEFLOW_TWOP_PORES_EVAPORATION_PROBLEM_HH
#define DUMUX_TEST_FREEFLOW_TWOP_PORES_EVAPORATION_PROBLEM_HH

#include <dune/common/float_cmp.hh>
#include <algorithm>
#include <cmath>
#include <vector>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/geometry/diameter.hh>

namespace Dumux {

/*!
 * \brief Test problem for two-phase flow in porous media with evaporation
 *
 * Water fills the pores (below pillars) and evaporates from the free surface.
 * Hele-Shaw drag acts around pillar surfaces.
 */
template <class TypeTag, class BaseProblem>
class TwoPhaseEvaporationProblem : public BaseProblem
{
    using ParentType = BaseProblem;
    using BoundaryTypes = typename ParentType::BoundaryTypes;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using BoundaryFluxes = typename ParentType::BoundaryFluxes;
    using DirichletValues = typename ParentType::DirichletValues;
    using InitialValues = typename ParentType::InitialValues;
    using Sources = typename ParentType::Sources;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVariables>::GridVolumeVariables::LocalView;

    static constexpr int dim = GridGeometry::GridView::dimension;
    static constexpr auto dimWorld = GridGeometry::GridView::dimensionworld;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using Vertex = typename GridGeometry::GridView::template Codim<dim>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using Extrusion = Extrusion_t<GridGeometry>;

    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

public:
    TwoPhaseEvaporationProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                                std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry, couplingManager, ParentType::isMomentumProblem() ? "Momentum" : "Mass")
    {

        surfaceTension_ = getParam<Scalar>("Problem.SurfaceTension", 24.5);
        interfaceThickness_ = getParam<Scalar>("Problem.InterfaceThickness", 0.08);

        // Mobility M controls diffusion rate in CH equation: flux = -M∇μ
        // Typical scaling: M ~ ε² for interface-controlled dynamics
        mobility_ = getParam<Scalar>("Problem.Mobility", 1.0);

        // Surface tension energy γ appears in: μ = -γε∇²φ + (γ/ε)f'(φ)
        // For f(φ) = (1/4)(φ²-1)², the surface energy integral gives σ_eff = γ·(2√2/3),
        // so to recover the physical surface tension σ, we need γ = 3σ/(2√2).
        scaledSurfaceTension_ = surfaceTension_ * 3.0 / (2.0 * std::sqrt(2.0)) * getParam<Scalar>("Problem.SurfaceTensionFactor", 1.0);

        // Energy scale (γ/ε) for the double-well potential derivative: f'(φ) = (γ/ε)φ(φ²-1)
        energyScale_ = scaledSurfaceTension_ / interfaceThickness_ * getParam<Scalar>("Problem.EnergyScaleFactor", 1.0);

        // Read parameters
        rho1_ = getParam<Scalar>("Problem.Density1", 998.0);   // water
        rho2_ = getParam<Scalar>("Problem.Density2", 1.225);   // air
        eta1_ = getParam<Scalar>("Problem.Viscosity1", 0.001); // water
        eta2_ = getParam<Scalar>("Problem.Viscosity2", 1.8e-5); // air

        evaporationRate_ = getParam<Scalar>("Problem.EvaporationRate", 1.0e-5);
        inletVelocity_   = getParam<Scalar>("Problem.InletVelocity", 0.0);
        enableKortewegForce_ = getParam<bool>("Problem.EnableKortewegForce", true);

        // Contact angle on pillar surfaces (in radians). 90° gives the natural BC (∂φ/∂n = 0).
        contactAngle_ = getParam<Scalar>("Problem.ContactAngleDegrees", 90.0) * M_PI / 180.0;

        // Compute pillar centers using the same geometry as channel_with_pillars.geo
        pillarRadius_ = getParam<Scalar>("Problem.PillarRadius", 0.025);
        {
            const int nx = getParam<int>("Problem.PillarNx", 6);
            const int ny = getParam<int>("Problem.PillarNy", 4);
            const Scalar channelL = 1.0;
            const Scalar boxW = 0.5 * channelL;
            const Scalar boxH = 0.4;
            const Scalar margin = 0.1; // distance from top pillar row to gas channel

            const Scalar xBoxLeft  = 0.5*(channelL - boxW);
            const Scalar xBoxRight = xBoxLeft + boxW;
            const Scalar yBoxTop   = 0.0;
            const Scalar yBoxBot   = -boxH;

            const Scalar xInnerLeft  = xBoxLeft  + pillarRadius_;
            const Scalar xInnerRight = xBoxRight - pillarRadius_;
            const Scalar yInnerTop   = yBoxTop - pillarRadius_ - margin;
            const Scalar yInnerBot   = yBoxBot + pillarRadius_;

            const Scalar dx = 3.0 * pillarRadius_;
            const Scalar dy = std::sqrt(3.0) / 2.0 * dx;

            const Scalar widthNeeded  = (nx - 1)*dx + 0.5*dx;
            const Scalar heightNeeded = (ny - 1)*dy;
            const Scalar xOffset = (xInnerRight - xInnerLeft - widthNeeded) / 2.0;
            const Scalar yOffset = (yInnerTop   - yInnerBot  - heightNeeded) / 2.0;

            pillarCenters_.resize(nx * ny);
            for (int j = 0; j < ny; ++j)
                for (int i = 0; i < nx; ++i)
                    pillarCenters_[j*nx + i] = {
                        xInnerLeft + xOffset + i*dx + (j%2)*0.5*dx,
                        yInnerBot  + yOffset + j*dy
                    };
        }
    }

    Scalar mobility() const
    {
        return mobility_;
    }

    Scalar surfaceTension() const
    {
        // Returns γε, the coefficient in front of the Laplacian in the chem. potential
        return scaledSurfaceTension_ * interfaceThickness_;
    }

    /*!
     * \brief Specifies Dirichlet boundary values.
     */
    DirichletValues dirichletAtPos(const GlobalPosition& globalPos) const
    {
        DirichletValues values(0.0);
        if constexpr (ParentType::isMomentumProblem())
        {
            if (isInlet_(globalPos))
                values[0] = parabolicProfile_(globalPos[dimWorld-1]);
        }
        // mass problem: no Dirichlet BCs used (all Neumann)

        return values;
    }

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param element The finite element
     * \param scv The sub control volume
     */
    BoundaryTypes boundaryTypes(const Element& element,
                                const SubControlVolume& scv) const
    {
        BoundaryTypes values;

        if constexpr (ParentType::isMomentumProblem())
        {
            // Inlet: Dirichlet (prescribed velocity).
            // Outlet: Neumann (free outflow).
            // Everything else: Dirichlet no-slip (all solid walls including pillars).
            if (isOutlet_(scv.dofPosition()))
                values.setAllNeumann();
            else
                values.setAllDirichlet();
        }
        else
        {
            values.setAllNeumann();
        }

        return values;
    }

    /*!
     * \brief Return the sources within the domain.
     *
     * \param globalPos The global position
     */
    Sources source(const Element &element,
                   const FVElementGeometry& fvGeometry,
                   const ElementVolumeVariables& elemVolVars,
                   const SubControlVolume &scv) const

    {
        Sources source(0.0);
        if constexpr (!ParentType::isMomentumProblem())
        {
            const auto& volVars = elemVolVars[scv];
            const auto phi = volVars.phaseField();

            // Double-well potential derivative: f'(φ) = (γ/ε)·φ(φ²-1)
            // f(φ) = (γ/ε)·(1/4)(φ²-1)² → f'(φ) = (γ/ε)·(1/2)·2φ(φ²-1)
            source[Indices::chemicalPotentialEqIdx] = volVars.chemicalPotential() - energyScale_*phi*(phi*phi - 1.0);

            // Evaporation drives φ → -1 (gas) at the interface.
            // The density change ∂ρ/∂t = (dρ/dφ)∂φ/∂t is handled by the ∂ρ/∂t
            // storage term in the continuity equation — no explicit mass source needed there.
            // Adding a mass sink in conti would create convergent inflow from boundaries.
            source[Indices::phaseFieldEqIdx] -= evaporationRate_/rho1_ * 3.0/2.0*(1.0 - phi*phi)/interfaceThickness_;
        }
        else
        {
            if (enableKortewegForce_)
            {
                const auto ip = ipData(fvGeometry, scv);
                const auto grads = this->couplingManager().gradients(element, fvGeometry, ip);
                const auto mu = this->couplingManager().values(element, fvGeometry, ip)[CouplingManager::chemicalPotentialIdx];
                source = mu * grads[CouplingManager::phaseFieldIdx];
            }
        }

        // if constexpr (ParentType::isMomentumProblem())
        // {
        //     if (scv.dofPosition()[1] < 0.0)
        //     {
        //         const auto v = elemVolVars[scv].velocity();
        //         source -= dampingForceMagnitude_ * v.two_norm() * v;
        //     }

        //     source -= rho2_ * this->gravity();
        // }

        return source;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Neumann control volume.
     *
     * \param element The element for which the Neumann boundary condition is set
     * \param fvGeometry The fvGeometry
     * \param elemVolVars The element volume variables
     * \param elemFaceVars The element face variables
     * \param scvf The boundary sub control volume face
     */
    template<class ElementVolumeVariables, class ElementFluxVariablesCache>
    BoundaryFluxes neumann(const Element& element,
                           const FVElementGeometry& fvGeometry,
                           const ElementVolumeVariables& elemVolVars,
                           const ElementFluxVariablesCache& elemFluxVarsCache,
                           const SubControlVolumeFace& scvf) const
    {
        BoundaryFluxes values(0.0);
        if constexpr (!ParentType::isMomentumProblem())
        {
            // Continuity (volume flux) at every boundary face.
            // At no-slip walls faceVelocity = 0, so those faces contribute nothing.
            // At the inlet (momentum Dirichlet) this provides the inflow term; without
            // it the continuity equation has no inflow signal and the solver finds u=0.
            const auto rho = elemVolVars[scvf.insideScvIdx()].density();
            values[Indices::conti0EqIdx] =
                rho * this->faceVelocity(element, fvGeometry, scvf) * scvf.unitOuterNormal();

            // Advective outflow for phase field: let φ exit with the flow.
            // For any backflow at the outlet, impose gas phase (φ = -1).
            const auto phi = (values[Indices::conti0EqIdx] >= 0.0)
                ? elemVolVars[scvf.insideScvIdx()].phaseField()
                : -1.0;
            values[Indices::phaseFieldEqIdx] = phi * values[Indices::conti0EqIdx] / rho;

            if (isOnPillarSurface_(scvf.ipGlobal()))
            {
                // Contact angle BC on pillar walls (Jacqmin 2000):
                //   γε (∂φ/∂n) = f_w'(φ)   with f_w(φ) = -(γ cos θ)(3/4)(φ - φ³/3)
                //   → f_w'(φ) = -(γ cos θ)(3/4)(1 - φ²)
                // Natural BC (θ = 90°) gives f_w' = 0, recovering ∂φ/∂n = 0.
                const auto phi = elemVolVars[scvf.insideScvIdx()].phaseField();
                values[Indices::chemicalPotentialEqIdx] =
                    -scaledSurfaceTension_ * std::cos(contactAngle_) * (3.0/4.0) * (1.0 - phi*phi);
            }
        }

        // For Stokes (no inertia): zero Neumann = do-nothing (zero traction) at outlet.

        return values;
    }

    /*!
     * \brief Initial values for momentum DOFs (velocity = 0).
     */
    InitialValues initial(const SubControlVolume& scv) const
    {
        return InitialValues(0.0);
    }

    /*!
     * \brief Initial values for mass DOFs (vertices): tanh profile for phase field.
     *
     * The interface is placed at y = 0 (channel/box boundary).
     * phi = +1 (liquid/water) for y < 0, phi = -1 (gas/vapor) for y > 0.
     */
    InitialValues initial(const Vertex& vertex) const
    {
        if constexpr (ParentType::isMomentumProblem())
            return InitialValues(0.0);
        else
        {
            InitialValues values(0.0);
            const auto pos = vertex.geometry().corner(0)[dimWorld-1] + 0.1; // shift slightly so channel is free
            values[Indices::pressureIdx] = 0.0;
            values[Indices::phaseFieldIdx] = std::tanh(-pos / (std::sqrt(2.0)*interfaceThickness_));
            values[Indices::chemicalPotentialIdx] = 0.0;
            return values;
        }
    }

    /*!
     * \brief Returns the interface flux contribution due to the capillary forces
     */
    template<class Context>
    auto interfaceFlux(const Context& context) const
    {
        const auto& element = context.element();
        const auto& fvGeometry = context.fvGeometry();
        const auto& elemVolVars = context.elemVolVars();
        const auto& scvf = context.scvFace();
        const auto& fluxVarCache = context.elemFluxVarsCache()[scvf];
        const auto& shapeValues = fluxVarCache.shapeValues();

        const auto grads = this->couplingManager().gradients(element, fvGeometry, fluxVarCache.ipData());
        // const auto phi = this->couplingManager().values(element, fvGeometry, fluxVarCache.ipData())[CouplingManager::phaseFieldIdx];
        const auto& gradPhi = grads[CouplingManager::phaseFieldIdx];
        const auto& gradChemicalPotential = grads[CouplingManager::chemicalPotentialIdx];
        const auto normal = scvf.unitOuterNormal();
        std::decay_t<decltype(gradPhi)> v(0.0);
        for (const auto& localDof : localDofs(fvGeometry))
            v.axpy(shapeValues[localDof.index()][0], elemVolVars[localDof.index()].velocity());

        // Korteweg traction (AGG07 form):
        // T_cap n = γε[(∇φ ⊗ ∇φ)] n --> instead we do this as source term
        // T_diff n = v ⊗ (M ∇μ) 1/2 (ρ1 - ρ2) n if considering conservative form and inertia

        if (this->enableInertiaTerms())
        {
            return (
                // (scaledSurfaceTension_ * interfaceThickness_) * gradPhi * (gradPhi * normal)
                // - (0.5*gradPhi.two_norm2() + 0.25*(phi*phi - 1.0)*(phi*phi - 1.0)) * normal
                + this->mobility() * v * (gradChemicalPotential * normal) * this->mixtureDensityDerivative()
            ) * Extrusion::area(fvGeometry, scvf) * elemVolVars[scvf.insideScvIdx()].extrusionFactor();
        }
        // else
        // {
        //     return (
        //         (scaledSurfaceTension_ * interfaceThickness_) * gradPhi * (gradPhi * normal)
        //         - (0.5*gradPhi.two_norm2() + 0.25*(phi*phi - 1.0)*(phi*phi - 1.0)) * normal
        //     ) * Extrusion::area(fvGeometry, scvf) * elemVolVars[scvf.insideScvIdx()].extrusionFactor();
        // }

        return std::decay_t<decltype(gradPhi)>(0.0);
    }

    Scalar mixtureViscosity(const Scalar phaseField) const
    {
        const auto phi = std::clamp(phaseField, -1.0, 1.0);
        return 0.5 * ((1.0 + phi) * eta1_ + (1.0 - phi) * eta2_);
    }

    Scalar mixtureDensity(const Scalar phaseField) const
    {
        const auto phi = std::clamp(phaseField, -1.0, 1.0);
        return 0.5 * ((1.0 + phi) * rho1_ + (1.0 - phi) * rho2_);
    }

    Scalar mixtureDensityDerivative() const
    {
        return 0.5 * (rho1_  - rho2_);
    }

    template<class Solution, class GridVariables>
    void printMassBalanceSummary(const Solution& sol, const GridVariables& gridVariables) const
    {
        Scalar totalMassLiquid = 0.0;
        Scalar totalMassGas = 0.0;
        const auto& gg = this->gridGeometry();
        auto fvGeometry = localView(gg);
        for (const auto& element : elements(gg.gridView()))
        {
            fvGeometry.bind(element);
            for (const auto& scv : scvs(fvGeometry))
            {
                const auto phi = sol[scv.dofIndex()][Indices::phaseFieldIdx];
                const auto volume = scv.volume();
                // per-phase mass = phase density × volume fraction × cell volume
                totalMassLiquid += rho1_ * 0.5 * (1.0 + phi) * volume;
                totalMassGas    += rho2_ * 0.5 * (1.0 - phi) * volume;
            }
        }

        std::cout << "\033[1;36m[mass] total = " << totalMassLiquid + totalMassGas
                  << " kg/m  |  liquid = " << totalMassLiquid
                  << " kg/m  |  gas = " << totalMassGas
                  << " kg/m\033[0m" << std::endl;
    }

private:

    // Parabolic (Poiseuille) inlet profile with mean velocity inletVelocity_.
    // Channel occupies y ∈ [0, channelTop]; u(y) = 6·U_mean·(y - y_bot)·(y_top - y) / H².
    Scalar parabolicProfile_(const Scalar y) const
    {
        const Scalar yBot = 0.0; // channel floor (gas/liquid interface level)
        const Scalar yTop = this->gridGeometry().bBoxMax()[dimWorld-1];
        const Scalar H = yTop - yBot;
        return 6.0 * inletVelocity_ * (y - yBot) * (yTop - y) / (H * H);
    }

    bool isTop_(const GlobalPosition& pos, const Scalar eps = eps_) const
    { return pos[dimWorld-1] > this->gridGeometry().bBoxMax()[dimWorld-1] - eps; }

    bool isBottom_(const GlobalPosition& pos, const Scalar eps = eps_) const
    { return pos[dimWorld-1] < this->gridGeometry().bBoxMin()[dimWorld-1] + eps; }

    bool isInlet_(const GlobalPosition& pos) const
    { return pos[0] < this->gridGeometry().bBoxMin()[0] + eps_; }

    bool isOutlet_(const GlobalPosition& pos) const
    { return pos[0] > this->gridGeometry().bBoxMax()[0] - eps_; }

    bool isSide_(const GlobalPosition& pos) const
    { return isInlet_(pos) || isOutlet_(pos); }

    // A point is on a pillar surface if it lies within pillarRadius_ of any pillar center.
    // Tolerance slightly larger than the mesh size near pillars (lc_pillar = 0.001 in the .geo).
    bool isOnPillarSurface_(const GlobalPosition& pos, const Scalar tol = 2e-3) const
    {
        for (const auto& c : pillarCenters_)
            if ((pos - c).two_norm() < pillarRadius_ + tol)
                return true;
        return false;
    }

    static constexpr Scalar eps_ = 1e-6;

    Scalar energyScale_, mobility_, surfaceTension_, interfaceThickness_, scaledSurfaceTension_;

    // Fluid properties
    Scalar rho1_, rho2_;  // densities
    Scalar eta1_, eta2_;  // viscosities

    Scalar evaporationRate_;
    Scalar inletVelocity_;
    bool enableKortewegForce_;
    Scalar contactAngle_;
    Scalar pillarRadius_;
    std::vector<GlobalPosition> pillarCenters_;
};

} // end namespace Dumux

#endif
