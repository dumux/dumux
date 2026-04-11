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
    using Tensor = Dune::FieldMatrix<Scalar, dim, dimWorld>;

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
     * \brief Specifies the type of boundary condition for a degree of freedom.
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;

        // Inlet and outlet: no-flow for phase field, Dirichlet velocity for momentum
        values.setAllNeumann();

        // Top boundary: free surface with evaporation
        if (globalPos[1] > 0.08) // near top
            values.setAllNeumann();

        return values;
    }

    /*!
     * \brief Specifies Dirichlet boundary values.
     */
    DirichletValues dirichletAtPos(const GlobalPosition &globalPos) const
    {
        return DirichletValues(0.0);
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
            if (isTop_(scv.dofPosition()))
                values.setAllDirichlet();
            else
                values.setAllNeumann();
        }
        else
            values.setAllNeumann();

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
        }
        else
        {
            const auto ip = ipData(fvGeometry, scv);
            const auto grads = this->couplingManager().gradients(element, fvGeometry, ip);
            const auto mu = this->couplingManager().values(element, fvGeometry, ip)[CouplingManager::chemicalPotentialIdx];
            source = mu * grads[CouplingManager::phaseFieldIdx];
            source -= grads[CouplingManager::pressureIdx];
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
            if (isSide_(scvf.ipGlobal()))
            {
                values[Indices::conti0EqIdx] = 0.0;
                values[Indices::phaseFieldEqIdx] = 0.0;
            }

            else if (isBottom_(scvf.ipGlobal()))
            {
                values[Indices::conti0EqIdx] = this->faceVelocity(element, fvGeometry, scvf)*scvf.unitOuterNormal();
                values[Indices::phaseFieldEqIdx] = 0.0;
            }
        }
        else
        {
            if (isSide_(scvf.ipGlobal()))
            {
                const auto& fluxVarCache = elemFluxVarsCache[scvf];
                if (this->enableInertiaTerms())
                {
                    // advective term: vv*n
                    const auto elemSol = elementSolution(element, elemVolVars, fvGeometry);
                    const auto v = evalSolution(element, element.geometry(), fvGeometry.gridGeometry(), elemSol, scvf.ipGlobal());
                    const auto rho = this->couplingManager().density(element, fvGeometry, fluxVarCache.ipData());
                    values.axpy(rho*(v*scvf.unitOuterNormal()), v);
                }

                // stress tensor
                Tensor gradV(0.0);
                for (const auto& localDof : localDofs(fvGeometry))
                {
                    const auto& volVars = elemVolVars[localDof];
                    for (int dir = 0; dir < dim; ++dir)
                        gradV[dir].axpy(volVars.velocity(dir), fluxVarCache.gradN(localDof.index()));
                }

                // get viscosity from the problem
                const auto mu = this->couplingManager().effectiveViscosity(element, fvGeometry, fluxVarCache.ipData());
                // compute -mu*gradV*n*dA
                BoundaryFluxes momFlux = mv(gradV + getTransposed(gradV), scvf.unitOuterNormal());
                momFlux *= -mu;

                const auto pressure = this->couplingManager().pressure(element, fvGeometry, fluxVarCache.ipData());
                momFlux += pressure * scvf.unitOuterNormal();

                const auto normal = scvf.unitOuterNormal();
                values += (momFlux * normal) * normal;
            }

            else if (isBottom_(scvf.ipGlobal()))
            {
                const auto& fluxVarCache = elemFluxVarsCache[scvf];
                if (this->enableInertiaTerms())
                {
                    // advective term: vv*n
                    const auto elemSol = elementSolution(element, elemVolVars, fvGeometry);
                    const auto v = evalSolution(element, element.geometry(), fvGeometry.gridGeometry(), elemSol, scvf.ipGlobal());
                    const auto rho = this->couplingManager().density(element, fvGeometry, fluxVarCache.ipData());
                    values.axpy(rho*(v*scvf.unitOuterNormal()), v);
                }

                // zero pressure at bottom
            }
        }

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

        InitialValues values(0.0);
        const auto& pos = vertex.geometry().corner(0);
        values[Indices::phaseFieldIdx] = std::tanh(-pos[dimWorld-1] / (std::sqrt(2.0)*interfaceThickness_));
        values[Indices::chemicalPotentialIdx] = 0.0;
        return values;
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

private:

    bool isTop_(const GlobalPosition& pos, const Scalar eps = eps_) const
    { return pos[dimWorld-1] > this->gridGeometry().bBoxMax()[dimWorld-1] - eps; }

    bool isBottom_(const GlobalPosition& pos, const Scalar eps = eps_) const
    { return pos[dimWorld-1] < this->gridGeometry().bBoxMin()[dimWorld-1] + eps; }

    bool isSide_(const GlobalPosition& pos) const
    { return pos[0] < this->gridGeometry().bBoxMin()[0] + eps_ ||
             pos[0] > this->gridGeometry().bBoxMax()[0] - eps_; }

    static constexpr Scalar eps_ = 1e-6;

    Scalar energyScale_, mobility_, surfaceTension_, interfaceThickness_, scaledSurfaceTension_;

    // Fluid properties
    Scalar rho1_, rho2_;  // densities
    Scalar eta1_, eta2_;  // viscosities

    Scalar evaporationRate_;
};

} // end namespace Dumux

#endif
