// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Saffman-Taylor fingering test for the Hele-Shaw Darcy-Cahn-Hilliard model.
 *
 * Domain: [0, L] × [0, H], flow from left (high pressure) to right (low pressure).
 * Less viscous fluid (φ=+1) displaces more viscous fluid (φ=-1).
 * The interface is initially perturbed to trigger fingering.
 */
#ifndef DUMUX_TEST_FREEFLOW_HELESHAW_2P_PROBLEM_HH
#define DUMUX_TEST_FREEFLOW_HELESHAW_2P_PROBLEM_HH

#include <cmath>
#include <algorithm>

#include <dune/common/fvector.hh>

#include <dumux/common/boundarytypes.hh>
#include <dumux/discretization/scvandscvfiterators.hh>
#include <dumux/common/fvproblem.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

namespace Dumux {

template<class TypeTag>
class HeleShawTwoPTestProblem : public FVProblem<TypeTag>
{
    using ParentType = FVProblem<TypeTag>;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    using GravityVector = Dune::FieldVector<Scalar, GridView::dimensionworld>;

public:
    HeleShawTwoPTestProblem(std::shared_ptr<const GridGeometry> gg)
    : ParentType(gg)
    {
        pInlet_ = getParam<Scalar>("Problem.PressureInlet");
        pOutlet_ = getParam<Scalar>("Problem.PressureOutlet", 0.0);

        rho1_ = getParam<Scalar>("Problem.DensityInvading", 1000.0);
        rho2_ = getParam<Scalar>("Problem.DensityReceding",  1000.0);
        eta1_ = getParam<Scalar>("Problem.ViscosityInvading", 1.0);
        eta2_ = getParam<Scalar>("Problem.ViscosityReceding", 10.0);

        permeability_ = getParam<Scalar>("Problem.Permeability");

        const auto sigma = getParam<Scalar>("Problem.SurfaceTension");
        epsilon_  = getParam<Scalar>("Problem.InterfaceThickness");
        chMobility_ = getParam<Scalar>("Problem.CHMobility", epsilon_ * epsilon_);

        // γ = 3σ/(2√2),  surfaceTensionCoeff = γε,  energyScale = γ/ε
        const auto gamma = sigma * 3.0 / (2.0 * std::sqrt(2.0));
        surfaceTensionCoeff_ = gamma * epsilon_;
        energyScale_ = gamma / epsilon_;

        x0_        = getParam<Scalar>("Problem.InterfacePosition");
        delta_     = getParam<Scalar>("Problem.InterfacePerturbation", 0.0);
        nModes_    = getParam<int>("Problem.PerturbationModes", 1);
        baseMode_  = getParam<int>("Problem.PerturbationBaseMode", 1);

        const auto g = getParam<Scalar>("Problem.Gravity", 0.0);
        gravity_ = { 0.0, -g };
    }

    // ------- boundary conditions -------

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& pos) const
    {
        BoundaryTypes values;
        if (isInlet_(pos) || isOutlet_(pos))
        {
            // pressure is Dirichlet, phase field and chemical potential are free (Neumann)
            values.setDirichlet(Indices::pressureIdx);
            values.setNeumann(Indices::phaseFieldEqIdx);
            values.setNeumann(Indices::chemPotEqIdx);
        }
        else
            values.setAllNeumann();  // top/bottom walls: no-flow for all equations

        return values;
    }

    PrimaryVariables dirichletAtPos(const GlobalPosition& pos) const
    {
        PrimaryVariables values(0.0);
        values[Indices::pressureIdx] = isInlet_(pos) ? pInlet_ : pOutlet_;
        return values;
    }

    // Neumann fluxes: the key is the advective outflow of φ at inlet/outlet.
    // Without this, the outlet DOF has no outflux and φ creeps toward wrong values.
    template<class EVV, class EFVC, class SCVF>
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const EVV& elemVolVars,
                        const EFVC& elemFluxVarsCache,
                        const SCVF& scvf) const
    {
        NumEqVector flux(0.0);
        const auto& pos = scvf.ipGlobal();

        // Top/bottom walls: no-flow for all equations
        if (!isInlet_(pos) && !isOutlet_(pos))
            return flux;

        // Reconstruct the Darcy velocity at this boundary face.
        // elemFluxVarsCache[scvf] holds shape-function gradients even for boundary faces.
        static constexpr int dw = GridView::dimensionworld;
        const auto& cache = elemFluxVarsCache[scvf];
        Dune::FieldVector<Scalar, dw> gradP(0.0), gradPhi(0.0);
        for (const auto& scv : scvs(fvGeometry))
        {
            const auto& vv = elemVolVars[scv];
            gradP.axpy(vv.pressure(),    cache.gradN(scv.indexInElement()));
            gradPhi.axpy(vv.phaseField(), cache.gradN(scv.indexInElement()));
        }

        const auto& ins = elemVolVars[scvf.insideScvIdx()];
        const Scalar lambda = permeability_ / ins.mixtureViscosity();
        const auto   n      = scvf.unitOuterNormal();
        const Scalar muIP   = ins.chemicalPotential();
        const Scalar rhoIP  = ins.mixtureDensity();

        // Normal Darcy velocity (per unit area; assembly multiplies by area)
        const Scalar vn = -lambda * (gradP*n - muIP*(gradPhi*n) - rhoIP*(gravity_*n));

        // Phase field: upwind from inside for outflow, prescribed BC value for inflow.
        // Inlet carries φ=+1 (invading fluid), outlet carries whatever is inside.
        const Scalar phiBC    = isInlet_(pos) ? Scalar(1.0) : Scalar(-1.0);
        const Scalar phiUpwind = (vn < 0.0) ? phiBC : ins.phaseField();

        // Returns flux per unit area (no area factor here; assembly handles that)
        flux[Indices::phaseFieldEqIdx] = phiUpwind * vn;
        // μ: natural BC ∂μ/∂n = 0 → zero CH diffusion through walls  ✓
        // p: Dirichlet at inlet/outlet, so this Neumann value is never used for eq 0  ✓

        return flux;
    }

    // ------- initial conditions -------

    PrimaryVariables initialAtPos(const GlobalPosition& pos) const
    {
        PrimaryVariables values(0.0);

        // Sinusoidally perturbed interface position x_int(y)
        const auto H = this->gridGeometry().bBoxMax()[1];
        const auto L = this->gridGeometry().bBoxMax()[0];
        auto xInt = x0_;
        for (int k = baseMode_; k < baseMode_ + nModes_; ++k)
            xInt += delta_ * std::cos(k * 2.0 * M_PI * pos[1] / H);

        // Smooth tanh profile: φ=+1 left of interface (invading), φ=-1 right (receding)
        values[Indices::phaseFieldIdx] = std::tanh((xInt - pos[0]) / (std::sqrt(2.0) * epsilon_));
        values[Indices::chemPotIdx]    = 0.0;

        // Linear pressure from inlet to outlet as initial guess
        values[Indices::pressureIdx] = pInlet_ + (pOutlet_ - pInlet_) * (pos[0] / L);

        return values;
    }

    // ------- source term -------

    template<class ElementVolumeVariables>
    NumEqVector source(const Element&,
                       const FVElementGeometry&,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolume& scv) const
    {
        NumEqVector src(0.0);
        const auto& vv = elemVolVars[scv];
        const auto phi = vv.phaseField();
        // chemical potential equation: μ - (γ/ε)φ(φ²-1) = 0
        src[Indices::chemPotEqIdx] = vv.chemicalPotential() - energyScale_ * phi * (phi*phi - 1.0);
        return src;
    }

    // ------- fluid properties (called from VolumeVariables::update and LocalResidual) -------

    Scalar mixtureDensity(Scalar phi) const
    {
        phi = std::clamp(phi, Scalar(-1.0), Scalar(1.0));
        return 0.5 * ((1.0 + phi) * rho1_ + (1.0 - phi) * rho2_);
    }

    Scalar mixtureViscosity(Scalar phi) const
    {
        phi = std::clamp(phi, Scalar(-1.0), Scalar(1.0));
        return 0.5 * ((1.0 + phi) * eta1_ + (1.0 - phi) * eta2_);
    }

    Scalar permeability() const      { return permeability_; }
    Scalar chMobility() const        { return chMobility_; }
    Scalar surfaceTensionCoeff() const { return surfaceTensionCoeff_; }
    Scalar energyScale() const       { return energyScale_; }
    GravityVector gravity() const    { return gravity_; }

private:
    bool isInlet_(const GlobalPosition& pos) const
    { return pos[0] < this->gridGeometry().bBoxMin()[0] + 1e-8; }

    bool isOutlet_(const GlobalPosition& pos) const
    { return pos[0] > this->gridGeometry().bBoxMax()[0] - 1e-8; }

    Scalar pInlet_, pOutlet_;
    Scalar rho1_, rho2_, eta1_, eta2_;
    Scalar permeability_;
    Scalar epsilon_, chMobility_, surfaceTensionCoeff_, energyScale_;
    Scalar x0_, delta_;
    int nModes_, baseMode_;
    GravityVector gravity_;
};

} // end namespace Dumux

#endif
