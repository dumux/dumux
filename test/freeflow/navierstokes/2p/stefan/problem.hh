// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_TEST_FREEFLOW_TWOP_STEFAN_EVAPORATION_PROBLEM_HH
#define DUMUX_TEST_FREEFLOW_TWOP_STEFAN_EVAPORATION_PROBLEM_HH

#include <algorithm>
#include <cmath>

#include <dune/common/fvector.hh>

#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

namespace Dumux {

/*!
 * \brief Case A: planar phase-field evaporation against the Stefan sharp-interface limit.
 *
 * The strip is used as a 1D column: liquid is on the left (phi=+1), gas is on
 * the right (phi=-1). Evaporation is a diffuse-interface source that removes
 * liquid phase and adds vapor mass. The right boundary fixes the sub-saturated
 * far-field vapor concentration.
 */
template <class TypeTag, class BaseProblem>
class StefanEvaporationProblem : public BaseProblem
{
    using ParentType = BaseProblem;
    using BoundaryTypes = typename ParentType::BoundaryTypes;
    using BoundaryFluxes = typename ParentType::BoundaryFluxes;
    using DirichletValues = typename ParentType::DirichletValues;
    using InitialValues = typename ParentType::InitialValues;
    using Sources = typename ParentType::Sources;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVariables>::GridVolumeVariables::LocalView;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    static constexpr int dim = GridGeometry::GridView::dimension;
    static constexpr int dimWorld = GridGeometry::GridView::dimensionworld;

    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using Vertex = typename GridGeometry::GridView::template Codim<dim>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

public:
    StefanEvaporationProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                             std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry, couplingManager, ParentType::isMomentumProblem() ? "Momentum" : "Mass")
    {
        surfaceTension_ = getParam<Scalar>("Problem.SurfaceTension", 1e-3);
        interfaceThickness_ = getParam<Scalar>("Problem.InterfaceThickness", 1e-2);
        mobility_ = getParam<Scalar>("Problem.Mobility", 1e-5);

        const auto gamma = surfaceTension_ * 3.0/(2.0*std::sqrt(2.0));
        surfaceTensionCoeff_ = gamma*interfaceThickness_;
        energyScale_ = gamma/interfaceThickness_;

        rhoLiquid_ = getParam<Scalar>("Problem.Density1", 1.0);
        rhoGas_ = getParam<Scalar>("Problem.Density2", 1e-2);
        etaLiquid_ = getParam<Scalar>("Problem.Viscosity1", 1.0);
        etaGas_ = getParam<Scalar>("Problem.Viscosity2", 1e-2);

        interfacePosition0_ = getParam<Scalar>("Problem.InterfacePosition", 0.45);
        enableEvaporation_ = getParam<bool>("Problem.EnableEvaporation", true);
        evaporationRateCoeff_ = getParam<Scalar>("Problem.EvaporationRateCoeff", 5.0);
        vaporSatConc_ = getParam<Scalar>("Problem.VaporSatConc", 0.2);
        vaporFarFieldConc_ = getParam<Scalar>("Problem.VaporFarFieldConc", 0.05);
        vaporDiffusivity_ = getParam<Scalar>("Problem.VaporDiffusivity", 1e-2);

        // phi storage is not liquid volume fraction storage. Since alpha_l=(1+phi)/2,
        // dphi/dt = -2*m_dot/rho_l is the volume-consistent phase-change source.
        phaseChangeFactor_ = getParam<Scalar>("Problem.PhaseChangeFactor", 2.0);
        analyticTimeOffset_ = getParam<Scalar>("Problem.AnalyticTimeOffset", 0.0);
        stefanLambda_ = solveStefanLambda_();
    }

    BoundaryTypes boundaryTypes(const Element& element,
                                const SubControlVolume& scv) const
    {
        BoundaryTypes values;

        if constexpr (ParentType::isMomentumProblem())
            values.setAllDirichlet(); // Case A is diffusion-limited: keep the fluid at rest.
        else
        {
            values.setAllNeumann();
            const auto& pos = scv.dofPosition();

            if (isLeft_(pos) || isRight_(pos))
                values.setDirichlet(Indices::phaseFieldIdx);

            if (isRight_(pos))
                values.setDirichlet(Indices::vaporIdx);

            if (isPressureAnchor_(pos))
                values.setDirichlet(Indices::pressureIdx);
        }

        return values;
    }

    DirichletValues dirichletAtPos(const GlobalPosition& globalPos) const
    {
        DirichletValues values(0.0);

        if constexpr (!ParentType::isMomentumProblem())
        {
            if (isLeft_(globalPos))
                values[Indices::phaseFieldIdx] = 1.0;
            else if (isRight_(globalPos))
                values[Indices::phaseFieldIdx] = -1.0;

            if (isRight_(globalPos))
                values[Indices::vaporIdx] = vaporFarFieldConc_;
        }

        return values;
    }

    Sources source(const Element& element,
                   const FVElementGeometry& fvGeometry,
                   const ElementVolumeVariables& elemVolVars,
                   const SubControlVolume& scv) const
    {
        Sources source(0.0);

        if constexpr (!ParentType::isMomentumProblem())
        {
            const auto& volVars = elemVolVars[scv];
            const auto phi = volVars.phaseField();

            source[Indices::chemicalPotentialEqIdx] =
                volVars.chemicalPotential() - energyScale_*phi*(phi*phi - 1.0);

            if (enableEvaporation_)
            {
                const auto mDot = evaporationMassSource(phi, volVars.vaporConcentration());
                source[Indices::phaseFieldEqIdx] -= phaseChangeFactor_*mDot/rhoLiquid_;
                source[Indices::vaporEqIdx] += mDot;
            }
        }

        return source;
    }

    template<class ElementVolumeVariables, class ElementFluxVariablesCache>
    BoundaryFluxes neumann(const Element& element,
                           const FVElementGeometry& fvGeometry,
                           const ElementVolumeVariables& elemVolVars,
                           const ElementFluxVariablesCache& elemFluxVarsCache,
                           const typename FVElementGeometry::SubControlVolumeFace& scvf) const
    { return BoundaryFluxes(0.0); }

    InitialValues initialAtPos(const GlobalPosition& pos) const
    {
        InitialValues values(0.0);

        if constexpr (!ParentType::isMomentumProblem())
        {
            values[Indices::pressureIdx] = 0.0;
            values[Indices::phaseFieldIdx] =
                std::tanh((analyticalInterfacePosition(0.0) - pos[0])/(std::sqrt(2.0)*interfaceThickness_));
            values[Indices::chemicalPotentialIdx] = 0.0;
            values[Indices::vaporIdx] = analyticalVaporConcentration(pos, 0.0);
        }

        return values;
    }

    InitialValues initial(const SubControlVolume& scv) const
    { return initialAtPos(scv.dofPosition()); }

    InitialValues initial(const Vertex& vertex) const
    { return initialAtPos(vertex.geometry().corner(0)); }

    template<class Context>
    auto interfaceFlux(const Context& context) const
    {
        Dune::FieldVector<Scalar, dimWorld> zero(0.0);
        return zero;
    }

    Scalar mobility() const
    { return mobility_; }

    Scalar surfaceTension() const
    { return surfaceTensionCoeff_; }

    Scalar vaporDiffusivity() const
    { return vaporDiffusivity_; }

    Scalar mixtureViscosity(const Scalar phaseField) const
    {
        const auto phi = std::clamp(phaseField, Scalar(-1.0), Scalar(1.0));
        return 0.5*((1.0 + phi)*etaLiquid_ + (1.0 - phi)*etaGas_);
    }

    Scalar mixtureDensity(const Scalar phaseField) const
    {
        const auto phi = std::clamp(phaseField, Scalar(-1.0), Scalar(1.0));
        return 0.5*((1.0 + phi)*rhoLiquid_ + (1.0 - phi)*rhoGas_);
    }

    Scalar mixtureDensityDerivative() const
    { return 0.5*(rhoLiquid_ - rhoGas_); }

    Scalar liquidVolumeFraction(const Scalar phaseField) const
    { return 0.5*(1.0 + std::clamp(phaseField, Scalar(-1.0), Scalar(1.0))); }

    Scalar gasVolumeFraction(const Scalar phaseField) const
    { return 1.0 - liquidVolumeFraction(phaseField); }

    Scalar interfacialDelta(const Scalar phaseField) const
    {
        const auto g = std::max(Scalar(0), 1.0 - phaseField*phaseField);
        return g/(2.0*std::sqrt(2.0)*interfaceThickness_);
    }

    Scalar evaporationMassSource(const Scalar phaseField, const Scalar vaporConcentration) const
    {
        if (!enableEvaporation_)
            return 0.0;

        const auto drivingForce = std::max(Scalar(0), vaporSatConc_ - vaporConcentration);
        return evaporationRateCoeff_*drivingForce*interfacialDelta(phaseField);
    }

    Scalar analyticalInterfacePosition(const Scalar time) const
    {
        const auto t = time + analyticTimeOffset_;
        if (t <= 0.0)
            return interfacePosition0_;
        return interfacePosition0_ - 2.0*stefanLambda_*std::sqrt(vaporDiffusivity_*t);
    }

    Scalar analyticalMassFlux(const Scalar time) const
    {
        const auto t = time + analyticTimeOffset_;
        if (t <= 0.0)
            return 0.0;
        return rhoLiquid_*stefanLambda_*std::sqrt(vaporDiffusivity_/t);
    }

    Scalar analyticalVaporConcentration(const GlobalPosition& pos, const Scalar time) const
    {
        const auto t = time + analyticTimeOffset_;
        if (t <= 0.0)
            return pos[0] >= interfacePosition0_ ? vaporFarFieldConc_ : vaporSatConc_;

        const auto eta = (pos[0] - interfacePosition0_)/(2.0*std::sqrt(vaporDiffusivity_*t));
        return vaporFarFieldConc_
            + (vaporSatConc_ - vaporFarFieldConc_)*(1.0 - std::erf(eta))/(1.0 + std::erf(stefanLambda_));
    }

    Scalar liquidDensity() const
    { return rhoLiquid_; }

    Scalar vaporSatConc() const
    { return vaporSatConc_; }

    Scalar vaporFarFieldConc() const
    { return vaporFarFieldConc_; }

    Scalar stefanLambda() const
    { return stefanLambda_; }

private:
    bool isLeft_(const GlobalPosition& pos) const
    { return pos[0] < this->gridGeometry().bBoxMin()[0] + eps_; }

    bool isRight_(const GlobalPosition& pos) const
    { return pos[0] > this->gridGeometry().bBoxMax()[0] - eps_; }

    bool isBottom_(const GlobalPosition& pos) const
    { return pos[dimWorld-1] < this->gridGeometry().bBoxMin()[dimWorld-1] + eps_; }

    bool isPressureAnchor_(const GlobalPosition& pos) const
    { return isLeft_(pos) && isBottom_(pos); }

    Scalar solveStefanLambda_() const
    {
        const auto beta = (vaporSatConc_ - vaporFarFieldConc_)/(rhoLiquid_*std::sqrt(M_PI));
        if (beta <= 0.0)
            return 0.0;

        auto residual = [beta](const Scalar lambda)
        {
            return lambda*(1.0 + std::erf(lambda))*std::exp(lambda*lambda) - beta;
        };

        Scalar lo = 0.0;
        Scalar hi = 1.0;
        while (residual(hi) < 0.0 && hi < 16.0)
            hi *= 2.0;

        for (int i = 0; i < 80; ++i)
        {
            const auto mid = 0.5*(lo + hi);
            if (residual(mid) < 0.0)
                lo = mid;
            else
                hi = mid;
        }

        return 0.5*(lo + hi);
    }

    static constexpr Scalar eps_ = 1e-8;

    Scalar surfaceTension_;
    Scalar interfaceThickness_;
    Scalar mobility_;
    Scalar surfaceTensionCoeff_;
    Scalar energyScale_;

    Scalar rhoLiquid_;
    Scalar rhoGas_;
    Scalar etaLiquid_;
    Scalar etaGas_;

    Scalar interfacePosition0_;
    bool enableEvaporation_;
    Scalar evaporationRateCoeff_;
    Scalar vaporSatConc_;
    Scalar vaporFarFieldConc_;
    Scalar vaporDiffusivity_;
    Scalar phaseChangeFactor_;
    Scalar analyticTimeOffset_;
    Scalar stefanLambda_;
};

} // end namespace Dumux

#endif
