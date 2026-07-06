// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Hele-Shaw (Darcy) Cahn-Hilliard two-phase model for Saffman-Taylor instability.
 *
 * Three coupled equations for pressure p, phase field φ, chemical potential μ:
 *
 *  (1) Darcy continuity:    ∇·(λ(φ)(∇p − μ∇φ − ρ_mix(φ)g)) = 0
 *                           λ(φ) = K/η_mix(φ),  K = b²/12 (Hele-Shaw permeability)
 *                           μ∇φ is the Korteweg capillary body force
 *
 *  (2) Phase field (CH):    ∂φ/∂t + ∇·(v φ) = ∇·(M∇μ)
 *                           v = -λ(φ)(∇p − μ∇φ − ρ_mix g)
 *
 *  (3) Chemical potential:  μ = (γ/ε)φ(φ²-1) - γε∇²φ
 *
 * Discretized with CVFE (Box method). No coupling manager required.
 */
#ifndef DUMUX_PHASEFIELD_HELESHAW_2P_MODEL_HH
#define DUMUX_PHASEFIELD_HELESHAW_2P_MODEL_HH

#include <algorithm>

#include <dune/common/fvector.hh>
#include <dumux/common/math.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/defaultlocaloperator.hh>

namespace Dumux {

// ============================================================
// 1. Volume variables
// ============================================================

template<class Traits>
class HeleShawTwoPVolumeVariables
{
    using Scalar = typename Traits::PrimaryVariables::value_type;
    static_assert(Traits::PrimaryVariables::dimension == Traits::ModelTraits::numEq());
public:
    using PrimaryVariables = typename Traits::PrimaryVariables;
    using Indices = typename Traits::ModelTraits::Indices;

    template<class ElementSolution, class Problem, class Element, class SubControlVolume>
    void update(const ElementSolution& elemSol,
                const Problem& problem,
                const Element& element,
                const SubControlVolume& scv)
    {
        priVars_ = elemSol[scv.indexInElement()];
        const auto phi = priVars_[Indices::phaseFieldIdx];
        rhoMix_ = problem.mixtureDensity(phi);
        etaMix_ = problem.mixtureViscosity(phi);
    }

    Scalar pressure() const        { return priVars_[Indices::pressureIdx]; }
    Scalar phaseField() const      { return priVars_[Indices::phaseFieldIdx]; }
    Scalar chemicalPotential() const { return priVars_[Indices::chemPotIdx]; }
    Scalar mixtureDensity() const  { return rhoMix_; }
    Scalar mixtureViscosity() const{ return etaMix_; }

    Scalar priVar(int pvIdx) const { return priVars_[pvIdx]; }
    const PrimaryVariables& priVars() const { return priVars_; }
    Scalar extrusionFactor() const { return 1.0; }

private:
    PrimaryVariables priVars_;
    Scalar rhoMix_, etaMix_;
};

// ============================================================
// 2. Local residual
// ============================================================

template<class TypeTag>
class HeleShawTwoPLocalResidual
: public DiscretizationDefaultLocalOperator<TypeTag>
{
    using ParentType = DiscretizationDefaultLocalOperator<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using NumEqVector = Dumux::NumEqVector<GetPropType<TypeTag, Properties::PrimaryVariables>>;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using VolumeVariables = typename GridVariables::GridVolumeVariables::VolumeVariables;
    using ElementVolumeVariables = typename GridVariables::GridVolumeVariables::LocalView;
    using ElementFluxVariablesCache = typename GridVariables::GridFluxVariablesCache::LocalView;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename ModelTraits::Indices;

    static constexpr int dimWorld = GridView::dimensionworld;

public:
    using ParentType::ParentType;

    NumEqVector computeStorage(const Problem&,
                               const SubControlVolume&,
                               const VolumeVariables& volVars) const
    {
        NumEqVector storage(0.0);
        // pressure equation is elliptic (incompressibility constraint): no storage
        storage[Indices::continuityEqIdx] = 0.0;
        // phase field is transported: ∂φ/∂t
        storage[Indices::phaseFieldEqIdx] = volVars.phaseField();
        // chemical potential equation is algebraic: no storage
        storage[Indices::chemPotEqIdx] = 0.0;
        return storage;
    }

    NumEqVector computeFlux(const Problem& problem,
                            const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolumeFace& scvf,
                            const ElementFluxVariablesCache& elemFluxVarsCache) const
    {
        static_assert(DiscretizationMethods::isCVFE<typename GridGeometry::DiscretizationMethod>,
            "HeleShawTwoPLocalResidual requires a CVFE discretization");

        NumEqVector flux(0.0);

        const auto& fluxVarCache = elemFluxVarsCache[scvf];
        Dune::FieldVector<Scalar, dimWorld> gradP(0.0);
        Dune::FieldVector<Scalar, dimWorld> gradPhi(0.0);
        Dune::FieldVector<Scalar, dimWorld> gradMu(0.0);

        for (const auto& scv : scvs(fvGeometry))
        {
            const auto& vv = elemVolVars[scv];
            const auto& gradN = fluxVarCache.gradN(scv.indexInElement());
            gradP.axpy(vv.pressure(), gradN);
            gradPhi.axpy(vv.phaseField(), gradN);
            gradMu.axpy(vv.chemicalPotential(), gradN);
        }

        // Mixture properties at the face integration point (arithmetic average)
        const auto& inside  = elemVolVars[scvf.insideScvIdx()];
        const auto& outside = elemVolVars[scvf.outsideScvIdx()];
        const auto rhoMixIP = 0.5 * (inside.mixtureDensity()   + outside.mixtureDensity());
        const auto etaMixIP = 0.5 * (inside.mixtureViscosity()  + outside.mixtureViscosity());
        const auto muIP     = 0.5 * (inside.chemicalPotential() + outside.chemicalPotential());

        // Darcy mobility λ = K / η_mix
        const auto lambda = problem.permeability() / etaMixIP;
        const auto g = problem.gravity();
        const auto n = scvf.unitOuterNormal();
        const auto A = scvf.area();

        // Darcy volume flux: Q = -λ(∇p − μ∇φ − ρ_mix g)·n * A
        // The μ∇φ term is the Korteweg capillary body force (surface tension contribution).
        const Scalar darcyFlux = -lambda * (gradP * n - muIP*(gradPhi * n) - rhoMixIP*(g * n)) * A;

        // (1) Darcy continuity: ∇·v = 0  →  flux = Darcy volume flux
        flux[Indices::continuityEqIdx] = darcyFlux;

        // (2) Phase field: upwind advection + Cahn-Hilliard diffusion -M∇μ·n
        const auto phiUpwind = (darcyFlux >= 0.0) ? inside.phaseField() : outside.phaseField();
        flux[Indices::phaseFieldEqIdx] = phiUpwind * darcyFlux
                                       - problem.chMobility() * (gradMu * n) * A;

        // (3) Chemical potential: -γε ∇φ·n  (the Laplacian term in the CH equation)
        flux[Indices::chemPotEqIdx] = -problem.surfaceTensionCoeff() * (gradPhi * n) * A;

        return flux;
    }

    NumEqVector computeSource(const Problem& problem,
                              const Element& element,
                              const FVElementGeometry& fvGeometry,
                              const ElementVolumeVariables& elemVolVars,
                              const SubControlVolume& scv) const
    {
        NumEqVector source(0.0);
        source += problem.source(element, fvGeometry, elemVolVars, scv);
        return source;
    }
};

// ============================================================
// 3. IOFields
// ============================================================

struct HeleShawTwoPIOFields
{
    template<class OutputModule>
    static void initOutputModule(OutputModule& out)
    {
        out.addVolumeVariable([](const auto& v){ return v.pressure(); },        "p");
        out.addVolumeVariable([](const auto& v){ return v.phaseField(); },      "phi");
        out.addVolumeVariable([](const auto& v){ return v.chemicalPotential();}, "mu");
        out.addVolumeVariable([](const auto& v){ return v.mixtureDensity(); },  "rho");
        out.addVolumeVariable([](const auto& v){ return v.mixtureViscosity(); },"eta");
    }
};

} // end namespace Dumux

// ============================================================
// 4. TypeTag and property specializations
// ============================================================

namespace Dumux::Properties {

namespace TTag {
struct HeleShawTwoP {};
} // end namespace TTag

template<class TypeTag>
struct LocalResidual<TypeTag, TTag::HeleShawTwoP>
{ using type = HeleShawTwoPLocalResidual<TypeTag>; };

template<class TypeTag>
struct Scalar<TypeTag, TTag::HeleShawTwoP>
{ using type = double; };

template<class TypeTag>
struct ModelTraits<TypeTag, TTag::HeleShawTwoP>
{
    struct type
    {
        struct Indices
        {
            // primary variable indices
            static constexpr int pressureIdx   = 0;
            static constexpr int phaseFieldIdx = 1;
            static constexpr int chemPotIdx    = 2;

            // equation indices (1:1 with primary variables here)
            static constexpr int continuityEqIdx = 0;
            static constexpr int phaseFieldEqIdx = 1;
            static constexpr int chemPotEqIdx    = 2;
        };

        static constexpr int numEq() { return 3; }
    };
};

template<class TypeTag>
struct PrimaryVariables<TypeTag, TTag::HeleShawTwoP>
{
    using type = Dune::FieldVector<
        GetPropType<TypeTag, Properties::Scalar>,
        GetPropType<TypeTag, Properties::ModelTraits>::numEq()
    >;
};

template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::HeleShawTwoP>
{
    struct Traits
    {
        using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
        using ModelTraits      = GetPropType<TypeTag, Properties::ModelTraits>;
    };
    using type = HeleShawTwoPVolumeVariables<Traits>;
};

template<class TypeTag>
struct IOFields<TypeTag, TTag::HeleShawTwoP>
{ using type = HeleShawTwoPIOFields; };

} // end namespace Dumux::Properties

#endif
