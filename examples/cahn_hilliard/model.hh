// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//

#ifndef DUMUX_EXAMPLES_CAHN_HILLIARD_MODEL_HH
#define DUMUX_EXAMPLES_CAHN_HILLIARD_MODEL_HH

// # Cahn-Hilliard equation model definition
//
// In the file `model.hh`, we define the model equations and
// set all default model properties. The setup consist of four steps:
// 1. Create a model type tag (used to specialize properties)
// 2. Define the volume variables class computing and storing variables for a control volume
// 3. Define the local residual class implementing the discrete equation
// 4. Specialize important properties of the model such that Dumux knows how to assemble the system matrix
//
// We implement the classes
//
// * `CahnHilliardModel` (step 1),
// * `CahnHilliardModelVolumeVariables` (step 2) and
// * `CahnHilliardModelLocalResidual` (step 3).
//
// __Table of contents__
//
// [TOC]
//
// We start in `model.hh` with the necessary header includes:
// [[details]] includes
#include <dune/common/fvector.hh>
#include <dumux/common/math.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/discretization/method.hh>
// [[/details]]
//
// ## 1. Property Tag
//
// The property tag is simply an empty struct with the name `CahnHilliardModel`.
//
// [[content]]
// [[codeblock]]
namespace Dumux::Properties::TTag {
struct CahnHilliardModel {};
} // end namespace Dumux::Properties::TTag
// [[/codeblock]]
// [[/content]]
//
// ## 2. The volume variables
//
// The volume variables store the control volume variables, both primary and secondary.
// Let's have a look at the class implementation.
//
// [[content]]
//
// We write out a full volume variable class compatible with the Dumux assembly.
// We have to fulfill the same interface as the [`BasicVolumeVariables`](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/blob/master/dumux/common/volumevariables.hh).
// To shorten the code, we could also inherit from `BasicVolumeVariables`. However,
// we want to show the entire interface here as a basis for the implementation of more complex variables.
// Note that in `update` we could, for example, compute some secondary variable that is a
// possibly nonlinear function of the primary variables. That secondary variable could be exposed
// via an interface and then used below in the `CahnHilliardModelLocalResidual` class implementation.
//
// [[codeblock]]
namespace Dumux {
template <class Traits>
class CahnHilliardModelVolumeVariables
{
    using Scalar = typename Traits::PrimaryVariables::value_type;
    static_assert(Traits::PrimaryVariables::dimension == Traits::ModelTraits::numEq());
public:
    //! export the type used for the primary variables
    using PrimaryVariables = typename Traits::PrimaryVariables;
    //! export the indices type
    using Indices = typename Traits::ModelTraits::Indices;
// [[/codeblock]]
    //
    // **Update variables:** The `update` function stores the local primary variables of the current solution and
    // potentially computes secondary variables. Secondary variables can be nonlinear functions
    // of the primary variables.
    //
    // [[codeblock]]
    //! Update all quantities for a given control volume
    template<class ElementSolution, class Problem, class Element, class SubControlVolume>
    void update(const ElementSolution& elemSol,
                const Problem& problem,
                const Element& element,
                const SubControlVolume& scv)
    {
        priVars_ = elemSol[scv.indexInElement()];
    }
    // [[/codeblock]]
    //
    // **Access functions:** Named and generic functions to access different primary variables.
    // The volumevariables are also expected to return an extrusion factor.
    // Here we don't extrude our domain and therefore return $1.0$.
    // [[codeblock]]
    Scalar concentration() const
    { return priVars_[Indices::concentrationIdx]; }

    Scalar chemicalPotential() const
    { return priVars_[Indices::chemicalPotentialIdx]; }

    Scalar priVar(const int pvIdx) const
    { return priVars_[pvIdx]; }

    const PrimaryVariables& priVars() const
    { return priVars_; }

    Scalar extrusionFactor() const
    { return 1.0; }

private:
    PrimaryVariables priVars_;
};
} // end namespace Dumux
// [[/codeblock]]
// [[/content]]
//
// ## 3. The local residual
//
// The local residual assembles the contribution to the residual for
// all degrees of freedom associated with an element.
// See [examples/diffusion/doc/model.hh](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/blob/master/examples/diffusion/doc/main.md)
// for a more detailed explanation for control-volume finite element method local residuals.
// Let's have a look at the class implementation.
//
// [[content]]
//
// The class `CahnHilliardModelLocalResidual` inherits from a base class set in
// the model properties, depending on the discretization scheme.
// See [examples/diffusion/doc/model.hh](https://git.iws.uni-stuttgart.de/dumux-repositories/dumux/-/blob/master/examples/diffusion/doc/main.md)
// for details on the `BaseLocalResidual`.
namespace Dumux {
template<class TypeTag>
class CahnHilliardModelLocalResidual
: public GetPropType<TypeTag, Properties::BaseLocalResidual>
{
    // [[exclude]]
    // the base local residual is selected depending on the chosen discretization scheme
    using ParentType = GetPropType<TypeTag, Properties::BaseLocalResidual>;
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
    // [[/exclude]]
    //
    // **Storage term:** The function `computeStorage` receives the volume variables
    // at the previous or current time step and computes the value of the storage terms.
    // In this case the mass balance equation is a conservation equation of the concentration and
    // the equation for the chemical potential does not have a storage term.
    //
    // [[codeblock]]
    NumEqVector computeStorage(const Problem& problem,
                               const SubControlVolume& scv,
                               const VolumeVariables& volVars) const
    {
        NumEqVector storage;
        storage[Indices::massBalanceEqIdx] = volVars.concentration();
        storage[Indices::chemicalPotentialEqIdx] = 0.0;
        return storage;
    }
    // [[/codeblock]]

    // **Flux term:** The function `computeFlux` computes the integrated
    // fluxes over a sub control volume face.
    //
    // ```math
    // \begin{aligned}
    // F_{K,\sigma,0} &= -M \vert \sigma \vert \sum_{B \in \mathcal{B}_K} \mu_{h,B} \nabla N_B \cdot\boldsymbol{n} \cr
    // F_{K,\sigma,1} &= -\gamma \vert \sigma \vert \sum_{B \in \mathcal{B}_K} c_{h,B} \nabla N_B \cdot\boldsymbol{n}
    // \end{aligned}
    // ````
    //
    // [[codeblock]]
    NumEqVector computeFlux(const Problem& problem,
                            const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolumeFace& scvf,
                            const ElementFluxVariablesCache& elemFluxVarsCache) const
    {
        static_assert(DiscretizationMethods::isCVFE<typename GridGeometry::DiscretizationMethod>,
            "This local residual is hard-coded to control-volume finite element schemes");

        const auto& fluxVarCache = elemFluxVarsCache[scvf];
        Dune::FieldVector<Scalar, dimWorld> gradConcentration(0.0);
        Dune::FieldVector<Scalar, dimWorld> gradChemicalPotential(0.0);
        for (const auto& scv : scvs(fvGeometry))
        {
            const auto& volVars = elemVolVars[scv];
            // v.axpy(a, w) means v += a*w
            gradConcentration.axpy(
                volVars.concentration(),
                fluxVarCache.gradN(scv.indexInElement())
            );
            gradChemicalPotential.axpy(
                volVars.chemicalPotential(),
                fluxVarCache.gradN(scv.indexInElement())
            );
        }

        NumEqVector flux;
        // Compute the flux with `vtmv` (vector transposed times matrix times vector).
        // The mobility and surface tension coefficients comes from the `problem` (see Part 2 of the example).
        flux[Indices::massBalanceEqIdx] = -1.0*vtmv(
            scvf.unitOuterNormal(), problem.mobility(), gradChemicalPotential
        )*scvf.area();
        flux[Indices::chemicalPotentialEqIdx] = -1.0*vtmv(
            scvf.unitOuterNormal(), problem.surfaceTension(), gradConcentration
        )*scvf.area();
        return flux;
    }
    // [[/codeblock]]

    // **Source term:** The function `computeSource` computes the source terms for a sub control volume.
    // We implement a model-specific source term for the chemical potential equation before
    // deferring further implementation to the problem where we add the derivative of the free
    // energy.
    //
    // [[codeblock]]
    NumEqVector computeSource(const Problem& problem,
                              const Element& element,
                              const FVElementGeometry& fvGeometry,
                              const ElementVolumeVariables& elemVolVars,
                              const SubControlVolume &scv) const
    {
        NumEqVector source(0.0);
        source[Indices::massBalanceEqIdx] = 0.0;
        source[Indices::chemicalPotentialEqIdx] = elemVolVars[scv].chemicalPotential();
        // add contributions from problem (e.g. double well potential)
        source += problem.source(element, fvGeometry, elemVolVars, scv);
        return source;
    }
};
} // end namespace Dumux
// [[/codeblock]]
// [[/content]]
//
// ## 4. The model properties/traits
//
// In the `Dumux::Properties` namespace, we specialize properties for
// the created type tag `CahnHilliardModel`.
//
// [[content]]
namespace Dumux::Properties {

// The type of the local residual is the class defined above.
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::CahnHilliardModel>
{ using type = CahnHilliardModelLocalResidual<TypeTag>; };

// We compute with double precision floating point numbers.
template<class TypeTag>
struct Scalar<TypeTag, TTag::CahnHilliardModel>
{ using type = double; };

// The model traits specify some information about our equation system.
// Here we have two equations. The indices allow to access primary variables
// and equations with named indices.
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::CahnHilliardModel>
{
    struct type
    {
        struct Indices
        {
            static constexpr int concentrationIdx = 0;
            static constexpr int chemicalPotentialIdx = 1;

            static constexpr int massBalanceEqIdx = 0;
            static constexpr int chemicalPotentialEqIdx = 1;
        };

        static constexpr int numEq() { return 2; }
    };
};

// The primary variable vector has entries of type `Scalar` and is
// as large as the number of equations (here 2) but we keep it general
// here by obtaining the number of equations from the `ModelTraits`.
template<class TypeTag>
struct PrimaryVariables<TypeTag, TTag::CahnHilliardModel>
{
    using type = Dune::FieldVector<
        GetPropType<TypeTag, Properties::Scalar>,
        GetPropType<TypeTag, Properties::ModelTraits>::numEq()
    >;
};

// Finally, the type of the volume variables is the class defined above.
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::CahnHilliardModel>
{
    struct Traits
    {
        using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
        using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    };

    using type = CahnHilliardModelVolumeVariables<Traits>;
};

} // end namespace Dumux::Properties
// [[/content]]
// [[exclude]]
#endif
// [[/exclude]]
