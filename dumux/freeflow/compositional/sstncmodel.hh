// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FreeflowNCModel
 *
 * \brief A single-phase, multi-component SST model
 *
 * For a detailed model description see dumux/freeflow/compositional/navierstokesncmodel.hh
 */

#ifndef DUMUX_SST_NC_MODEL_HH
#define DUMUX_SST_NC_MODEL_HH

#include <dumux/common/properties.hh>
#include <dumux/freeflow/compositional/navierstokesncmodel.hh>
#include <dumux/freeflow/nonisothermal/iofields.hh>
#include <dumux/freeflow/compositional/iofields.hh>
#include <dumux/freeflow/rans/twoeq/sst/model.hh>

namespace Dumux {

///////////////////////////////////////////////////////////////////////////
// properties for the single-phase, multi-component SST model
///////////////////////////////////////////////////////////////////////////
namespace Properties {

// Create new type tags
namespace TTag {
//! The type tags for the single-phase, multi-component isothermal SST model
struct SSTNC { using InheritsFrom = std::tuple<NavierStokesNC>; };
} // end namespace TTag

///////////////////////////////////////////////////////////////////////////
// default property values
///////////////////////////////////////////////////////////////////////////

/*!
 * \ingroup FreeflowNCModel
 * \brief Traits for the SST multi-component model
 *
 * \tparam dimension The dimension of the problem
 * \tparam nComp The number of components to be considered
 * \tparam useM Use molar or mass balances
 * \tparam replaceCompEqIdx The index of the component balance equation that should be replaced by a total mass/mole balance
 */
template<int dimension, int nComp, bool useMoles, int replaceCompEqIdx>
struct SSTNCModelTraits : NavierStokesNCModelTraits<dimension, nComp, useMoles, replaceCompEqIdx>
{
    //! There are as many momentum balance equations as dimensions
    //! and as many balance equations as components, plus two turbulence equations
    static constexpr int numEq() { return dimension+nComp+2; }

    //! The model does include a turbulence model
    static constexpr bool usesTurbulenceModel() { return true; }

    //! return the type of turbulence model used
    static constexpr auto turbulenceModel()
    { return TurbulenceModel::sst; }

    //! the indices
    using Indices = RANSTwoEqIndices<dimension, nComp>;
};

//!< states some specifics of the isothermal multi-component sst model
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::SSTNC>
{
private:
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    static constexpr int dimension = GridView::dimension;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    static constexpr int numComponents = FluidSystem::numComponents;
    static constexpr bool useMoles = getPropValue<TypeTag, Properties::UseMoles>();
    static constexpr int replaceCompEqIdx = getPropValue<TypeTag, Properties::ReplaceCompEqIdx>();
public:
    using type = SSTNCModelTraits<dimension, numComponents, useMoles, replaceCompEqIdx>;
};

//! Set the volume variables property
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::SSTNC>
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FSY = GetPropType<TypeTag, Properties::FluidSystem>;
    using FST = GetPropType<TypeTag, Properties::FluidState>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
    static_assert(FSY::numComponents == MT::numFluidComponents(), "Number of components mismatch between model and fluid system");
    static_assert(FST::numComponents == MT::numFluidComponents(), "Number of components mismatch between model and fluid state");
    static_assert(FSY::numPhases == MT::numFluidPhases(), "Number of phases mismatch between model and fluid system");
    static_assert(FST::numPhases == MT::numFluidPhases(), "Number of phases mismatch between model and fluid state");
    using BaseTraits = NavierStokesVolumeVariablesTraits<PV, FSY, FST, MT>;

    using DT = GetPropType<TypeTag, Properties::MolecularDiffusionType>;
    template<class BaseTraits, class DT>
    struct NCTraits : public BaseTraits { using DiffusionType = DT; };
    using NCVolVars = FreeflowNCVolumeVariables<NCTraits<BaseTraits, DT>>;
public:
    using type = SSTVolumeVariables<NCTraits<BaseTraits, DT>, NCVolVars>;
};

//! The local residual
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::SSTNC>
{
private:
    using BaseLocalResidual = FreeflowNCResidual<TypeTag>;
public:
    using type = SSTResidual<TypeTag, BaseLocalResidual>;
};

//! The flux variables
template<class TypeTag>
struct FluxVariables<TypeTag, TTag::SSTNC>
{
private:
    using BaseFluxVariables = FreeflowNCFluxVariables<TypeTag>;
public:
    using type = SSTFluxVariables<TypeTag, BaseFluxVariables>;
};

//! The specific I/O fields
template<class TypeTag>
struct IOFields<TypeTag, TTag::SSTNC> { using type = FreeflowNCIOFields<SSTIOFields, true/*turbulenceModel*/>; };

//////////////////////////////////////////////////////////////////////////
// Property values for non-isothermal multi-component SST model
//////////////////////////////////////////////////////////////////////////

// Create new type tags
namespace TTag {
//! The type tags for the single-phase, multi-component non-isothermal SST models
struct SSTNCNI { using InheritsFrom = std::tuple<SSTNC, NavierStokesNCNI>; };
} // end namespace TTag

//! The model traits of the non-isothermal model
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::SSTNCNI>
{
private:
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    static constexpr int dim = GridView::dimension;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    static constexpr int numComponents = FluidSystem::numComponents;
    static constexpr bool useMoles = getPropValue<TypeTag, Properties::UseMoles>();
    static constexpr int replaceCompEqIdx = getPropValue<TypeTag, Properties::ReplaceCompEqIdx>();
    using IsothermalModelTraits = SSTNCModelTraits<dim, numComponents, useMoles, replaceCompEqIdx>;

public:
    using type = FreeflowNIModelTraits<IsothermalModelTraits>;
};

//! Set the volume variables property
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::SSTNCNI>
{
private:
    using PV = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using FSY = GetPropType<TypeTag, Properties::FluidSystem>;
    using FST = GetPropType<TypeTag, Properties::FluidState>;
    using MT = GetPropType<TypeTag, Properties::ModelTraits>;
    static_assert(FSY::numComponents == MT::numFluidComponents(), "Number of components mismatch between model and fluid system");
    static_assert(FST::numComponents == MT::numFluidComponents(), "Number of components mismatch between model and fluid state");
    static_assert(FSY::numPhases == MT::numFluidPhases(), "Number of phases mismatch between model and fluid system");
    static_assert(FST::numPhases == MT::numFluidPhases(), "Number of phases mismatch between model and fluid state");
    using BaseTraits = NavierStokesVolumeVariablesTraits<PV, FSY, FST, MT>;

    using DT = GetPropType<TypeTag, Properties::MolecularDiffusionType>;
    template<class BaseTraits, class DT>
    struct NCNITraits : public BaseTraits { using DiffusionType = DT; };
    using NCNIVolVars = FreeflowNCVolumeVariables<NCNITraits<BaseTraits, DT>>;
public:
    using type = SSTVolumeVariables<NCNITraits<BaseTraits, DT>, NCNIVolVars>;
};

//! The local residual
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::SSTNCNI>
{
private:
    using BaseLocalResidual = FreeflowNCResidual<TypeTag>;
public:
    using type = SSTResidual<TypeTag, BaseLocalResidual>;
};

//! The flux variables
template<class TypeTag>
struct FluxVariables<TypeTag, TTag::SSTNCNI>
{
private:
    using BaseFluxVariables = FreeflowNCFluxVariables<TypeTag>;
public:
    using type = SSTFluxVariables<TypeTag, BaseFluxVariables>;
};

//! The specific I/O fields
template<class TypeTag>
struct IOFields<TypeTag, TTag::SSTNCNI>
{
private:
    using IsothermalIOFields = FreeflowNCIOFields<SSTIOFields, true/*turbulenceModel*/>;
public:
    using type = FreeflowNonIsothermalIOFields<IsothermalIOFields, true/*turbulenceModel*/>;
};

} // end namespace Properties
} // end namespace Dumux

#endif
