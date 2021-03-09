// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \ingroup FreeflowNCModel
 *
 * \brief A single-phase, multi-component one-equation model
 *
 * For a detailed model decription see dumux/freeflow/compositional/navierstokesncmodel.hh
 */

#ifndef DUMUX_ONEEQ_NC_MODEL_HH
#define DUMUX_ONEEQ_NC_MODEL_HH

#include <dumux/common/properties.hh>
#include <dumux/freeflow/compositional/navierstokesncmodel.hh>
#include <dumux/freeflow/nonisothermal/iofields.hh>
#include <dumux/freeflow/rans/oneeq/model.hh>

#include "iofields.hh"

namespace Dumux {

///////////////////////////////////////////////////////////////////////////
// properties for the single-phase, multi-component one-equation model
///////////////////////////////////////////////////////////////////////////
namespace Properties {

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

// Create new type tags
namespace TTag {
//! The type tags for the single-phase, multi-component isothermal one-equation model
struct OneEqNC { using InheritsFrom = std::tuple<NavierStokesNC>; };
} // end namespace TTag

///////////////////////////////////////////////////////////////////////////
// default property values
///////////////////////////////////////////////////////////////////////////

/*!
 * \ingroup FreeflowNCModel
 * \brief Traits for the one-equation multi-component model
 *
 * \tparam dimension The dimension of the problem
 * \tparam nComp The number of components to be considered
 * \tparam useM Use molar or mass balances
 * \tparam replaceCompEqIdx The index of the component balance equation that should be replaced by a total mass/mole balance
 */
template<int dimension, int nComp, bool useMoles, int replaceCompEqIdx>
struct OneEqNCModelTraits : NavierStokesNCModelTraits<dimension, nComp, useMoles, replaceCompEqIdx>
{
    //! There are as many momentum balance equations as dimensions
    //! and as many balance equations as components.
    static constexpr int numEq() { return dimension+nComp+1; }

    //! The model does include a turbulence model
    static constexpr bool usesTurbulenceModel() { return true; }

    //! return the type of turbulence model used
    static constexpr auto turbulenceModel()
    { return TurbulenceModel::oneeq; }

    //! the indices
    using Indices = OneEqIndices<dimension, nComp>;
};

//!< states some specifics of the isothermal multi-component one-equation model
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::OneEqNC>
{
private:
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    static constexpr int dimension = GridView::dimension;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    static constexpr int numComponents = FluidSystem::numComponents;
    static constexpr bool useMoles = getPropValue<TypeTag, Properties::UseMoles>();
    static constexpr int replaceCompEqIdx = getPropValue<TypeTag, Properties::ReplaceCompEqIdx>();
public:
    using type = OneEqNCModelTraits<dimension, numComponents, useMoles, replaceCompEqIdx>;
};

//! Set the volume variables property
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::OneEqNC>
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
    using type = OneEqVolumeVariables<NCTraits<BaseTraits, DT>, NCVolVars>;
};

//! The local residual
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::OneEqNC>
{
private:
    using BaseLocalResidual = FreeflowNCResidual<TypeTag>;
public:
    using type = OneEqResidual<TypeTag, BaseLocalResidual>;
};

//! The flux variables
template<class TypeTag>
struct FluxVariables<TypeTag, TTag::OneEqNC>
{
private:
    using BaseFluxVariables = FreeflowNCFluxVariables<TypeTag>;
public:
    using type = OneEqFluxVariables<TypeTag, BaseFluxVariables>;
};

//! The specific I/O fields
template<class TypeTag>
struct IOFields<TypeTag, TTag::OneEqNC> { using type = FreeflowNCIOFields<OneEqIOFields, true/*turbulenceModel*/>; };

//////////////////////////////////////////////////////////////////////////
// Property values for non-isothermal multi-component one-equation model
//////////////////////////////////////////////////////////////////////////

// Create new type tags
namespace TTag {
//! The type tags for the single-phase, multi-component non-isothermal one-equation models
struct OneEqNCNI { using InheritsFrom = std::tuple<OneEqNC, NavierStokesNCNI>; };
} // end namespace TTag

//! The model traits of the non-isothermal model
template<class TypeTag>
struct ModelTraits<TypeTag, TTag::OneEqNCNI>
{
private:
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    static constexpr int dim = GridView::dimension;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    static constexpr int numComponents = FluidSystem::numComponents;
    static constexpr bool useMoles = getPropValue<TypeTag, Properties::UseMoles>();
    static constexpr int replaceCompEqIdx = getPropValue<TypeTag, Properties::ReplaceCompEqIdx>();
    using IsothermalModelTraits = OneEqNCModelTraits<dim, numComponents, useMoles, replaceCompEqIdx>;
public:
    using type = FreeflowNIModelTraits<IsothermalModelTraits>;
};

//! Set the volume variables property
template<class TypeTag>
struct VolumeVariables<TypeTag, TTag::OneEqNCNI>
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
    using type = OneEqVolumeVariables<NCNITraits<BaseTraits, DT>, NCNIVolVars>;
};

//! The local residual
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::OneEqNCNI>
{
private:
    using BaseLocalResidual = FreeflowNCResidual<TypeTag>;
public:
    using type = OneEqResidual<TypeTag, BaseLocalResidual>;
};

//! The flux variables
template<class TypeTag>
struct FluxVariables<TypeTag, TTag::OneEqNCNI>
{
private:
    using BaseFluxVariables = FreeflowNCFluxVariables<TypeTag>;
public:
    using type = OneEqFluxVariables<TypeTag, BaseFluxVariables>;
};

//! The specific I/O fields
template<class TypeTag>
struct IOFields<TypeTag, TTag::OneEqNCNI>
{
private:
    using IsothermalIOFields = FreeflowNCIOFields<OneEqIOFields, true/*turbulenceModel*/>;
public:
    using type = FreeflowNonIsothermalIOFields<IsothermalIOFields, true/*turbulenceModel*/>;
};

} // end namespace Properties
} // end namespace Dumux

#endif
