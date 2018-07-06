// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 * \brief A single-phase, multi-component k-omega model
 *
 * \copydoc Dumux::FreeflowNCModel
 */

#ifndef DUMUX_KOMEGA_NC_MODEL_HH
#define DUMUX_KOMEGA_NC_MODEL_HH

#include <dumux/common/properties.hh>
#include <dumux/freeflow/compositional/navierstokesncmodel.hh>
#include <dumux/freeflow/nonisothermal/vtkoutputfields.hh>
#include <dumux/freeflow/rans/twoeq/komega/model.hh>

#include "vtkoutputfields.hh"

namespace Dumux {

///////////////////////////////////////////////////////////////////////////
// properties for the single-phase, multi-component k-omega model
///////////////////////////////////////////////////////////////////////////
namespace Properties {

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tags for the single-phase, multi-component isothermal k-omega model
NEW_TYPE_TAG(KOmegaNC, INHERITS_FROM(NavierStokesNC));

///////////////////////////////////////////////////////////////////////////
// default property values
///////////////////////////////////////////////////////////////////////////

/*!
 * \ingroup FreeflowNCModel
 * \brief Traits for the low-Reynolds k-epsilon multi-component model
 *
 * \tparam dimension The dimension of the problem
 * \tparam nComp The number of components to be considered
 * \tparam fluidSystemPhaseIdx The the index of the phase used for the fluid system
 * \tparam replaceCompEqIdx The index of the component balance equation that should be replaced by a total mass/mole balance
 * \tparam useM Use molar or mass balances
 */
template<int dimension, int nComp, int fluidSystemPhaseIdx, int replaceCompEqIdx, bool useMoles>
struct KOmegaNCModelTraits : NavierStokesNCModelTraits<dimension, nComp, fluidSystemPhaseIdx, replaceCompEqIdx, useMoles>
{
    //! There are as many momentum balance equations as dimensions
    //! and as many balance equations as components.
    static constexpr int numEq() { return dimension+nComp+2; }

    //! The model does include a turbulence model
    static constexpr bool usesTurbulenceModel() { return true; }

    //! the indices
    using Indices = FreeflowNCIndices<dimension, numEq(), fluidSystemPhaseIdx, replaceCompEqIdx, KOmegaIndices<dimension, nComp, fluidSystemPhaseIdx>>;
};

//!< states some specifics of the isothermal multi-component low-Reynolds k-epsilon model
SET_PROP(KOmegaNC, ModelTraits)
{
private:
    using GridView = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::GridView;
    static constexpr int dimension = GridView::dimension;
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    static constexpr int numComponents = FluidSystem::numComponents;
    static constexpr int phaseIdx = GET_PROP_VALUE(TypeTag, PhaseIdx);
    static constexpr int replaceCompEqIdx = GET_PROP_VALUE(TypeTag, ReplaceCompEqIdx);
    static constexpr bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);
public:
    using type = KOmegaNCModelTraits<dimension, numComponents, phaseIdx, replaceCompEqIdx, useMoles>;
};

//! Set the volume variables property
SET_PROP(KOmegaNC, VolumeVariables)
{
private:
    using PV = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using FSY = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using FST = typename GET_PROP_TYPE(TypeTag, FluidState);
    using MT = typename GET_PROP_TYPE(TypeTag, ModelTraits);

    using Traits = NavierStokesVolumeVariablesTraits<PV, FSY, FST, MT>;
    using NCVolVars = FreeflowNCVolumeVariables<Traits>;
public:
    using type = KOmegaVolumeVariables<Traits, NCVolVars>;
};

//! The local residual
SET_PROP(KOmegaNC, LocalResidual)
{
private:
    using BaseLocalResidual = FreeflowNCResidual<TypeTag>;
public:
    using type = KOmegaResidual<TypeTag, BaseLocalResidual>;
};

//! The flux variables
SET_PROP(KOmegaNC, FluxVariables)
{
private:
    using BaseFluxVariables = FreeflowNCFluxVariables<TypeTag>;
public:
    using type = KOmegaFluxVariables<TypeTag, BaseFluxVariables>;
};

//! The specific vtk output fields
SET_PROP(KOmegaNC, VtkOutputFields)
{
private:
    using ModelTraits = typename GET_PROP_TYPE(TypeTag, ModelTraits);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    static constexpr int phaseIdx = GET_PROP_VALUE(TypeTag, PhaseIdx);
    using SinglePhaseVtkOutputFields = KOmegaVtkOutputFields<FVGridGeometry>;
public:
    using type = FreeflowNCVtkOutputFields<SinglePhaseVtkOutputFields, ModelTraits, FVGridGeometry, FluidSystem, phaseIdx>;
};

//////////////////////////////////////////////////////////////////////////
// Property values for non-isothermal multi-component k-omega model
//////////////////////////////////////////////////////////////////////////

//! The type tags for the single-phase, multi-component non-isothermal k-omega models
NEW_TYPE_TAG(KOmegaNCNI, INHERITS_FROM(NavierStokesNCNI));

//! The model traits of the non-isothermal model
SET_PROP(KOmegaNCNI, ModelTraits)
{
private:
    using GridView = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::GridView;
    static constexpr int dim = GridView::dimension;
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    static constexpr int numComponents = FluidSystem::numComponents;
    static constexpr int phaseIdx = GET_PROP_VALUE(TypeTag, PhaseIdx);
    static constexpr int replaceCompEqIdx = GET_PROP_VALUE(TypeTag, ReplaceCompEqIdx);
    static constexpr bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);
    using IsothermalModelTraits = KOmegaNCModelTraits<dim, numComponents, phaseIdx, replaceCompEqIdx, useMoles>;
public:
    using type = FreeflowNIModelTraits<IsothermalModelTraits>;
};

//! Set the volume variables property
SET_PROP(KOmegaNCNI, VolumeVariables)
{
private:
    using PV = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using FSY = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using FST = typename GET_PROP_TYPE(TypeTag, FluidState);
    using MT = typename GET_PROP_TYPE(TypeTag, ModelTraits);

    using Traits = NavierStokesVolumeVariablesTraits<PV, FSY, FST, MT>;
    using NCVolVars = FreeflowNCVolumeVariables<Traits>;
public:
    using type = KOmegaVolumeVariables<Traits, NCVolVars>;
};

//! The local residual
SET_PROP(KOmegaNCNI, LocalResidual)
{
private:
    using BaseLocalResidual = FreeflowNCResidual<TypeTag>;
public:
    using type = KOmegaResidual<TypeTag, BaseLocalResidual>;
};

//! The flux variables
SET_PROP(KOmegaNCNI, FluxVariables)
{
private:
    using BaseFluxVariables = FreeflowNCFluxVariables<TypeTag>;
public:
    using type = KOmegaFluxVariables<TypeTag, BaseFluxVariables>;
};

//! The specific vtk output fields
SET_PROP(KOmegaNCNI, VtkOutputFields)
{
private:
    using ModelTraits = typename GET_PROP_TYPE(TypeTag, ModelTraits);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    static constexpr int phaseIdx = GET_PROP_VALUE(TypeTag, PhaseIdx);
    using BaseVtkOutputFields = KOmegaVtkOutputFields<FVGridGeometry>;
    using NonIsothermalFields = FreeflowNonIsothermalVtkOutputFields<BaseVtkOutputFields, ModelTraits>;
public:
    using type = FreeflowNCVtkOutputFields<NonIsothermalFields, ModelTraits, FVGridGeometry, FluidSystem, phaseIdx>;
};

// \}
} // end namespace Properties
} // end namespace Dumux

#endif
