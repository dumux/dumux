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
 * \brief A single-phase, multi-component Reynolds-Averaged Navier-Stokes 0-Eq. model
 *
 * \copydoc Dumux::FreeflowNCModel
 */

#ifndef DUMUX_LOWREKEPSILON_NC_MODEL_HH
#define DUMUX_LOWREKEPSILON_NC_MODEL_HH

#include <dumux/common/properties.hh>
#include <dumux/freeflow/compositional/navierstokesncmodel.hh>
#include <dumux/freeflow/nonisothermal/vtkoutputfields.hh>
#include <dumux/freeflow/rans/twoeq/lowrekepsilon/model.hh>

#include "vtkoutputfields.hh"

#include <dumux/freeflow/ransnc/volumevariables.hh>

namespace Dumux {

///////////////////////////////////////////////////////////////////////////
// properties for the single-phase, multi-component low-Re k-epsilon model
///////////////////////////////////////////////////////////////////////////
namespace Properties {

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tags for the single-phase, multi-component isothermal low-Re k-epsilon model
NEW_TYPE_TAG(LowReKEpsilonNC, INHERITS_FROM(NavierStokesNC));

///////////////////////////////////////////////////////////////////////////
// default property values
///////////////////////////////////////////////////////////////////////////

/*!
 * \ingroup FreeflowNCModel
 * \brief Traits for the low-Reynolds k-epsilon multi-component model
 */
template<int dimension, int nComp, int phaseIdx, int replaceCompEqIdx, bool useMoles>
struct LowReKEpsilonNCModelTraits : NavierStokesNCModelTraits<dimension, nComp, phaseIdx, replaceCompEqIdx, useMoles>
{
    //! There are as many momentum balance equations as dimensions
    //! and as many balance equations as components.
    static constexpr int numEq() { return dimension+nComp+2; }

    //! The model does include a turbulence model
    static constexpr bool usesTurbulenceModel() { return true; }

    //! the indices
    using Indices = FreeflowNCIndices<dimension, numEq(), phaseIdx, replaceCompEqIdx, LowReKEpsilonIndices<dimension, nComp>>;
};

//!< states some specifics of the isothermal multi-component low-Reynolds k-epsilon model
SET_PROP(LowReKEpsilonNC, ModelTraits)
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
    using type = LowReKEpsilonNCModelTraits<dimension, numComponents, phaseIdx, replaceCompEqIdx, useMoles>;
};

//! Set the volume variables property
SET_PROP(LowReKEpsilonNC, VolumeVariables)
{
private:
    using PV = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using FSY = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using FST = typename GET_PROP_TYPE(TypeTag, FluidState);
    using MT = typename GET_PROP_TYPE(TypeTag, ModelTraits);

    using Traits = NavierStokesVolumeVariablesTraits<PV, FSY, FST, MT>;
    using NCVolVars = FreeflowNCVolumeVariables<Traits>;
    using RANSVolVars = LowReKEpsilonVolumeVariables<Traits, NCVolVars>;
public:
    using type = RANSNCVolumeVariables<Traits, RANSVolVars>;
};

//! The local residual
SET_PROP(LowReKEpsilonNC, LocalResidual)
{
private:
    using BaseLocalResidual = LowReKEpsilonResidual<TypeTag>;
public:
    using type = FreeflowNCResidual<TypeTag, BaseLocalResidual>;
};

//! The flux variables
SET_PROP(LowReKEpsilonNC, FluxVariables)
{
private:
    using BaseFluxVariables = LowReKEpsilonFluxVariables<TypeTag>;
public:
    using type = FreeflowNCFluxVariables<TypeTag, BaseFluxVariables>;
};

//! The specific vtk output fields
SET_PROP(LowReKEpsilonNC, VtkOutputFields)
{
private:
    using ModelTraits = typename GET_PROP_TYPE(TypeTag, ModelTraits);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    static constexpr int phaseIdx = GET_PROP_VALUE(TypeTag, PhaseIdx);
    using SinglePhaseVtkOutputFields = LowReKEpsilonVtkOutputFields<FVGridGeometry>;
public:
    using type = FreeflowNCVtkOutputFields<SinglePhaseVtkOutputFields, ModelTraits, FVGridGeometry, FluidSystem, phaseIdx>;
};

//////////////////////////////////////////////////////////////////////////
// Property values for non-isothermal multi-component low-Re k-epsilon model
//////////////////////////////////////////////////////////////////////////

//! The type tags for the single-phase, multi-component non-isothermal low-Re k-epsilon models
NEW_TYPE_TAG(LowReKEpsilonNCNI, INHERITS_FROM(NavierStokesNCNI));

//! The model traits of the non-isothermal model
SET_PROP(LowReKEpsilonNCNI, ModelTraits)
{
private:
    using GridView = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::GridView;
    static constexpr int dim = GridView::dimension;
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    static constexpr int numComponents = FluidSystem::numComponents;
    static constexpr int phaseIdx = GET_PROP_VALUE(TypeTag, PhaseIdx);
    static constexpr int replaceCompEqIdx = GET_PROP_VALUE(TypeTag, ReplaceCompEqIdx);
    static constexpr bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);
    using IsothermalModelTraits = LowReKEpsilonNCModelTraits<dim, numComponents, phaseIdx, replaceCompEqIdx, useMoles>;
public:
    using type = FreeflowNIModelTraits<IsothermalModelTraits>;
};

//! Set the volume variables property
SET_PROP(LowReKEpsilonNCNI, VolumeVariables)
{
private:
    using PV = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using FSY = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using FST = typename GET_PROP_TYPE(TypeTag, FluidState);
    using MT = typename GET_PROP_TYPE(TypeTag, ModelTraits);

    using Traits = NavierStokesVolumeVariablesTraits<PV, FSY, FST, MT>;
    using NCVolVars = FreeflowNCVolumeVariables<Traits>;
    using RANSVolVars = LowReKEpsilonVolumeVariables<Traits, NCVolVars>;
public:
    using type = RANSNCVolumeVariables<Traits, RANSVolVars>;
};

//! The local residual
SET_PROP(LowReKEpsilonNCNI, LocalResidual)
{
private:
    using BaseLocalResidual = LowReKEpsilonResidual<TypeTag>;
public:
    using type = FreeflowNCResidual<TypeTag, BaseLocalResidual>;
};

//! The flux variables
SET_PROP(LowReKEpsilonNCNI, FluxVariables)
{
private:
    using BaseFluxVariables = LowReKEpsilonFluxVariables<TypeTag>;
public:
    using type = FreeflowNCFluxVariables<TypeTag, BaseFluxVariables>;
};

//! The specific vtk output fields
SET_PROP(LowReKEpsilonNCNI, VtkOutputFields)
{
private:
    using ModelTraits = typename GET_PROP_TYPE(TypeTag, ModelTraits);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    static constexpr int phaseIdx = GET_PROP_VALUE(TypeTag, PhaseIdx);
    using BaseVtkOutputFields = LowReKEpsilonVtkOutputFields<FVGridGeometry>;
    using NonIsothermalFields = FreeflowNonIsothermalVtkOutputFields<BaseVtkOutputFields, ModelTraits>;
public:
    using type = FreeflowNCVtkOutputFields<NonIsothermalFields, ModelTraits, FVGridGeometry, FluidSystem, phaseIdx>;
};

// \}
} // end namespace Properties
} // end namespace Dumux

#endif
