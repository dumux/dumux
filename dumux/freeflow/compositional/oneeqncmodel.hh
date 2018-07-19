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
 * \brief A single-phase, multi-component one-equation model
 *
 * \copydoc Dumux::FreeflowNCModel
 */

#ifndef DUMUX_ONEEQ_NC_MODEL_HH
#define DUMUX_ONEEQ_NC_MODEL_HH

#include <dumux/common/properties.hh>
#include <dumux/freeflow/compositional/navierstokesncmodel.hh>
#include <dumux/freeflow/nonisothermal/vtkoutputfields.hh>
#include <dumux/freeflow/rans/oneeq/model.hh>

#include "vtkoutputfields.hh"

namespace Dumux {

///////////////////////////////////////////////////////////////////////////
// properties for the single-phase, multi-component one-equation model
///////////////////////////////////////////////////////////////////////////
namespace Properties {

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tags for the single-phase, multi-component isothermal one-equation model
NEW_TYPE_TAG(OneEqNC, INHERITS_FROM(NavierStokesNC));

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

    //! the indices
    using Indices = OneEqIndices<dimension, nComp>;
};

//!< states some specifics of the isothermal multi-component one-equation model
SET_PROP(OneEqNC, ModelTraits)
{
private:
    using GridView = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::GridView;
    static constexpr int dimension = GridView::dimension;
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    static constexpr int numComponents = FluidSystem::numComponents;
    static constexpr bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);
    static constexpr int replaceCompEqIdx = GET_PROP_VALUE(TypeTag, ReplaceCompEqIdx);
public:
    using type = OneEqNCModelTraits<dimension, numComponents, useMoles, replaceCompEqIdx>;
};

//! Set the volume variables property
SET_PROP(OneEqNC, VolumeVariables)
{
private:
    using PV = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using FSY = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using FST = typename GET_PROP_TYPE(TypeTag, FluidState);
    using MT = typename GET_PROP_TYPE(TypeTag, ModelTraits);

    using Traits = NavierStokesVolumeVariablesTraits<PV, FSY, FST, MT>;
    using NCVolVars = FreeflowNCVolumeVariables<Traits>;
public:
    using type = OneEqVolumeVariables<Traits, NCVolVars>;
};

//! The local residual
SET_PROP(OneEqNC, LocalResidual)
{
private:
    using BaseLocalResidual = FreeflowNCResidual<TypeTag>;
public:
    using type = OneEqResidual<TypeTag, BaseLocalResidual>;
};

//! The flux variables
SET_PROP(OneEqNC, FluxVariables)
{
private:
    using BaseFluxVariables = FreeflowNCFluxVariables<TypeTag>;
public:
    using type = OneEqFluxVariables<TypeTag, BaseFluxVariables>;
};

//! The specific vtk output fields
SET_PROP(OneEqNC, VtkOutputFields)
{
private:
    using ModelTraits = typename GET_PROP_TYPE(TypeTag, ModelTraits);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using SinglePhaseVtkOutputFields = OneEqVtkOutputFields<FVGridGeometry>;
public:
    using type = FreeflowNCVtkOutputFields<SinglePhaseVtkOutputFields, ModelTraits, FVGridGeometry, FluidSystem>;
};

//////////////////////////////////////////////////////////////////////////
// Property values for non-isothermal multi-component one-equation model
//////////////////////////////////////////////////////////////////////////

//! The type tags for the single-phase, multi-component non-isothermal one-equation models
NEW_TYPE_TAG(OneEqNCNI, INHERITS_FROM(NavierStokesNCNI));

//! The model traits of the non-isothermal model
SET_PROP(OneEqNCNI, ModelTraits)
{
private:
    using GridView = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::GridView;
    static constexpr int dim = GridView::dimension;
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    static constexpr int numComponents = FluidSystem::numComponents;
    static constexpr bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);
    static constexpr int replaceCompEqIdx = GET_PROP_VALUE(TypeTag, ReplaceCompEqIdx);
    using IsothermalModelTraits = OneEqNCModelTraits<dim, numComponents, useMoles, replaceCompEqIdx>;
public:
    using type = FreeflowNIModelTraits<IsothermalModelTraits>;
};

//! Set the volume variables property
SET_PROP(OneEqNCNI, VolumeVariables)
{
private:
    using PV = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using FSY = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using FST = typename GET_PROP_TYPE(TypeTag, FluidState);
    using MT = typename GET_PROP_TYPE(TypeTag, ModelTraits);

    using Traits = NavierStokesVolumeVariablesTraits<PV, FSY, FST, MT>;
    using NCVolVars = FreeflowNCVolumeVariables<Traits>;
public:
    using type = OneEqVolumeVariables<Traits, NCVolVars>;
};

//! The local residual
SET_PROP(OneEqNCNI, LocalResidual)
{
private:
    using BaseLocalResidual = FreeflowNCResidual<TypeTag>;
public:
    using type = OneEqResidual<TypeTag, BaseLocalResidual>;
};

//! The flux variables
SET_PROP(OneEqNCNI, FluxVariables)
{
private:
    using BaseFluxVariables = FreeflowNCFluxVariables<TypeTag>;
public:
    using type = OneEqFluxVariables<TypeTag, BaseFluxVariables>;
};

//! The specific vtk output fields
SET_PROP(OneEqNCNI, VtkOutputFields)
{
private:
    using ModelTraits = typename GET_PROP_TYPE(TypeTag, ModelTraits);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using BaseVtkOutputFields = OneEqVtkOutputFields<FVGridGeometry>;
    using NonIsothermalFields = FreeflowNonIsothermalVtkOutputFields<BaseVtkOutputFields, ModelTraits>;
public:
    using type = FreeflowNCVtkOutputFields<NonIsothermalFields, ModelTraits, FVGridGeometry, FluidSystem>;
};

// \}
} // end namespace Properties
} // end namespace Dumux

#endif
