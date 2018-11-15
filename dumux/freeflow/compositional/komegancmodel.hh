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
#include <dumux/freeflow/nonisothermal/iofields.hh>
#include <dumux/freeflow/rans/twoeq/komega/model.hh>

#include "iofields.hh"

namespace Dumux {

///////////////////////////////////////////////////////////////////////////
// properties for the single-phase, multi-component k-omega model
///////////////////////////////////////////////////////////////////////////
namespace Properties {

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

// Create new type tags
namespace TTag {
//! The type tags for the single-phase, multi-component isothermal k-omega model
struct KOmegaNC { using InheritsFrom = std::tuple<NavierStokesNC>; };
} // end namespace TTag

///////////////////////////////////////////////////////////////////////////
// default property values
///////////////////////////////////////////////////////////////////////////

/*!
 * \ingroup FreeflowNCModel
 * \brief Traits for the low-Reynolds k-epsilon multi-component model
 *
 * \tparam dimension The dimension of the problem
 * \tparam nComp The number of components to be considered
 * \tparam useM Use molar or mass balances
 * \tparam replaceCompEqIdx The index of the component balance equation that should be replaced by a total mass/mole balance
 */
template<int dimension, int nComp, bool useMoles, int replaceCompEqIdx>
struct KOmegaNCModelTraits : NavierStokesNCModelTraits<dimension, nComp, useMoles, replaceCompEqIdx>
{
    //! There are as many momentum balance equations as dimensions
    //! and as many balance equations as components.
    static constexpr int numEq() { return dimension+nComp+2; }

    //! The model does include a turbulence model
    static constexpr bool usesTurbulenceModel() { return true; }

    //! the indices
    using Indices = KOmegaIndices<dimension, nComp>;

    //! return the names of the primary variables in cells
    template <class FluidSystem>
    static std::string primaryVariableNameCell(int pvIdx, int state = 0)
    {
        using ParentType = NavierStokesNCModelTraits<dimension, nComp, useMoles, replaceCompEqIdx>;
        if (pvIdx < nComp)
            return ParentType::template primaryVariableNameCell<FluidSystem>(pvIdx, state);
        else if (pvIdx == nComp)
            return "k";
        else
            return "omega";
    }
};

//!< states some specifics of the isothermal multi-component low-Reynolds k-epsilon model
SET_PROP(KOmegaNC, ModelTraits)
{
private:
    using GridView = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::GridView;
    static constexpr int dimension = GridView::dimension;
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    static constexpr int numComponents = FluidSystem::numComponents;
    static constexpr bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);
    static constexpr int replaceCompEqIdx = GET_PROP_VALUE(TypeTag, ReplaceCompEqIdx);
public:
    using type = KOmegaNCModelTraits<dimension, numComponents, useMoles, replaceCompEqIdx>;
};

//! Set the volume variables property
SET_PROP(KOmegaNC, VolumeVariables)
{
private:
    using PV = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using FSY = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using FST = typename GET_PROP_TYPE(TypeTag, FluidState);
    using MT = typename GET_PROP_TYPE(TypeTag, ModelTraits);

    static_assert(FSY::numComponents == MT::numComponents(), "Number of components mismatch between model and fluid system");
    static_assert(FST::numComponents == MT::numComponents(), "Number of components mismatch between model and fluid state");
    static_assert(FSY::numPhases == MT::numPhases(), "Number of phases mismatch between model and fluid system");
    static_assert(FST::numPhases == MT::numPhases(), "Number of phases mismatch between model and fluid state");

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

//! The specific I/O fields
SET_TYPE_PROP(KOmegaNC, IOFields, FreeflowNCIOFields<KOmegaIOFields, true/*turbulenceModel*/>);

//////////////////////////////////////////////////////////////////////////
// Property values for non-isothermal multi-component k-omega model
//////////////////////////////////////////////////////////////////////////

// Create new type tags
namespace TTag {
//! The type tags for the single-phase, multi-component non-isothermal k-omega models
struct KOmegaNCNI { using InheritsFrom = std::tuple<NavierStokesNCNI>; };
} // end namespace TTag

//! The model traits of the non-isothermal model
SET_PROP(KOmegaNCNI, ModelTraits)
{
private:
    using GridView = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::GridView;
    static constexpr int dim = GridView::dimension;
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    static constexpr int numComponents = FluidSystem::numComponents;
    static constexpr bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);
    static constexpr int replaceCompEqIdx = GET_PROP_VALUE(TypeTag, ReplaceCompEqIdx);
    using IsothermalModelTraits = KOmegaNCModelTraits<dim, numComponents, useMoles, replaceCompEqIdx>;
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

    static_assert(FSY::numComponents == MT::numComponents(), "Number of components mismatch between model and fluid system");
    static_assert(FST::numComponents == MT::numComponents(), "Number of components mismatch between model and fluid state");
    static_assert(FSY::numPhases == MT::numPhases(), "Number of phases mismatch between model and fluid system");
    static_assert(FST::numPhases == MT::numPhases(), "Number of phases mismatch between model and fluid state");

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

//! The specific I/O fields
SET_PROP(KOmegaNCNI, IOFields)
{
private:
    using IsothermalIOFields = FreeflowNCIOFields<KOmegaIOFields, true/*turbulenceModel*/>;
public:
    using type = FreeflowNonIsothermalIOFields<IsothermalIOFields, true/*turbulenceModel*/>;
};

// \}
} // end namespace Properties
} // end namespace Dumux

#endif
