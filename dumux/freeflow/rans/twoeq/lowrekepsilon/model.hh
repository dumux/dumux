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
 * \ingroup LowReKEpsilonModel
 *
 * \brief A single-phase, isothermal low-Reynolds k-epsilon model
 *
 * \copydoc RANSModel
 *
 * These models calculate the eddy viscosity with two additional PDEs,
 * one for the turbulent kinetic energy (k) and for the dissipation (epsilon).
 */

#ifndef DUMUX_LOWREKEPSILON_MODEL_HH
#define DUMUX_LOWREKEPSILON_MODEL_HH

#include <dumux/common/properties.hh>
#include <dumux/freeflow/properties.hh>
#include <dumux/freeflow/rans/model.hh>

#include "fluxvariables.hh"
#include "indices.hh"
#include "localresidual.hh"
#include "volumevariables.hh"
#include "vtkoutputfields.hh"

namespace Dumux
{
namespace Properties {

/*!
 * \ingroup LowReKEpsilonModel
 * \brief Traits for the low-Reynolds k-epsilon model
 */
template<int dimension>
struct LowReKEpsilonModelTraits
{
    //! The dimension of the model
    static constexpr int dim() { return dimension; }

    //! There are as many momentum balance equations as dimensions,
    //! one mass balance equation and two turbulent transport equations
    static constexpr int numEq() { return dim()+1+2; }

    //! The number of phases is always 1
    static constexpr int numPhases() { return 1; }

    //! The number of components
    static constexpr int numComponents() { return 1; }

    //! Enable advection
    static constexpr bool enableAdvection() { return true; }

    //! The one-phase model has no molecular diffusion
    static constexpr bool enableMolecularDiffusion() { return true; }

    //! The model is isothermal
    static constexpr bool enableEnergyBalance() { return false; }

    //! the indices
    using Indices = LowReKEpsilonIndices<dim()>;
};

///////////////////////////////////////////////////////////////////////////
// default property values for the isothermal low-Reynolds k-epsilon model
///////////////////////////////////////////////////////////////////////////

//! The type tag for the single-phase, isothermal low-Reynolds k-epsilon model
NEW_TYPE_TAG(LowReKEpsilon, INHERITS_FROM(RANS));

//!< states some specifics of the isothermal low-Reynolds k-epsilon model
SET_PROP(LowReKEpsilon, ModelTraits)
{
private:
    using GridView = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::GridView;
    static constexpr int dim = GridView::dimension;
public:
    using type = LowReKEpsilonModelTraits<dim>;
};

//! The flux variables
SET_TYPE_PROP(LowReKEpsilon, FluxVariables, LowReKEpsilonFluxVariables<TypeTag>);

//! The local residual
SET_TYPE_PROP(LowReKEpsilon, LocalResidual, LowReKEpsilonResidual<TypeTag>);

//! Set the volume variables property
SET_PROP(LowReKEpsilon, VolumeVariables)
{
private:
    using PV = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using FSY = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using FST = typename GET_PROP_TYPE(TypeTag, FluidState);
    using MT = typename GET_PROP_TYPE(TypeTag, ModelTraits);

    using Traits = NavierStokesVolumeVariablesTraits<PV, FSY, FST, MT>;
    using NSVolVars = NavierStokesVolumeVariables<Traits>;
public:
    using type = LowReKEpsilonVolumeVariables<Traits, NSVolVars>;
};

//! The specific vtk output fields
SET_PROP(LowReKEpsilon, VtkOutputFields)
{
private:
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
public:
    using type = LowReKEpsilonVtkOutputFields<FVGridGeometry>;
};

//////////////////////////////////////////////////////////////////
// default property values for the non-isothermal low-Reynolds k-epsilon model
//////////////////////////////////////////////////////////////////

//! The type tag for the single-phase, isothermal low-Reynolds k-epsilon model
NEW_TYPE_TAG(LowReKEpsilonNI, INHERITS_FROM(RANSNI));

//! Set the volume variables property
SET_PROP(LowReKEpsilonNI, VolumeVariables)
{
private:
    using PV = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using FSY = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using FST = typename GET_PROP_TYPE(TypeTag, FluidState);
    using MT = typename GET_PROP_TYPE(TypeTag, ModelTraits);

    using Traits = NavierStokesVolumeVariablesTraits<PV, FSY, FST, MT>;
    using NSVolVars = NavierStokesVolumeVariables<Traits>;
public:
    using type = LowReKEpsilonVolumeVariables<Traits, NSVolVars>;
};

// \}
}

} // end namespace

#endif // DUMUX_LOWREKEPSILON_MODEL_HH
