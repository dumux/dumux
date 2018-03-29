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
 * \ingroup RANSModel
 *
 * \brief A single-phase, multi-component isothermal Navier-Stokes model
 *
 * \copydoc Dumux::RANSModel
 *
 * The system is closed by a <B> component mass/mole balance equation </B> for each component \f$\kappa\f$:
 * \f[
 *    \frac{\partial \left(\varrho X^\kappa\right)}{\partial t}
 *    + \nabla \cdot \left( \varrho {\boldsymbol{v}} X^\kappa
 *    - (D^\kappa + D_\text{t}) \varrho \frac{M^\kappa}{M} \textbf{grad}\, x^\kappa \right)
 *    - q^\kappa = 0
 * \f]
 *
 * Alternatively, one component balance equation can be replace by a <B> total mass/mole balance equation </B>:
 * \f[
 *    \frac{\partial \varrho_g}{\partial t}
 *    + \nabla \cdot \left(
 *        \varrho {\boldsymbol{v}}
 *        - \sum_\kappa (D^\kappa + D_\text{t}) \varrho \frac{M^\kappa}{M} \textbf{grad}\, x^\kappa
 *      \right)
 *    - q = 0
 * \f]
 *
 * The eddy diffusivity \f$ D_\text{t} \f$ is related to the eddy viscosity \f$ \nu_\text{t} \f$
 * by the turbulent Schmidt number:
 * \f[ D_\text{t} = \frac{\nu_\text{t}}{\mathrm{Sc}_\text{t}} \f]
 *
 * So far, only the staggered grid spatial discretization (for structured grids) is available.
 */

#ifndef DUMUX_RANS_NC_MODEL_HH
#define DUMUX_RANS_NC_MODEL_HH

#include <dumux/common/properties.hh>

#include <dumux/freeflow/rans/model.hh>
#include <dumux/freeflow/rans/zeroeq/model.hh>
#include <dumux/freeflow/navierstokesnc/model.hh>
#include <dumux/freeflow/nonisothermal/model.hh>
#include <dumux/freeflow/nonisothermal/indices.hh>
#include <dumux/freeflow/nonisothermal/vtkoutputfields.hh>

#include "indices.hh"
#include "volumevariables.hh"
#include "vtkoutputfields.hh"

namespace Dumux {

/*!
 * \ingroup RANSModel
 * \brief Traits for the RANS multi-component model
 */
template<int dim, int nComp>
struct RANSNCModelTraits : NavierStokesNCModelTraits<dim, nComp>
{ };

///////////////////////////////////////////////////////////////////////////
// properties for the single-phase, multi-component RANS model
///////////////////////////////////////////////////////////////////////////
namespace Properties {

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tags for the single-phase, multi-component isothermal RANS model
NEW_TYPE_TAG(RANSNC, INHERITS_FROM(RANS, NavierStokesNC));
NEW_TYPE_TAG(ZeroEqNC, INHERITS_FROM(ZeroEq, RANSNC));

///////////////////////////////////////////////////////////////////////////
// default property values
///////////////////////////////////////////////////////////////////////////

//! The volume variables
SET_TYPE_PROP(RANSNC, VolumeVariables, RANSNCVolumeVariables<TypeTag>);
SET_TYPE_PROP(ZeroEqNC, SinglePhaseVolumeVariables, ZeroEqVolumeVariables<TypeTag>);

//! The indices
SET_PROP(RANSNC, Indices)
{
private:
    static constexpr int numEq = GET_PROP_TYPE(TypeTag, ModelTraits)::numEq();
    static constexpr int dim = GET_PROP_TYPE(TypeTag, GridView)::dimension;
    static constexpr int phaseIdx = GET_PROP_VALUE(TypeTag, PhaseIdx);
    static constexpr int replaceCompEqIdx = GET_PROP_VALUE(TypeTag, ReplaceCompEqIdx);
public:
    using type = RANSNCIndices<dim, numEq, phaseIdx, replaceCompEqIdx>;
};

//! The specific vtk output fields
SET_PROP(ZeroEqNC, VtkOutputFields)
{
private:
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    static constexpr int phaseIdx = GET_PROP_VALUE(TypeTag, PhaseIdx);
    using SinglePhaseVtkOutputFields = RANSVtkOutputFields<FVGridGeometry>;
public:
    using type = RANSNCVtkOutputFields<FVGridGeometry, FluidSystem, phaseIdx, SinglePhaseVtkOutputFields>;
};

//////////////////////////////////////////////////////////////////////////
// Property values for non-isothermal multi-component RANS model
//////////////////////////////////////////////////////////////////////////

//! The type tags for the single-phase, multi-component non-isothermal RANS models
NEW_TYPE_TAG(RANSNCNI, INHERITS_FROM(RANSNI, NavierStokesNC));
NEW_TYPE_TAG(ZeroEqNCNI, INHERITS_FROM(ZeroEqNI, RANSNCNI));

//! The volume variables
SET_TYPE_PROP(RANSNCNI, VolumeVariables, RANSNCVolumeVariables<TypeTag>);
SET_TYPE_PROP(ZeroEqNCNI, SinglePhaseVolumeVariables, ZeroEqVolumeVariables<TypeTag>);

//! The model traits of the non-isothermal model
SET_PROP(RANSNCNI, ModelTraits)
{
private:
    using GridView = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::GridView;
    static constexpr int dim = GridView::dimension;
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    static constexpr int numComponents = FluidSystem::numComponents;
    using IsothermalModelTraits = NavierStokesNCModelTraits<dim, numComponents>;
public:
    using type = NavierStokesNIModelTraits<IsothermalModelTraits>;
};

//! The indices
SET_PROP(RANSNCNI, Indices)
{
private:
    static constexpr int numEq = GET_PROP_TYPE(TypeTag, ModelTraits)::numEq();
    static constexpr int dim = GET_PROP_TYPE(TypeTag, GridView)::dimension;
    static constexpr int phaseIdx = GET_PROP_VALUE(TypeTag, PhaseIdx);
    static constexpr int replaceCompEqIdx = GET_PROP_VALUE(TypeTag, ReplaceCompEqIdx);
    using IsothermalIndices = NavierStokesNCIndices<dim, numEq, phaseIdx, replaceCompEqIdx>;
public:
    using type = NavierStokesNonIsothermalIndices<dim, numEq, IsothermalIndices>;
};

//! The specific vtk output fields
SET_PROP(ZeroEqNCNI, VtkOutputFields)
{
private:
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    static constexpr int phaseIdx = GET_PROP_VALUE(TypeTag, PhaseIdx);
    using NavierStokesFields = NavierStokesVtkOutputFields<FVGridGeometry>;
    using RANSFields = RANSVtkOutputFields<FVGridGeometry>;
    using SinglePhaseVtkOutputFields = RANSNonIsothermalVtkOutputFields<NavierStokesFields, RANSFields>;
public:
    using type = RANSNCVtkOutputFields<FVGridGeometry, FluidSystem, phaseIdx, SinglePhaseVtkOutputFields>;
};

// \}
} // end namespace Properties
} // end namespace Dumux


#endif
