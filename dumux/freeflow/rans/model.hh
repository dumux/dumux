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
 * \brief A single-phase, isothermal Reynolds-Averaged Navier-Stokes model
 *
 * This model implements a single-phase, isothermal Reynolds-Averaged
 * Navier-Stokes model, solving the <B> momentum balance equation </B>
 * \f[
 \frac{\partial (\varrho \textbf{v})}{\partial t} + \nabla \cdot (\varrho \textbf{v} \textbf{v}^{\textup{T}}) = \nabla \cdot (\mu_\textrm{eff} (\nabla \textbf{v} + \nabla \textbf{v}^{\textup{T}}))
   - \nabla p + \varrho \textbf{g} - \textbf{f}
 * \f]
 * The effective viscosity is composed of the fluid and the eddy viscosity:
 * \f[
 *    \mu_\textrm{eff} = \mu + \mu_\textrm{t}
 * \f].
 */

#ifndef DUMUX_RANS_MODEL_HH
#define DUMUX_RANS_MODEL_HH

#include <dumux/common/properties.hh>
#include <dumux/discretization/methods.hh>
#include <dumux/freeflow/properties.hh>
#include <dumux/freeflow/navierstokes/model.hh>
#include <dumux/freeflow/navierstokes/indices.hh>
#include <dumux/freeflow/nonisothermal/indices.hh>
#include <dumux/material/fluidstates/immiscible.hh>

#include "volumevariables.hh"
#include "vtkoutputfields.hh"

namespace Dumux
{

// \{
///////////////////////////////////////////////////////////////////////////
// properties for the single-phase Reynolds-Averaged Navier-Stokes model
///////////////////////////////////////////////////////////////////////////
namespace Properties {

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tag for the single-phase, isothermal Reynolds-Averaged Navier-Stokes model
NEW_TYPE_TAG(RANS, INHERITS_FROM(NavierStokes));

///////////////////////////////////////////////////////////////////////////
// default property values for the isothermal single phase model
///////////////////////////////////////////////////////////////////////////
SET_BOOL_PROP(RANS, EnableInertiaTerms, true); //!< Explicitly force the consideration of inertia terms by default

//! The volume variables
SET_TYPE_PROP(RANS, VolumeVariables, RANSVolumeVariables<TypeTag>);

//! The specific vtk output fields
SET_PROP(RANS, VtkOutputFields)
{
private:
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
public:
     using type = RANSVtkOutputFields<FVGridGeometry>;
};

//////////////////////////////////////////////////////////////////
// Property values for non-isothermal Reynolds-averaged Navier-Stokes model
//////////////////////////////////////////////////////////////////

//! The type tag for the single-phase, isothermal Reynolds-Averaged Navier-Stokes model
NEW_TYPE_TAG(RANSNI, INHERITS_FROM(RANS));

//! The model traits of the non-isothermal model
SET_PROP(RANSNI, ModelTraits)
{
private:
    using GridView = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::GridView;
    static constexpr int dim = GridView::dimension;
    using IsothermalTraits = NavierStokesModelTraits<dim>;
public:
    using type = NavierStokesNIModelTraits<IsothermalTraits>;
};

//! The indices required by the non-isothermal single-phase model
SET_PROP(RANSNI, Indices)
{
private:
    static constexpr int numEq = GET_PROP_TYPE(TypeTag, ModelTraits)::numEq();
    static constexpr int dim = GET_PROP_TYPE(TypeTag, GridView)::dimension;
    using IsothermalIndices = NavierStokesIndices<dim, numEq>;
public:
    using type = NavierStokesNonIsothermalIndices<dim, numEq, IsothermalIndices>;
};

//! The specific non-isothermal vtk output fields
// SET_PROP(RANSNI, VtkOutputFields)
// {
// private:
//      using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
//      using IsothermalFields = NavierStokesVtkOutputFields<FVGridGeometry>;
// public:
//      using type = NavierStokesNonIsothermalVtkOutputFields<IsothermalFields>;
// };

//! Use Fourier's Law as default heat conduction type
SET_TYPE_PROP(RANSNI, HeatConductionType, FouriersLaw<TypeTag>);

// \}
}

} // end namespace

#endif // DUMUX_RANS_MODEL_HH
