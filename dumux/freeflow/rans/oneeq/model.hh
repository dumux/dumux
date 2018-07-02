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
 * \ingroup OneEqModel
 *
 * \brief A single-phase, isothermal one-equation turbulence model by Spalart-Allmaras
 *
 * \copydoc RANSModel
 *
 * \todo please doc me
 */

#ifndef DUMUX_ONEEQ_MODEL_HH
#define DUMUX_ONEEQ_MODEL_HH

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
 * \ingroup OneEqModel
 * \brief Traits for the Spalart-Allmaras model
 */
template<int dimension>
struct OneEqModelTraits : RANSModelTraits<dimension>
{
    //! The dimension of the model
    static constexpr int dim() { return dimension; }

    //! There are as many momentum balance equations as dimensions,
    //! one mass balance equation and two turbulent transport equations
    static constexpr int numEq() { return dim()+1+1; }

    //! The number of components
    static constexpr int numComponents() { return 1; }

    //! the indices
    using Indices = OneEqIndices<dim(), numComponents()>;
};

///////////////////////////////////////////////////////////////////////////
// default property values for the isothermal Spalart-Allmaras model
///////////////////////////////////////////////////////////////////////////

//! The type tag for the single-phase, isothermal Spalart-Allmaras model
NEW_TYPE_TAG(OneEq, INHERITS_FROM(RANS));

//!< states some specifics of the isothermal Spalart-Allmaras model
SET_PROP(OneEq, ModelTraits)
{
private:
    using GridView = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::GridView;
    static constexpr int dim = GridView::dimension;
public:
    using type = OneEqModelTraits<dim>;
};

//! The flux variables
SET_PROP(OneEq, FluxVariables)
{
private:
    using BaseFluxVariables = NavierStokesFluxVariables<TypeTag>;
public:
    using type = OneEqFluxVariables<TypeTag, BaseFluxVariables>;
};

//! The local residual
SET_PROP(OneEq, LocalResidual)
{
private:
    using BaseLocalResidual = NavierStokesResidual<TypeTag>;
public:
    using type = OneEqResidual<TypeTag, BaseLocalResidual>;
};

//! Set the volume variables property
SET_PROP(OneEq, VolumeVariables)
{
private:
    using PV = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using FSY = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using FST = typename GET_PROP_TYPE(TypeTag, FluidState);
    using MT = typename GET_PROP_TYPE(TypeTag, ModelTraits);

    using Traits = NavierStokesVolumeVariablesTraits<PV, FSY, FST, MT>;
    using NSVolVars = NavierStokesVolumeVariables<Traits>;
public:
    using type = OneEqVolumeVariables<Traits, NSVolVars>;
};

//! The specific vtk output fields
SET_PROP(OneEq, VtkOutputFields)
{
private:
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
public:
    using type = OneEqVtkOutputFields<FVGridGeometry>;
};

//////////////////////////////////////////////////////////////////
// default property values for the non-isothermal Spalart-Allmaras model
//////////////////////////////////////////////////////////////////

//! The type tag for the single-phase, isothermal Spalart-Allmaras model
NEW_TYPE_TAG(OneEqNI, INHERITS_FROM(RANSNI));

//! The model traits of the non-isothermal model
SET_PROP(OneEqNI, ModelTraits)
{
private:
    using GridView = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::GridView;
    static constexpr int dim = GridView::dimension;
    using IsothermalTraits = OneEqModelTraits<dim>;
public:
    using type = FreeflowNIModelTraits<IsothermalTraits>;
};

//! Set the volume variables property
SET_PROP(OneEqNI, VolumeVariables)
{
private:
    using PV = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using FSY = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using FST = typename GET_PROP_TYPE(TypeTag, FluidState);
    using MT = typename GET_PROP_TYPE(TypeTag, ModelTraits);

    using Traits = NavierStokesVolumeVariablesTraits<PV, FSY, FST, MT>;
    using NSVolVars = NavierStokesVolumeVariables<Traits>;
public:
    using type = OneEqVolumeVariables<Traits, NSVolVars>;
};

//! The specific non-isothermal vtk output fields
SET_PROP(OneEqNI, VtkOutputFields)
{
private:
    using ModelTraits = typename GET_PROP_TYPE(TypeTag, ModelTraits);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using IsothermalFields = OneEqVtkOutputFields<FVGridGeometry>;
public:
    using type = FreeflowNonIsothermalVtkOutputFields<IsothermalFields, ModelTraits>;
};

// \}
}

} // end namespace

#endif
