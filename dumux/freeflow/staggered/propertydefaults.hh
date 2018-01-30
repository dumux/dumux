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
 * \ingroup Properties
 * \ingroup ImplicitProperties
 * \ingroup NavierStokesModel
 * \file
 *
 * \brief Defines the properties required for the one-phase fully implicit model.
 */
#ifndef DUMUX_NAVIERSTOKES_PROPERTY_DEFAULTS_HH
#define DUMUX_NAVIERSTOKES_PROPERTY_DEFAULTS_HH

#include "properties.hh"

#include "model.hh"
#include "volumevariables.hh"
#include "indices.hh"
#include "problem.hh"
#include "localresidual.hh"
#include "fluxvariables.hh"
#include "fluxvariablescache.hh"
#include "velocityoutput.hh"
#include "vtkoutputmodule.hh"
#include "boundarytypes.hh"

#include <dumux/implicit/staggered/localresidual.hh>
#include <dumux/freeflow/staggeredni/localresidual.hh>
#include <dumux/freeflow/staggeredni/fluxvariables.hh>

#include <dumux/material/fluidsystems/gasphase.hh>
#include <dumux/material/fluidsystems/liquidphase.hh>
#include <dumux/material/components/nullcomponent.hh>
#include <dumux/material/fluidsystems/1p.hh>

#include <dumux/discretization/staggered/freeflow/staggeredgeometryhelper.hh>
#include <dumux/discretization/staggered/freeflow/subcontrolvolumeface.hh>
#include <dumux/discretization/staggered/freeflow/facevariables.hh>



namespace Dumux
{

namespace Properties
{
// forward declaration
NEW_PROP_TAG(FluxVariables);
NEW_PROP_TAG(FluxVariablesCache);
NEW_PROP_TAG(StaggeredGeometryHelper);
}
// \{

///////////////////////////////////////////////////////////////////////////
// default property values for the isothermal single phase model
///////////////////////////////////////////////////////////////////////////
namespace Properties {
SET_INT_PROP(NavierStokes, NumEqCellCenter, 1); //!< set the number of equations to 1
SET_INT_PROP(NavierStokes, NumEqFace, 1); //!< set the number of equations to 1
SET_INT_PROP(NavierStokes, NumPhases, 1); //!< The number of phases in the 1p model is 1

//! The sub-controlvolume face
SET_PROP(NavierStokes, SubControlVolumeFace)
{
private:
    using Grid = typename GET_PROP_TYPE(TypeTag, Grid);
    using ScvfGeometry = typename Grid::template Codim<1>::Geometry;
    using IndexType = typename Grid::LeafGridView::IndexSet::IndexType;
public:
    typedef Dumux::StaggeredSubControlVolumeFace<ScvfGeometry, IndexType> type;
};

//! The geometry helper required for the stencils, etc.
SET_PROP(NavierStokes, StaggeredGeometryHelper)
{
private:
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
public:
    using type = StaggeredGeometryHelper<GridView>;
};

//! The variables living on the faces
SET_TYPE_PROP(NavierStokes, FaceVariables, StaggeredFaceVariables<TypeTag>);

//! The local residual function
SET_TYPE_PROP(NavierStokes, LocalResidual, StaggeredNavierStokesResidual<TypeTag>);

//! the Model property
SET_TYPE_PROP(NavierStokes, Model, NavierStokesModel<TypeTag>);

//! the VolumeVariables property
SET_TYPE_PROP(NavierStokes, VolumeVariables, NavierStokesVolumeVariables<TypeTag>);

//! The class that contains the different flux variables (i.e. darcy, diffusion, energy)
//! by default, we set the flux variables to ones for porous media
SET_TYPE_PROP(NavierStokes, FluxVariables, FreeFlowFluxVariables<TypeTag>);

//! The flux variables cache class, by default the one for porous media
SET_TYPE_PROP(NavierStokes, FluxVariablesCache, FreeFlowFluxVariablesCache<TypeTag>);

//! Enable advection
SET_BOOL_PROP(NavierStokes, EnableAdvection, true);

//! The one-phase model has no molecular diffusion
SET_BOOL_PROP(NavierStokes, EnableMolecularDiffusion, false);

//! The indices required by the isothermal single-phase model
SET_TYPE_PROP(NavierStokes, Indices, NavierStokesCommonIndices<TypeTag>);

//! The weight of the upwind control volume when calculating
//! fluxes. Use central differences by default.
SET_SCALAR_PROP(NavierStokes, ImplicitMassUpwindWeight, 0.5);

//! weight for the upwind mobility in the velocity calculation
//! fluxes. Use central differences by default.
SET_SCALAR_PROP(NavierStokes, ImplicitMobilityUpwindWeight, 0.5);

//! The fluid system to use by default
SET_TYPE_PROP(NavierStokes, FluidSystem, Dumux::FluidSystems::OneP<typename GET_PROP_TYPE(TypeTag, Scalar), typename GET_PROP_TYPE(TypeTag, Fluid)>);

SET_PROP(NavierStokes, Fluid)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef FluidSystems::LiquidPhase<Scalar, Dumux::NullComponent<Scalar> > type;
};

/*!
 * \brief The fluid state which is used by the volume variables to
 *        store the thermodynamic state. This should be chosen
 *        appropriately for the model ((non-)isothermal, equilibrium, ...).
 *        This can be done in the problem.
 */
SET_PROP(NavierStokes, FluidState){
    private:
        typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
        typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    public:
        typedef Dumux::ImmiscibleFluidState<Scalar, FluidSystem> type;
};

// set phase index
SET_INT_PROP(NavierStokes, PhaseIdx, 0);

// disable velocity output by default
SET_BOOL_PROP(NavierStokes, VtkAddVelocity, true);

// enable gravity by default
SET_BOOL_PROP(NavierStokes, ProblemEnableGravity, true);

SET_BOOL_PROP(NavierStokes, EnableInertiaTerms, true);

SET_BOOL_PROP(NavierStokes, EnableEnergyTransport, false);

SET_BOOL_PROP(NavierStokes, EnableComponentTransport, false);

//! Normalize the pressure term in the momentum balance or not
SET_BOOL_PROP(NavierStokes, NormalizePressure, true);

SET_PROP(NavierStokes, BoundaryValues)
{
private:
    using CellCenterBoundaryValues = typename GET_PROP_TYPE(TypeTag, CellCenterPrimaryVariables);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FaceBoundaryValues = Dune::FieldVector<typename GET_PROP_TYPE(TypeTag, Scalar),
                                                 GridView::dimension>;
public:
    using type = StaggeredPrimaryVariables<TypeTag, CellCenterBoundaryValues, FaceBoundaryValues>;
};

//! Boundary types at a single degree of freedom
SET_TYPE_PROP(NavierStokes,
              BoundaryTypes,
              StaggeredFreeFlowBoundaryTypes<GET_PROP_VALUE(TypeTag, NumEq)>);

SET_TYPE_PROP(NavierStokes, VtkOutputModule, FreeFlowStaggeredVtkOutputModule<TypeTag>);

SET_TYPE_PROP(NavierStokes, VelocityOutput, StaggeredFreeFlowVelocityOutput<TypeTag>);

SET_TYPE_PROP(NavierStokes, EnergyLocalResidual, FreeFlowEnergyLocalResidual<TypeTag>);

SET_TYPE_PROP(NavierStokes, EnergyFluxVariables, FreeFlowEnergyFluxVariables<TypeTag>);

//! average is used as default model to compute the effective thermal heat conductivity
// SET_PROP(NavierStokesNI, ThermalConductivityModel)
// { private :
//     typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
//   public:
//     typedef ThermalConductivityAverage<Scalar> type;
// };

//////////////////////////////////////////////////////////////////
// Property values for isothermal model required for the general non-isothermal model
//////////////////////////////////////////////////////////////////

// set isothermal Model
SET_TYPE_PROP(NavierStokesNI, IsothermalModel, NavierStokesModel<TypeTag>);

//set isothermal Indices
SET_TYPE_PROP(NavierStokesNI, IsothermalIndices, NavierStokesCommonIndices<TypeTag>);

//set isothermal NumEq
SET_INT_PROP(NavierStokesNI, IsothermalNumEqCellCenter, 1); //!< set the number of equations to 1
SET_INT_PROP(NavierStokesNI, IsothermalNumEqFace, 1); //!< set the number of equations to 1

// \}
} // end namespace Properties

} // end namespace Dumux

#endif
