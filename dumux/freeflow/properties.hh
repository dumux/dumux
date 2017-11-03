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
 * \file
 *
 * \brief Defines a type tag and some properties for free flow models.
 */

#ifndef DUMUX_FREE_FLOW_PROPERTIES_HH
#define DUMUX_FREE_FLOW_PROPERTIES_HH

#include <dumux/io/vtkoutputmodule.hh>

// #include <dumux/FreeFlow/implicit/fluxvariables.hh>
// #include <dumux/FreeFlow/implicit/fluxvariablescache.hh>
// #include <dumux/FreeFlow/nonisothermal/implicit/localresidual.hh>
// #include <dumux/FreeFlow/compositional/primaryvariableswitch.hh>
// #include <dumux/FreeFlow/implicit/velocityoutput.hh>

// #include <dumux/discretization/darcyslaw.hh>
// #include <dumux/discretization/fickslaw.hh>
// #include <dumux/discretization/fourierslaw.hh>
#include <dumux/discretization/staggered/freeflow/staggeredgeometryhelper.hh>
#include <dumux/discretization/staggered/freeflow/subcontrolvolumeface.hh>
#include <dumux/discretization/staggered/freeflow/facevariables.hh>
#include <dumux/implicit/staggered/primaryvariables.hh>

#include "./staggered/boundarytypes.hh"

namespace Dumux
{
namespace Properties
{
//! Type tag for models involving flow in porous media
NEW_TYPE_TAG(FreeFlow);

SET_INT_PROP(FreeFlow, NumEqFace, 1); //!< set the number of equations to 1



//! The sub-controlvolume face
SET_PROP(FreeFlow, SubControlVolumeFace)
{
private:
    using Grid = typename GET_PROP_TYPE(TypeTag, Grid);
    using ScvfGeometry = typename Grid::template Codim<1>::Geometry;
    using IndexType = typename Grid::LeafGridView::IndexSet::IndexType;
public:
    typedef Dumux::StaggeredSubControlVolumeFace<ScvfGeometry, IndexType> type;
};

//! The geometry helper required for the stencils, etc.
SET_PROP(FreeFlow, StaggeredGeometryHelper)
{
private:
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
public:
    using type = StaggeredGeometryHelper<GridView>;
};

//! The variables living on the faces
SET_TYPE_PROP(FreeFlow, FaceVariables, StaggeredFaceVariables<TypeTag>);


SET_PROP(FreeFlow, BoundaryValues)
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
SET_TYPE_PROP(FreeFlow,
              BoundaryTypes,
              StaggeredFreeFlowBoundaryTypes<GET_PROP_VALUE(TypeTag, NumEq)>);






// //! The flux variables for models involving flow in porous media
// SET_TYPE_PROP(FreeFlow, FluxVariables, PorousMediumFluxVariables<TypeTag>);
//
// //! The flux variables cache class for models involving flow in porous media
// SET_TYPE_PROP(FreeFlow, FluxVariablesCache, PorousMediumFluxVariablesCache<TypeTag>);
//
// //! By default, we use darcy's law for the advective fluxes
// SET_TYPE_PROP(FreeFlow, AdvectionType, DarcysLaw<TypeTag>);
//
// //! By default, we use fick's law for the diffusive fluxes
// SET_TYPE_PROP(FreeFlow, MolecularDiffusionType, FicksLaw<TypeTag>);
//
// //! By default, we use fourier's law as the default for heat conduction fluxes
// SET_TYPE_PROP(FreeFlow, HeatConductionType, FouriersLaw<TypeTag>);
//
// //! By default, parameters are solution-dependent
// SET_BOOL_PROP(FreeFlow, SolutionDependentAdvection, true);
// SET_BOOL_PROP(FreeFlow, SolutionDependentMolecularDiffusion, true);
// SET_BOOL_PROP(FreeFlow, SolutionDependentHeatConduction, true);
//
// //! By default, we evaluate the permeability in the volume
// SET_BOOL_PROP(FreeFlow, EvaluatePermeabilityAtScvfIP, false);
//
// //! The default implementation of the energy balance equation for flow problems in porous media.
// SET_TYPE_PROP(FreeFlow, EnergyLocalResidual, EnergyLocalResidual<TypeTag> );

//! Velocity output
// SET_TYPE_PROP(FreeFlow, VelocityOutput, ImplicitVelocityOutput<TypeTag>);

//! By default, we set an empty primary variables switch
// SET_TYPE_PROP(FreeFlow, PrimaryVariableSwitch, NoPrimaryVariableSwitch<TypeTag>);

} // namespace Properties
} // namespace Dumux

 #endif
