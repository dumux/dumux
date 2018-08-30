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
 * \ingroup IMPEC
 * \ingroup IMPETProperties
 *
 * \file
 *
 * \brief Defines the properties required for the sequential 2p2c models.
 */
#ifndef DUMUX_2P2CPROPERTIES_HH
#define DUMUX_2P2CPROPERTIES_HH

#include <dumux/porousmediumflow/sequential/pressureproperties.hh>
#include <dumux/porousmediumflow/sequential/transportproperties.hh>
#include <dumux/porousmediumflow/sequential/impetproperties.hh>
#include <dumux/porousmediumflow/2p/sequential/properties.hh>
#include <dumux/material/spatialparams/sequentialfv.hh>

namespace Dumux
{

//****** forward declarations  ******//
template<class TypeTag>
class VariableClass;

template<class TypeTag>
class CellData2P2C;

template<class Scalar, class FluidSystem>
class TwoPTwoCFluidState;

template <class TypeTag>
struct SequentialTwoPTwoCIndices;

////////////////////////////////
// properties
////////////////////////////////
namespace Properties
{

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////
//! The type tag for the compositional two-phase problems
NEW_TYPE_TAG(SequentialTwoPTwoC, INHERITS_FROM(Pressure, Transport, IMPET));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////
NEW_PROP_TAG( EnableCapillarity); //!< Returns whether capillarity is regarded
NEW_PROP_TAG( BoundaryMobility ); //!< Returns whether mobility or saturation is used for Dirichlet B.C.
//! A minimum permeability can be assigned via the runtime-Parameter SpatialParams.minBoundaryPermeability
NEW_PROP_TAG( RegulateBoundaryPermeability );
}}

//DUMUX includes
#include <dumux/porousmediumflow/2p/sequential/indices.hh>
#include <dumux/porousmediumflow/2p2c/sequential/celldata.hh>

namespace Dumux {
namespace Properties {
//////////////////////////////////////////////////////////////////
// Properties
//////////////////////////////////////////////////////////////////

SET_TYPE_PROP(SequentialTwoPTwoC, Indices,SequentialTwoPTwoCIndices<TypeTag>);

SET_INT_PROP(SequentialTwoPTwoC, NumEq, 3);

// set fluid/component information
SET_PROP(SequentialTwoPTwoC, NumPhases) //!< The number of phases in the 2p2c model is 2
{
    // the property is declared in dumux/porousmediumflow/sequential/properties.hh
private:
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);

public:
    static const int value = FluidSystem::numPhases;
    static_assert(value == 2,
                  "Only fluid systems with 2 phases are supported by the 2p2c model!");
};

SET_PROP(SequentialTwoPTwoC, NumComponents) //!< The number of components in the 2p2c model is 2
{
private:
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);

public:
    static const int value = FluidSystem::numComponents;
    static_assert(value == 2,
                  "Only fluid systems with 2 components are supported by the 2p2c model!");
};

//! Set the default formulation
SET_INT_PROP(SequentialTwoPTwoC,
        PressureFormulation,
        GET_PROP_TYPE(TypeTag, Indices)::pressureN);

SET_INT_PROP(SequentialTwoPTwoC,
        SaturationFormulation,
        GET_PROP_TYPE(TypeTag, Indices)::saturationW);

SET_INT_PROP(SequentialTwoPTwoC,
        VelocityFormulation,
        GET_PROP_TYPE(TypeTag, Indices)::velocityW);

SET_PROP(SequentialTwoPTwoC, TransportSolutionType)
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    //! type for vector of vector (of scalars)
    using type = Dune::BlockVector<Dune::BlockVector<Dune::FieldVector<Scalar,1> > >;

};

SET_BOOL_PROP(SequentialTwoPTwoC, EnableCompressibility, true); //!< Compositional models are very likely compressible
SET_BOOL_PROP(SequentialTwoPTwoC, VisitFacesOnlyOnce, false); //!< Faces are regarded from both sides
SET_BOOL_PROP(SequentialTwoPTwoC, EnableCapillarity, false); //!< Capillarity is enabled

SET_PROP(SequentialTwoPTwoC, BoundaryMobility) //!< Saturation scales flux on Dirichlet B.C.
{    static const int value = SequentialTwoPTwoCIndices<TypeTag>::satDependent;};

SET_TYPE_PROP(SequentialTwoPTwoC, Variables, VariableClass<TypeTag>);
SET_TYPE_PROP(SequentialTwoPTwoC, CellData, CellData2P2C<TypeTag>);
SET_TYPE_PROP(SequentialTwoPTwoC, FluidState, TwoPTwoCFluidState<typename GET_PROP_TYPE(TypeTag, Scalar), typename GET_PROP_TYPE(TypeTag, FluidSystem)>);


//! The spatial parameters to be employed.
SET_TYPE_PROP(SequentialTwoPTwoC, SpatialParams, SequentialFVSpatialParams<TypeTag>);
//! Switch off permeability regularization at Dirichlet boundaries by default.
SET_BOOL_PROP(SequentialTwoPTwoC, RegulateBoundaryPermeability, false);
}

/*!
 * \brief The common indices for the 2p2c model.
 *
 * The indices are all of the 2p model plus boundary condition flags
 * distinguishing between given composition or saturation on the boundary.
 * As we have 3 equations, one pressure and two component transport equations,
 * special equation indices have to be provided for boundary conditions.
 */
template <class TypeTag>
struct SequentialTwoPTwoCIndices : public SequentialTwoPCommonIndices
{
private:
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);

public:
    // Component indices
    static const int wPhaseIdx = FluidSystem::phase0Idx;
    static const int nPhaseIdx = FluidSystem::phase1Idx;

    // Component indices
    static const int wCompIdx = wPhaseIdx; //!< Component index equals phase index
    static const int nCompIdx = nPhaseIdx; //!< Component index equals phase index

    // Ensure pressure formulation index coincides with FluidSystem
    static const int pressureW = wPhaseIdx;
    static const int pressureN = nPhaseIdx;

    // Equation indices
    static const int pressureEqIdx = 0;
    static const int transportEqOffset = pressureEqIdx + 1; //!< Offset to access transport (mass conservation -) equations
    static const int contiWEqIdx = transportEqOffset + wCompIdx; //!< Index of the wetting component transport equation
    static const int contiNEqIdx = transportEqOffset + nCompIdx; //!< Index of the nonwetting component transport equation

    //! Type of value on the Boundary
    enum BoundaryFormulation
        {
            saturation=-1,
            concentration=-2
        };


    // BoundaryCondition flags
    static const int satDependent = 0;
    static const int permDependent = 1;
};

// \}

}

#endif
