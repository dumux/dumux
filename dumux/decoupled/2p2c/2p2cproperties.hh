// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 20010 by Benjamin Faigle                                  *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/

/*!
 * \ingroup IMPEC Properties
 * \file
 *
 * \brief Defines the properties required for the decoupled 2p2c models.
 */
#ifndef DUMUX_2P2CPROPERTIES_HH
#define DUMUX_2P2CPROPERTIES_HH

#include <dumux/decoupled/common/pressureproperties.hh>
#include <dumux/decoupled/common/transportproperties.hh>
#include <dumux/decoupled/common/impetproperties.hh>

namespace Dumux
{

//****** forward declarations  ******//
template<class TypeTag>
class VariableClass;

template<class TypeTag>
class CellData2P2C;

template<class TypeTag>
class DecoupledTwoPTwoCFluidState;

template <class TypeTag>
struct DecoupledTwoPTwoCIndices;

////////////////////////////////
// properties
////////////////////////////////
namespace Properties
{

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////
//! The type tag for the compositional two-phase problems
NEW_TYPE_TAG(DecoupledTwoPTwoC, INHERITS_FROM(Pressure, Transport, IMPET));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////
NEW_PROP_TAG( Indices );
NEW_PROP_TAG( SpatialParameters ); //!< The type of the soil properties object
NEW_PROP_TAG( EnableGravity); //!< Returns whether gravity is considered in the problem
NEW_PROP_TAG( PressureFormulation); //!< The formulation of the model
NEW_PROP_TAG( SaturationFormulation); //!< The formulation of the model
NEW_PROP_TAG( VelocityFormulation); //!< The formulation of the model
NEW_PROP_TAG( EnableCompressibility);// ! Returns whether compressibility is allowed
NEW_PROP_TAG( EnableCapillarity); //!< Returns whether capillarity is regarded
NEW_PROP_TAG( BoundaryMobility );
NEW_PROP_TAG( NumDensityTransport );
NEW_PROP_TAG( ErrorTermFactor );
NEW_PROP_TAG( ErrorTermLowerBound );
NEW_PROP_TAG( ErrorTermUpperBound );
NEW_PROP_TAG( FluidSystem );
NEW_PROP_TAG( FluidState );
NEW_PROP_TAG( EnableMultiPointFluxApproximationOnAdaptiveGrids ); // Two-point flux approximation (false) or mpfa (true)
NEW_PROP_TAG( EnableVolumeIntegral );
NEW_PROP_TAG( EnableSecondHalfEdge );

}}

//DUMUX includes
#include <dumux/decoupled/2p/2pindices.hh>
#include <dumux/decoupled/2p2c/cellData2p2c.hh>
#include <dumux/decoupled/2p2c/dec2p2cfluidstate.hh>

namespace Dumux {
namespace Properties {
//////////////////////////////////////////////////////////////////
// Properties
//////////////////////////////////////////////////////////////////
SET_TYPE_PROP(DecoupledTwoPTwoC, Indices,DecoupledTwoPTwoCIndices<TypeTag>);

SET_INT_PROP(DecoupledTwoPTwoC, NumEq, 3);

// set fluid/component information
SET_PROP(DecoupledTwoPTwoC, NumPhases) //!< The number of phases in the 2p model is 2
{
    // the property is created in decoupledproperties.hh
private:
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

public:
    static const int value = FluidSystem::numPhases;
    static_assert(value == 2,
                  "Only fluid systems with 2 phases are supported by the 2p2c model!");
};

SET_PROP(DecoupledTwoPTwoC, NumComponents) //!< The number of components in the 2p2c model is 2
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

public:
    static const int value = FluidSystem::numComponents;
    static_assert(value == 2,
                  "Only fluid systems with 2 components are supported by the 2p2c model!");
};

//! Set the default formulation
SET_INT_PROP(DecoupledTwoPTwoC,
        PressureFormulation,
        GET_PROP_TYPE(TypeTag, Indices)::pressureW);

SET_INT_PROP(DecoupledTwoPTwoC,
        SaturationFormulation,
        GET_PROP_TYPE(TypeTag, Indices)::saturationW);

SET_INT_PROP(DecoupledTwoPTwoC,
        VelocityFormulation,
        GET_PROP_TYPE(TypeTag, Indices)::velocityW);

SET_PROP(DecoupledTwoPTwoC, TransportSolutionType)
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename Dune::BlockVector<Dune::BlockVector<Dune::FieldVector<Scalar,1> > > type;//!<type for vector of vector (of scalars)

};

SET_BOOL_PROP(DecoupledTwoPTwoC, EnableCompressibility, true);
//! Faces are regarded from both sides
SET_BOOL_PROP(DecoupledTwoPTwoC, VisitFacesOnlyOnce, false);

SET_BOOL_PROP(DecoupledTwoPTwoC, EnableCapillarity, false);


SET_PROP(DecoupledTwoPTwoC, BoundaryMobility)
{
    static const int value = DecoupledTwoPTwoCIndices<TypeTag>::satDependent;
};
SET_PROP(DecoupledTwoPTwoC, NumDensityTransport)
{
    static const bool value = false;
};

SET_TYPE_PROP(DecoupledTwoPTwoC, Variables, VariableClass<TypeTag>);
SET_TYPE_PROP(DecoupledTwoPTwoC, CellData, CellData2P2C<TypeTag>);
SET_TYPE_PROP(DecoupledTwoPTwoC, FluidState, DecoupledTwoPTwoCFluidState<TypeTag>);

SET_BOOL_PROP(DecoupledTwoPTwoC, EnableMultiPointFluxApproximationOnAdaptiveGrids, false);
SET_BOOL_PROP(DecoupledTwoPTwoC, EnableVolumeIntegral, true);

SET_SCALAR_PROP(DecoupledTwoPTwoC, ErrorTermFactor, 0.5);
SET_SCALAR_PROP(DecoupledTwoPTwoC, ErrorTermLowerBound, 0.2);
SET_SCALAR_PROP(DecoupledTwoPTwoC, ErrorTermUpperBound, 0.9);
}

/*!
 * \brief The common indices for the 2p2c model.
 *
 * The indices are all of the 2p model plus boundary condition flags
 * distinguishing between given composition or saturation on the boundary.
 */
template <class TypeTag>
struct DecoupledTwoPTwoCIndices : DecoupledTwoPCommonIndices
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

public:
    // Component indices
    static const int wPhaseIdx = FluidSystem::wPhaseIdx;
    static const int nPhaseIdx = FluidSystem::nPhaseIdx;

    // Component indices
    static const int wCompIdx = wPhaseIdx; //!< Component index equals phase index
    static const int nCompIdx = nPhaseIdx; //!< Component index equals phase index

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
