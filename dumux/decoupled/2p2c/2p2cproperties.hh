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
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
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
class TwoPTwoCFluidState;

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
NEW_PROP_TAG( SpatialParams ); //!< The type of the soil properties object
NEW_PROP_TAG( SpatialParameters ); //!< DEPRECATED The old type of the soil properties object
NEW_PROP_TAG( ProblemEnableGravity); //!< Returns whether gravity is considered in the problem
NEW_PROP_TAG( EnableGravity); //!< DEPRECATED Returns whether gravity is considered in the problem
NEW_PROP_TAG( PressureFormulation); //!< The formulation of the model
NEW_PROP_TAG( SaturationFormulation); //!< The formulation of the model
NEW_PROP_TAG( VelocityFormulation); //!< The formulation of the model
NEW_PROP_TAG( EnableCompressibility); //!< Returns whether compressibility is allowed
NEW_PROP_TAG( EnableCapillarity); //!< Returns whether capillarity is regarded
NEW_PROP_TAG( BoundaryMobility ); //!< Returns whether mobility or saturation is used for Dirichlet B.C.
NEW_PROP_TAG( ImpetRestrictFluxInTransport ); //!< Restrict flux if direction reverses after pressure equation
NEW_PROP_TAG( RestrictFluxInTransport ); //!< DEPRECATED Restrict flux if direction reverses after pressure equation
NEW_PROP_TAG( ImpetErrorTermFactor ); //!< Damping factor \f$ \alpha \f$ in pressure equation
NEW_PROP_TAG( ErrorTermFactor ); //!< DEPRECATED Damping factor \f$ \alpha \f$ in pressure equation
NEW_PROP_TAG( ImpetErrorTermLowerBound ); //!< Upper bound for regularized error damping
NEW_PROP_TAG( ErrorTermLowerBound ); //!< DEPRECATED Upper bound for regularized error damping
NEW_PROP_TAG( ImpetErrorTermUpperBound ); //!< Lower bound where error is not corrected
NEW_PROP_TAG( ErrorTermUpperBound ); //!< DEPRECATED Lower bound where error is not corrected
NEW_PROP_TAG( FluidSystem ); //!< The fluid system
NEW_PROP_TAG( FluidState ); //!< The fluid state
NEW_PROP_TAG( ImpetEnableVolumeIntegral ); //!< Enables volume integral in the pressure equation (volume balance formulation)
NEW_PROP_TAG( EnableVolumeIntegral ); //!< DEPRECATED Enables volume integral in the pressure equation (volume balance formulation)
NEW_PROP_TAG( GridAdaptEnableMultiPointFluxApproximation); //!< HangingNode: Two-point flux approximation (false) or mpfa (true)
NEW_PROP_TAG( EnableMultiPointFluxApproximationOnAdaptiveGrids ); //Deprecated
NEW_PROP_TAG( EnableSecondHalfEdge ); //!< Uses second interaction volume for second half-edge in 2D
NEW_PROP_TAG( GridAdaptEnableSecondHalfEdge ); //!< Uses second interaction volume for second half-edge in 2D
}}

//DUMUX includes
#include <dumux/decoupled/2p/2pindices.hh>
#include <dumux/decoupled/2p2c/cellData2p2c.hh>
#include <dumux/decoupled/2p2c/2p2cfluidstate.hh>

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
        GET_PROP_TYPE(TypeTag, Indices)::pressureNW);

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

SET_BOOL_PROP(DecoupledTwoPTwoC, EnableCompressibility, true); //!< Compositional models are very likely compressible
SET_BOOL_PROP(DecoupledTwoPTwoC, VisitFacesOnlyOnce, false); //!< Faces are regarded from both sides
SET_BOOL_PROP(DecoupledTwoPTwoC, EnableCapillarity, false); //!< Capillarity is enabled
SET_INT_PROP(DecoupledTwoPTwoC, ImpetRestrictFluxInTransport, GET_PROP_VALUE(TypeTag, RestrictFluxInTransport)); //!< Restrict (no upwind) flux in transport step if direction reverses after pressure equation
SET_INT_PROP(DecoupledTwoPTwoC, RestrictFluxInTransport, 0); //!< DEPRECATED Restrict (no upwind) flux in transport step if direction reverses after pressure equation

SET_PROP(DecoupledTwoPTwoC, BoundaryMobility) //!< Saturation scales flux on Dirichlet B.C.
{    static const int value = DecoupledTwoPTwoCIndices<TypeTag>::satDependent;};

SET_TYPE_PROP(DecoupledTwoPTwoC, Variables, VariableClass<TypeTag>);
SET_TYPE_PROP(DecoupledTwoPTwoC, CellData, CellData2P2C<TypeTag>);
SET_TYPE_PROP(DecoupledTwoPTwoC, FluidState, TwoPTwoCFluidState<TypeTag>);


SET_TYPE_PROP(DecoupledTwoPTwoC, SpatialParameters, typename GET_PROP_TYPE(TypeTag, SpatialParams));//!< DEPRECATED SpatialParameters property

SET_BOOL_PROP(DecoupledTwoPTwoC, GridAdaptEnableMultiPointFluxApproximation,
        GET_PROP_VALUE(TypeTag,EnableMultiPointFluxApproximationOnAdaptiveGrids)); //!< MPFA disabled on adaptive grids
SET_BOOL_PROP(DecoupledTwoPTwoC, EnableMultiPointFluxApproximationOnAdaptiveGrids, false); //!<  DEPRECATED MPFA disabled on adaptive grids
SET_BOOL_PROP(DecoupledTwoPTwoC, ImpetEnableVolumeIntegral, GET_PROP_VALUE(TypeTag,EnableVolumeIntegral)); //!< Regard volume integral in pressure equation
SET_BOOL_PROP(DecoupledTwoPTwoC, EnableVolumeIntegral, true); //!< DEPRECATED Regard volume integral in pressure equation

SET_SCALAR_PROP(DecoupledTwoPTwoC, ImpetErrorTermFactor, GET_PROP_VALUE(TypeTag, ErrorTermFactor)); //!< Damping factor \f$ \alpha \f$ in pressure equation
SET_SCALAR_PROP(DecoupledTwoPTwoC, ErrorTermFactor, 0.5); //!< Damping factor \f$ \alpha \f$ in pressure equation
SET_SCALAR_PROP(DecoupledTwoPTwoC, ImpetErrorTermLowerBound, GET_PROP_VALUE(TypeTag, ErrorTermLowerBound)); //!< Lower bound where error is not corrected
SET_SCALAR_PROP(DecoupledTwoPTwoC, ErrorTermLowerBound, 0.2); //!< Lower bound where error is not corrected
SET_SCALAR_PROP(DecoupledTwoPTwoC, ImpetErrorTermUpperBound, GET_PROP_VALUE(TypeTag, ErrorTermUpperBound)); //!< Upper bound for regularized error damping
SET_SCALAR_PROP(DecoupledTwoPTwoC, ErrorTermUpperBound, 0.9); //!< Upper bound for regularized error damping

//Has to be removed if DEPRECATED EnableGravity is removed!
SET_BOOL_PROP(DecoupledTwoPTwoC, ProblemEnableGravity, GET_PROP_VALUE(TypeTag, EnableGravity));
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

    // Ensure pressure fomrulation index coincides with FluidSystem
    static const int pressureW = wPhaseIdx;
    static const int pressureNW = nPhaseIdx;

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
