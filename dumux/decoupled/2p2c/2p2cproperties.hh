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

//Dune-includes
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>

//DUMUX includes
#include <dumux/decoupled/2p/2pproperties.hh>

namespace Dumux
{

////////////////////////////////
// forward declarations
////////////////////////////////


template<class TypeTag>
class VariableClass2P2C;

template<class TypeTag>
class DecTwoPTwoCFluidState;

template <class TypeTag>
struct TwoPTwoCIndices;

////////////////////////////////
// properties
////////////////////////////////
namespace Properties
{

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////

//! The type tag for the two-phase problems
NEW_TYPE_TAG(DecoupledTwoPTwoC, INHERITS_FROM(IMPET));

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////

NEW_PROP_TAG ( TwoPTwoCIndices );
NEW_PROP_TAG( SpatialParameters ); //!< The type of the soil properties object
NEW_PROP_TAG( EnableGravity); //!< Returns whether gravity is considered in the problem
NEW_PROP_TAG( PressureFormulation); //!< The formulation of the model
NEW_PROP_TAG( SaturationFormulation); //!< The formulation of the model
NEW_PROP_TAG( VelocityFormulation); //!< The formulation of the model
NEW_PROP_TAG( EnableCompressibility);// ! Returns whether compressibility is allowed
NEW_PROP_TAG(EnableCapillarity); //!< Returns whether capillarity is regarded
NEW_PROP_TAG( BoundaryMobility );
NEW_PROP_TAG( NumDensityTransport );
NEW_PROP_TAG( FluidSystem );
NEW_PROP_TAG( FluidState );
NEW_PROP_TAG( mpfa ); // Two-point flux approximation (false) or mpfa (true)
//////////////////////////////////////////////////////////////////
// Properties
//////////////////////////////////////////////////////////////////
SET_PROP(DecoupledTwoPTwoC, TwoPTwoCIndices)
{
  typedef TwoPTwoCIndices<TypeTag> type;
};

SET_INT_PROP(DecoupledTwoP, NumEq, 2);

// set fluid/component information
SET_PROP(DecoupledTwoPTwoC, NumPhases) //!< The number of phases in the 2p model is 2
{
    // the property is created in decoupledproperties.hh
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;

public:
    static const int value = FluidSystem::numPhases;
    static_assert(value == 2,
                  "Only fluid systems with 2 phases are supported by the 2p2c model!");
};

SET_PROP(DecoupledTwoPTwoC, NumComponents) //!< The number of components in the 2p2c model is 2
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;

public:
    static const int value = FluidSystem::numComponents;
    static_assert(value == 2,
                  "Only fluid systems with 2 components are supported by the 2p2c model!");
};

//! Set the default formulation
SET_INT_PROP(DecoupledTwoPTwoC,
        PressureFormulation,
        GET_PROP_TYPE(TypeTag, PTAG(TwoPTwoCIndices))::pressureW);

SET_INT_PROP(DecoupledTwoPTwoC,
        SaturationFormulation,
        GET_PROP_TYPE(TypeTag, PTAG(TwoPTwoCIndices))::saturationW);

SET_INT_PROP(DecoupledTwoPTwoC,
        VelocityFormulation,
        GET_PROP_TYPE(TypeTag, PTAG(TwoPTwoCIndices))::velocityW);

SET_PROP(DecoupledTwoPTwoC, TransportSolutionType)
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename Dune::BlockVector<Dune::BlockVector<Dune::FieldVector<Scalar,1> > > type;//!<type for vector of vector (of scalars)

};

SET_BOOL_PROP(DecoupledTwoPTwoC, EnableCompressibility, true);

SET_BOOL_PROP(DecoupledTwoPTwoC, EnableCapillarity, false);


SET_PROP_DEFAULT(BoundaryMobility)
{
    static const int value = TwoPTwoCIndices<TypeTag>::satDependent;
};
SET_PROP_DEFAULT(NumDensityTransport)
{
    static const bool value = false;
};

SET_TYPE_PROP(DecoupledTwoPTwoC, Variables, VariableClass2P2C<TypeTag>);

SET_TYPE_PROP(DecoupledTwoPTwoC, FluidState, DecTwoPTwoCFluidState<TypeTag>);

SET_BOOL_PROP(DecoupledTwoPTwoC, mpfa, false);

}

/*!
 * \brief The common indices for the 2p2c model.
 *
 * The indices are all of the 2p model plus boundary condition flags
 * distinguishing between given composition or saturation on the boundary.
 */
template <class TypeTag>
struct TwoPTwoCIndices : TwoPCommonIndicesDecoupled<TypeTag>
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;

public:
    // BoundaryCondition flags
    static const int satDependent = 0;
    static const int permDependent = 1;
};

// \}

}

#endif
