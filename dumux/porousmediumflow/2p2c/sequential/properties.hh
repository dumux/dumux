// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \ingroup SequentialTwoPTwoCModel
 * \brief Defines the properties required for the sequential 2p2c models.
 */
#ifndef DUMUX_2P2CPROPERTIES_HH
#define DUMUX_2P2CPROPERTIES_HH

#include <dumux/material/fluidstates/compositional.hh>
#include <dumux/porousmediumflow/sequential/pressureproperties.hh>
#include <dumux/porousmediumflow/sequential/transportproperties.hh>
#include <dumux/porousmediumflow/sequential/impetproperties.hh>
#include <dumux/porousmediumflow/2p/sequential/properties.hh>
#include <dumux/material/spatialparams/sequentialfv.hh>

namespace Dumux {

//****** forward declarations  ******//
template<class TypeTag>
class VariableClass;

template<class TypeTag>
class CellData2P2C;

template <class TypeTag>
struct SequentialTwoPTwoCIndices;

////////////////////////////////
// properties
////////////////////////////////
namespace Properties {

//////////////////////////////////////////////////////////////////
// Type tags
//////////////////////////////////////////////////////////////////
//! The type tag for the compositional two-phase problems
// Create new type tags
namespace TTag {
struct SequentialTwoPTwoC { using InheritsFrom = std::tuple<Transport, IMPET, Pressure>; };
} // end namespace TTag

//////////////////////////////////////////////////////////////////
// Property tags
//////////////////////////////////////////////////////////////////
template<class TypeTag, class MyTypeTag>
struct  EnableCapillarity { using type = UndefinedProperty; }; //!< Returns whether capillarity is regarded
template<class TypeTag, class MyTypeTag>
struct  BoundaryMobility  { using type = UndefinedProperty; }; //!< Returns whether mobility or saturation is used for Dirichlet B.C.
//! A minimum permeability can be assigned via the runtime-Parameter SpatialParams.minBoundaryPermeability
template<class TypeTag, class MyTypeTag>
struct  RegulateBoundaryPermeability  { using type = UndefinedProperty; };
}}

//DUMUX includes
#include <dumux/porousmediumflow/2p/sequential/indices.hh>
#include <dumux/porousmediumflow/2p2c/sequential/celldata.hh>

namespace Dumux {
namespace Properties {
//////////////////////////////////////////////////////////////////
// Properties
//////////////////////////////////////////////////////////////////

template<class TypeTag>
struct Indices<TypeTag, TTag::SequentialTwoPTwoC> { using type = SequentialTwoPTwoCIndices<TypeTag>; };

template<class TypeTag>
struct NumEq<TypeTag, TTag::SequentialTwoPTwoC> { static constexpr int value = 3; };

// set fluid/component information
template<class TypeTag>
struct NumPhases<TypeTag, TTag::SequentialTwoPTwoC> //!< The number of phases in the 2p2c model is 2
{
    // the property is declared in dumux/porousmediumflow/sequential/properties.hh
private:
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

public:
    static const int value = FluidSystem::numPhases;
    static_assert(value == 2,
                  "Only fluid systems with 2 phases are supported by the 2p2c model!");
};

template<class TypeTag>
struct NumComponents<TypeTag, TTag::SequentialTwoPTwoC> //!< The number of components in the 2p2c model is 2
{
private:
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

public:
    static const int value = FluidSystem::numComponents;
    static_assert(value == 2,
                  "Only fluid systems with 2 components are supported by the 2p2c model!");
};

//! Set the default formulation
template<class TypeTag>
struct PressureFormulation<TypeTag, TTag::SequentialTwoPTwoC> { static constexpr int value = GetPropType<TypeTag, Properties::Indices>::pressureN; }; //!< Compositional models are very likely compressible

template<class TypeTag>
struct SaturationFormulation<TypeTag, TTag::SequentialTwoPTwoC> { static constexpr int value = GetPropType<TypeTag, Properties::Indices>::saturationW; }; //!< Compositional models are very likely compressible

template<class TypeTag>
struct VelocityFormulation<TypeTag, TTag::SequentialTwoPTwoC> { static constexpr int value = GetPropType<TypeTag, Properties::Indices>::velocityW; }; //!< Compositional models are very likely compressible

template<class TypeTag>
struct TransportSolutionType<TypeTag, TTag::SequentialTwoPTwoC>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    //! type for vector of vector (of scalars)
    using type = Dune::BlockVector<Dune::BlockVector<Dune::FieldVector<Scalar,1> > >;

};

template<class TypeTag>
struct EnableCompressibility<TypeTag, TTag::SequentialTwoPTwoC> { static constexpr bool value = true; }; //!< Compositional models are very likely compressible
template<class TypeTag>
struct VisitFacesOnlyOnce<TypeTag, TTag::SequentialTwoPTwoC> { static constexpr bool value = false; }; //!< Faces are regarded from both sides
template<class TypeTag>
struct EnableCapillarity<TypeTag, TTag::SequentialTwoPTwoC> { static constexpr bool value = false; }; //!< Capillarity is enabled

template<class TypeTag>
struct BoundaryMobility<TypeTag, TTag::SequentialTwoPTwoC> //!< Saturation scales flux on Dirichlet B.C.
{    static const int value = SequentialTwoPTwoCIndices<TypeTag>::satDependent;};

template<class TypeTag>
struct Variables<TypeTag, TTag::SequentialTwoPTwoC> { using type = VariableClass<TypeTag>; };
template<class TypeTag>
struct CellData<TypeTag, TTag::SequentialTwoPTwoC> { using type = CellData2P2C<TypeTag>; };
template<class TypeTag>
struct FluidState<TypeTag, TTag::SequentialTwoPTwoC> { using type = CompositionalFluidState<GetPropType<TypeTag, Properties::Scalar>, GetPropType<TypeTag, Properties::FluidSystem>>; };


//! The spatial parameters to be employed.
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::SequentialTwoPTwoC> { using type = SequentialFVSpatialParams<TypeTag>; };
//! Switch off permeability regularization at Dirichlet boundaries by default.
template<class TypeTag>
struct RegulateBoundaryPermeability<TypeTag, TTag::SequentialTwoPTwoC> { static constexpr bool value = false; };
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
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

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

} // end namespace Dumux

#endif
