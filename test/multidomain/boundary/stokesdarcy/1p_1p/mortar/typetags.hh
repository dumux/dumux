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
 * \ingroup BoundaryTests
 * \brief TODO doc me
 */
#include <dune/foamgrid/foamgrid.hh>
#include <dune/functions/functionspacebases/lagrangebasis.hh>

#include <dumux/common/properties.hh>

#include "problem_darcy.hh"
#include "problem_stokes.hh"

#include "fluxprojector.hh"
#include "pressureprojector.hh"

#include "fluxpreconditioner.hh"
#include "pressurepreconditioner.hh"

namespace Dumux {

// Traits class for a sub-domain
template<class SubDomainTypeTag>
struct SubDomainTraits
{
    using SubDomainSolutionVector = Dumux::GetPropType<SubDomainTypeTag, Dumux::Properties::SolutionVector>;
    using SubDomainGridGeometry = Dumux::GetPropType<SubDomainTypeTag, Dumux::Properties::FVGridGeometry>;
    using SubDomainGridVariables = Dumux::GetPropType<SubDomainTypeTag, Dumux::Properties::GridVariables>;
    using SubDomainFluxVariables = Dumux::GetPropType<SubDomainTypeTag, Dumux::Properties::FluxVariables>;
};

// Define grid and basis for mortar domain
struct MortarSpaceTraits
{
    using Scalar = double;
    using Grid = Dune::FoamGrid<1, 2>;
    using GridView = typename Grid::LeafGridView;
    using SolutionVector = Dune::BlockVector<Dune::FieldVector<Scalar, 1>>;
    using FEBasis = Dune::Functions::LagrangeBasis<GridView, 0>;
};

// Projector traits
template<class SubDomainTypeTag>
struct MortarProjectorTraits : public SubDomainTraits<SubDomainTypeTag>
{
    using MortarFEBasis = typename MortarSpaceTraits::FEBasis;
    using MortarSolutionVector = typename MortarSpaceTraits::SolutionVector;
};

// Projection operators
template<class SubDomainTypeTag>
struct FluxProjectorTraits { using type = Dumux::MortarFluxProjector< MortarProjectorTraits<SubDomainTypeTag> >; };
template<class SubDomainTypeTag>
struct PressureProjectorTraits { using type = Dumux::MortarPressureProjector< MortarProjectorTraits<SubDomainTypeTag> >; };

namespace Properties {

// create new type tag nodes for pressure/flux coupling
namespace TTag {
struct DarcyOnePTpfaPressure { using InheritsFrom = std::tuple<DarcyOnePTpfa>; };
struct DarcyOnePMpfaPressure { using InheritsFrom = std::tuple<DarcyOnePMpfa>; };
struct DarcyOnePBoxPressure { using InheritsFrom = std::tuple<DarcyOnePBox>; };

struct DarcyOnePTpfaFlux { using InheritsFrom = std::tuple<DarcyOnePTpfa>; };
struct DarcyOnePMpfaFlux { using InheritsFrom = std::tuple<DarcyOnePMpfa>; };
struct DarcyOnePBoxFlux { using InheritsFrom = std::tuple<DarcyOnePBox>; };

struct StokesOnePPressure { using InheritsFrom = std::tuple<StokesOneP>; };
struct StokesOnePFlux { using InheritsFrom = std::tuple<StokesOneP>; };
} // end namespace TTag

// set projector properties
template<class TypeTag>
struct MortarProjector<TypeTag, TTag::DarcyOnePTpfaPressure> { using type = typename PressureProjectorTraits<TypeTag>::type; };
template<class TypeTag>
struct MortarProjector<TypeTag, TTag::DarcyOnePMpfaPressure> { using type = typename PressureProjectorTraits<TypeTag>::type; };
template<class TypeTag>
struct MortarProjector<TypeTag, TTag::DarcyOnePBoxPressure> { using type = typename PressureProjectorTraits<TypeTag>::type; };

template<class TypeTag>
struct MortarProjector<TypeTag, TTag::DarcyOnePTpfaFlux> { using type = typename FluxProjectorTraits<TypeTag>::type; };
template<class TypeTag>
struct MortarProjector<TypeTag, TTag::DarcyOnePMpfaFlux> { using type = typename FluxProjectorTraits<TypeTag>::type; };
template<class TypeTag>
struct MortarProjector<TypeTag, TTag::DarcyOnePBoxFlux> { using type = typename FluxProjectorTraits<TypeTag>::type; };

template<class TypeTag>
struct MortarProjector<TypeTag, TTag::StokesOnePPressure> { using type = typename PressureProjectorTraits<TypeTag>::type; };
template<class TypeTag>
struct MortarProjector<TypeTag, TTag::StokesOnePFlux> { using type = typename FluxProjectorTraits<TypeTag>::type; };

} // end namespace Properties
} // end namespace Dumux
