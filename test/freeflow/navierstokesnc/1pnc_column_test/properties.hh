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
 * \ingroup NavierStokesNCTests
 * \brief The properties of the channel flow test for the multi-component staggered grid (Navier-)Stokes model.
 */
#ifndef DUMUX_COLUMN_1PNC_TEST_PROPERTIES_HH
#define DUMUX_COLUMN_1PNC_TEST_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/freeflow/navierstokes/problem.hh>
#include <dumux/freeflow/navierstokes/momentum/model.hh>
#include <dumux/freeflow/navierstokes/mass/1pnc/model.hh>

#include <dumux/discretization/staggered/freeflow/properties.hh>
#include <dumux/discretization/fcstaggered.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/multidomain/staggeredfreeflow/couplingmanager.hh>

//Uncomment if maxwell stefan should be used
//#include <dumux/flux/maxwellstefanslaw.hh>

//#include <dumux/freeflow/compositional/navierstokesncmodel.hh>

#include "co2tables.hh"
#include <dumux/material/fluidsystems/brineco2.hh>
//#include <dumux/material/fluidsystems/h2oco2.hh>
#include <dumux/material/fluidsystems/1padapter.hh>

#include "problem.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct ColumnNCTest { };
struct ColumnNCMomentum { using InheritsFrom = std::tuple<ColumnNCTest, NavierStokesMomentum, FaceCenteredStaggeredModel>; };
struct ColumnNCMass { using InheritsFrom = std::tuple<ColumnNCTest, NavierStokesMassOnePNC, CCTpfaModel>; };
//using InheritsFrom = std::tuple<NavierStokesNC, StaggeredFreeFlowModel>; };
} // end namespace TTag

// Set the problem property
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::ColumnNCTest>
{
private:
    using Traits = MultiDomainTraits<TTag::ColumnNCMomentum, TTag::ColumnNCMass>;
public:
    using type = StaggeredFreeFlowCouplingManager<Traits>;
};


template<class TypeTag>
struct ReplaceCompEqIdx<TypeTag, TTag::ColumnNCTest> { static constexpr int value = 0; };

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::ColumnNCTest> { using type = Dune::YaspGrid<2, Dune::TensorProductCoordinates<double, 2> >; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::ColumnNCTest> { using type = Dumux::ColumnNCTestProblem<TypeTag> ; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::ColumnNCTest> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::ColumnNCTest> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::ColumnNCTest> { static constexpr bool value = true; };

//here we can define if mass or molar balance is solved. You need to adapt the boundary conditions if you change that. Now it is in mole fraction
template<class TypeTag>
struct UseMoles<TypeTag, TTag::ColumnNCTest> { static constexpr bool value = true; };

//! Here we set FicksLaw or MaxwellStefansLaw, per default, Ficks law is used. Uncomment to use maxwell stefan. For 2 components they are the same though.
// template<class TypeTag>
// struct MolecularDiffusionType<TypeTag, TTag::ColumnNCTest> { using type = MaxwellStefansLaw<TypeTag>; };

//set the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::ColumnNCTest>
{
    using Scalar = GetPropType<TypeTag, Scalar>;
    using BrineCO2 = FluidSystems::BrineCO2<Scalar,
                                            HeterogeneousCO2Tables::CO2Tables,
                                            Components::TabulatedComponent<Components::H2O<Scalar>>,
                                            FluidSystems::BrineCO2DefaultPolicy</*constantSalinity=*/true, /*simpleButFast=*/true>>;
    static constexpr int phaseIdx = BrineCO2::liquidPhaseIdx;
    using type = FluidSystems::OnePAdapter<BrineCO2, phaseIdx>;
};


} // end namespace Dumux::Properties

#endif
