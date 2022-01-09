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
 * \brief The properties of the test for the compositional staggered grid (Navier-)Stokes model.
 */
#ifndef DUMUX_DENSITY_FLOW_NC_TEST_PROPERTIES_HH
#define DUMUX_DENSITY_FLOW_NC_TEST_PROPERTIES_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/freeflow/navierstokes/problem.hh>
#include <dumux/freeflow/navierstokes/momentum/model.hh>
#include <dumux/freeflow/navierstokes/mass/1pnc/model.hh>

#include <dumux/discretization/staggered/freeflow/properties.hh>
#include <dumux/discretization/fcstaggered.hh>
#include <dumux/discretization/cctpfa.hh>
#include <dumux/multidomain/staggeredfreeflow/couplingmanager.hh>

#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/1padapter.hh>
#include <dumux/material/fluidsystems/h2oair.hh>

#include "problem.hh"

namespace Dumux::Properties {

// Create new type tags
namespace TTag {
struct DensityDrivenFlowTest {};
struct DensityDrivenFlowMomentum { using InheritsFrom = std::tuple<DensityDrivenFlowTest, NavierStokesMomentum, FaceCenteredStaggeredModel>; };
#if !NONISOTHERMAL
struct DensityDrivenFlowMass { using InheritsFrom = std::tuple<DensityDrivenFlowTest, NavierStokesMassOnePNC, CCTpfaModel>; };
#else
struct DensityDrivenFlowMass { using InheritsFrom = std::tuple<DensityDrivenFlowTest, NavierStokesMassOnePNCNI, CCTpfaModel>; };
#endif
} // end namespace TTag

// Set the problem property
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::DensityDrivenFlowTest>
{
private:
    using Traits = MultiDomainTraits<TTag::DensityDrivenFlowMomentum, TTag::DensityDrivenFlowMass>;
public:
    using type = StaggeredFreeFlowCouplingManager<Traits>;
};

// Select the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::DensityDrivenFlowTest>
{
    using H2OAir = FluidSystems::H2OAir<GetPropType<TypeTag, Properties::Scalar>>;
    static constexpr int phaseIdx = H2OAir::liquidPhaseIdx;
    using type = FluidSystems::OnePAdapter<H2OAir, phaseIdx>;
};

template<class TypeTag>
struct ReplaceCompEqIdx<TypeTag, TTag::DensityDrivenFlowTest> { static constexpr int value = 0; };

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::DensityDrivenFlowTest> { using type = Dune::YaspGrid<2>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::DensityDrivenFlowTest> { using type = Dumux::DensityDrivenFlowProblem<TypeTag> ; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::DensityDrivenFlowTest> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::DensityDrivenFlowTest> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::DensityDrivenFlowTest> { static constexpr bool value = true; };

template<class TypeTag>
struct UseMoles<TypeTag, TTag::DensityDrivenFlowTest> { static constexpr bool value = true; };

} // end namespace Dumux::Properties

#endif
