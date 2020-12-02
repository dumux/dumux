// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
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
 * \brief test problem for diffusion models from the FVCA6 benchmark.
 */
#ifndef DUMUX_TEST_DIFFUSION_3D_PROBLEM_PROPERTIES_HH
#define DUMUX_TEST_DIFFUSION_3D_PROBLEM_PROPERTIES_HH

#include <dune/alugrid/grid.hh>

#include <dumux/material/components/constant.hh>

#include <dumux/common/properties.hh>
#include <dumux/porousmediumflow/2p/sequential/diffusion/cellcentered/pressureproperties.hh>
#include <dumux/porousmediumflow/2p/sequential/diffusion/mpfa/lmethod/3dpressureproperties.hh>
#include <dumux/porousmediumflow/2p/sequential/diffusion/mimetic/pressureproperties.hh>

#include <dumux/porousmediumflow/sequential/cellcentered/velocity.hh>

#include <dumux/linear/seqsolverbackend.hh>

#include "test_diffusionproblem3d.hh"
#include "test_diffusionspatialparams3d.hh"

namespace Dumux::Properties
{
// Create new type tags
namespace TTag {
struct DiffusionTest { using InheritsFrom = std::tuple<TestDiffusionSpatialParams3d, SequentialTwoP>; };

// set the types for the 2PFA FV method
struct FVTest { using InheritsFrom = std::tuple<DiffusionTest, FVPressureTwoP>; };

// set the types for the MPFA-L FV method
struct FVMPFAL3DTestTypeTag { using InheritsFrom = std::tuple<DiffusionTest, FvMpfaL3dPressureTwoP>; };

// set the types for the mimetic FD method
struct MimeticTest { using InheritsFrom = std::tuple<DiffusionTest, MimeticPressureTwoP>; };

} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::DiffusionTest> { using type = Dune::ALUGrid<3, 3, Dune::cube, Dune::nonconforming>; };

template<class TypeTag>
struct Problem<TypeTag, TTag::DiffusionTest> { using type = TestDiffusion3DProblem<TypeTag>; };

// Set the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::DiffusionTest>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using WettingPhase = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
    using NonwettingPhase = FluidSystems::OnePLiquid<Scalar, Components::Constant<1, Scalar> >;
    using type = FluidSystems::TwoPImmiscible<Scalar, WettingPhase, NonwettingPhase>;
};

template<class TypeTag>
struct LinearSolver<TypeTag, TTag::DiffusionTest> { using type = UMFPackBackend; };

} // end namespace Properties

#endif
