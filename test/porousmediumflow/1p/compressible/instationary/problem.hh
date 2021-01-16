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
 * \ingroup OnePTests
 * \brief The properties for the incompressible test
 */

#ifndef DUMUX_COMPRESSIBLE_ONEP_TEST_PROBLEM_HH
#define DUMUX_COMPRESSIBLE_ONEP_TEST_PROBLEM_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/material/components/h2o.hh>
#include <dumux/material/components/tabulatedcomponent.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include <dumux/common/boundarytypes.hh>

#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/ccmpfa.hh>
#include <dumux/discretization/box.hh>

#include <dumux/porousmediumflow/problem.hh>
#include <dumux/porousmediumflow/1p/model.hh>

#include "spatialparams.hh"

namespace Dumux {
// forward declarations
template<class TypeTag> class OnePTestProblem;

namespace Properties {
// create the type tag nodes
// Create new type tags
namespace TTag {
struct OnePCompressible { using InheritsFrom = std::tuple<OneP>; };
struct OnePCompressibleTpfa { using InheritsFrom = std::tuple<OnePCompressible, CCTpfaModel>; };
struct OnePCompressibleMpfa { using InheritsFrom = std::tuple<OnePCompressible, CCMpfaModel>; };
struct OnePCompressibleBox { using InheritsFrom = std::tuple<OnePCompressible, BoxModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::OnePCompressible> { using type = Dune::YaspGrid<2>; };

// Set the problem type
template<class TypeTag>
struct Problem<TypeTag, TTag::OnePCompressible> { using type = OnePTestProblem<TypeTag>; };

// set the spatial params
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::OnePCompressible>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = OnePTestSpatialParams<GridGeometry, Scalar>;
};

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::OnePCompressible>
{
private:
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
public:
    using type = FluidSystems::OnePLiquid<Scalar, Components::TabulatedComponent<Components::H2O<Scalar>>>;
};

// Disable caching (for testing purposes)
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::OnePCompressible> { static constexpr bool value = false; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::OnePCompressible> { static constexpr bool value = false; };
template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::OnePCompressible> { static constexpr bool value = false; };
} // end namespace Properties

/*!
 * \ingroup OnePTests
 * \brief  Test problem for the compressible one-phase model.
 *
 * Can be run as <tt>./test_box1pfv</tt> or
 * <tt>./test_cc1pfv</tt>
 */
template<class TypeTag>
class OnePTestProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    static constexpr int dimWorld = GridView::dimensionworld;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    OnePTestProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        Components::TabulatedComponent<Components::H2O<Scalar>>::init(272.15, 294.15, 10,
                                                      1.0e4, 1.0e6, 200);
    }

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary control volume.
     *
     * \param globalPos The position of the center of the finite volume
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;

        Scalar eps = 1.0e-6;
        if (globalPos[dimWorld-1] < eps || globalPos[dimWorld-1] > this->gridGeometry().bBoxMax()[dimWorld-1] - eps)
            values.setAllDirichlet();
        else
            values.setAllNeumann();

        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet control volume.
     *
     * \param globalPos The center of the finite volume which ought to be set.
     *
     * For this method, the \a values parameter stores primary variables.
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values(0);
        values[0] = 1.0e+5*(2.0 - globalPos[dimWorld-1]);
        return values;
    }

    /*!
     * \brief Evaluates the initial conditions
     *
     * \param globalPos The center of the finite volume which ought to be set.
     */
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    {
        return PrimaryVariables(1.0e5);
    }

    /*!
     * \brief Returns the temperature \f$\mathrm{[K]}\f$ for an isothermal problem.
     *
     * This is not specific to the discretization. By default it just
     * throws an exception so it must be overloaded by the problem if
     * no energy equation is used.
     */
    Scalar temperature() const
    {
        return 283.15; // 10°C
    }
};

} // end namespace Dumux

#endif
