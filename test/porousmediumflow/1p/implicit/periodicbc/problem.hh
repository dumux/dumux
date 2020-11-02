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

#ifndef DUMUX_INCOMPRESSIBLE_ONEP_TEST_PROBLEM_HH
#define DUMUX_INCOMPRESSIBLE_ONEP_TEST_PROBLEM_HH

#if HAVE_DUNE_SPGRID
#include <dune/grid/spgrid.hh>
#endif

#include <dumux/common/boundarytypes.hh>

#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/ccmpfa.hh>
#include <dumux/discretization/box.hh>

#include <dumux/porousmediumflow/problem.hh>
#include <dumux/porousmediumflow/1p/model.hh>
#include <dumux/porousmediumflow/1p/incompressiblelocalresidual.hh>

#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/1pliquid.hh>

#include "spatialparams.hh"

#ifndef FVGEOMCACHING
#define FVGEOMCACHING 0
#endif

namespace Dumux {
// forward declarations
template<class TypeTag> class OnePTestProblem;

namespace Properties {
// create the type tag nodes
// Create new type tags
namespace TTag {
struct OnePIncompressible { using InheritsFrom = std::tuple<OneP>; };
struct OnePIncompressibleTpfa { using InheritsFrom = std::tuple<OnePIncompressible, CCTpfaModel>; };
struct OnePIncompressibleMpfa { using InheritsFrom = std::tuple<OnePIncompressible, CCMpfaModel>; };
struct OnePIncompressibleBox { using InheritsFrom = std::tuple<OnePIncompressible, BoxModel>; };
} // end namespace TTag

// Set the grid type
#if HAVE_DUNE_SPGRID
template<class TypeTag>
struct Grid<TypeTag, TTag::OnePIncompressible> { using type = Dune::SPGrid<double, 2>; };
#endif

// Set the problem type
template<class TypeTag>
struct Problem<TypeTag, TTag::OnePIncompressible> { using type = OnePTestProblem<TypeTag>; };

// set the spatial params
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::OnePIncompressible> { using type = OnePTestSpatialParams<TypeTag>; };

// use the incompressible local residual (provides analytic jacobian)
template<class TypeTag>
struct LocalResidual<TypeTag, TTag::OnePIncompressible> { using type = OnePIncompressibleLocalResidual<TypeTag>; };

// the fluid system
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::OnePIncompressible>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = FluidSystems::OnePLiquid<Scalar, Components::SimpleH2O<Scalar> >;
};

// Enable caching
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::OnePIncompressible> { static constexpr bool value = false; };
template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::OnePIncompressible> { static constexpr bool value = false; };
template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::OnePIncompressible> { static constexpr bool value = FVGEOMCACHING; };
} // end namespace Properties

/*!
 * \ingroup OnePTests
 * \brief  Test problem for the incompressible one-phase model.
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
    using SourceValues = GetPropType<TypeTag, Properties::NumEqVector>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    static constexpr int dimWorld = GridView::dimensionworld;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    OnePTestProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        sourcePosition_ = getParam<GlobalPosition>("Source.Position");
        sourceValues_ = getParam<SourceValues>("Source.Values");
        sourceRadius_ = getParam<Scalar>("Source.Radius", 0.1);
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
        values.setAllDirichlet();
        return values;
    }

    /*!
     * \brief Applies a vector of point sources. The point sources
     *        are possibly solution dependent.
     *
     * \param pointSources A vector of PointSource s that contain
              source values for all phases and space positions.
     *
     * For this method, the values method of the point source
     * has to return the absolute rate values in units
     * \f$ [ \textnormal{unit of conserved quantity} / s ] \f$.
     * Positive values mean that the conserved quantity is created, negative ones mean that it vanishes.
     * E.g. for the mass balance that would be a mass rate in \f$ [ kg / s ] \f$.
     */
    template<class PointSource>
    void addPointSources(std::vector<PointSource>& pointSources) const
    {
        constexpr std::size_t numS = 20;
        pointSources.reserve(numS);
        for (int i = 0; i < numS; ++i)
        {
            const double angle = double(i)/numS*M_PI*2.0;
            const double weight = 1.0/double(numS);
            auto values = sourceValues_;
            values *= weight;
            auto pos = sourcePosition_;
            pos.axpy(sourceRadius_, GlobalPosition({std::cos(angle), std::sin(angle)}));
            pointSources.emplace_back(pos, values);
        }
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

private:
    GlobalPosition sourcePosition_;
    SourceValues sourceValues_;
    Scalar sourceRadius_;
};

} // end namespace Dumux

#endif
