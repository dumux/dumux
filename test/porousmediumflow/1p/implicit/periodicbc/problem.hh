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
#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#endif

#include <dumux/discretization/cellcentered/tpfa/properties.hh>
#include <dumux/discretization/cellcentered/mpfa/properties.hh>
#include <dumux/discretization/box/properties.hh>

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
SET_TYPE_PROP(OnePIncompressible, Grid, Dune::SPGrid<double, 2>);
// SET_TYPE_PROP(OnePIncompressible, Grid, Dune::ALUGrid<2, 2, Dune::cube, Dune::nonconforming>);

// Set the problem type
SET_TYPE_PROP(OnePIncompressible, Problem, OnePTestProblem<TypeTag>);

// set the spatial params
SET_TYPE_PROP(OnePIncompressible, SpatialParams, OnePTestSpatialParams<TypeTag>);

// use the incompressible local residual (provides analytic jacobian)
SET_TYPE_PROP(OnePIncompressible, LocalResidual, OnePIncompressibleLocalResidual<TypeTag>);

// the fluid system
SET_PROP(OnePIncompressible, FluidSystem)
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using type = FluidSystems::OnePLiquid<Scalar, Components::SimpleH2O<Scalar> >;
};

// Enable caching
SET_BOOL_PROP(OnePIncompressible, EnableGridVolumeVariablesCache, false);
SET_BOOL_PROP(OnePIncompressible, EnableGridFluxVariablesCache, false);
SET_BOOL_PROP(OnePIncompressible, EnableFVGridGeometryCache, FVGEOMCACHING);
} // end namespace Properties

/*!
 * \ingroup OnePTests
 * \brief  Test problem for the incompressible one-phase model:
 *\todo doc me!
 * <tt>./test_box1pfv</tt> or
 * <tt>./test_cc1pfv</tt>
 */
template<class TypeTag>
class OnePTestProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using SourceValues = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    static constexpr int dimWorld = GridView::dimensionworld;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    OnePTestProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
    {
        sourcePosition_ = getParam<GlobalPosition>("Source.Position");
        sourceValues_ = getParam<SourceValues>("Source.Values");
        sourceRadius_ = getParam<Scalar>("Source.Radius", 0.1);
    }

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary control volume.
     *
     * \param values The boundary types for the conservation equations
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
     * \brief Evaluate the boundary conditions for a Dirichlet
     *        control volume.
     *
     * \param values The Dirichlet values for the primary variables
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
