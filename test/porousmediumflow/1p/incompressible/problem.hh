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

#if HAVE_UG
#include <dune/grid/uggrid.hh>
#endif
#include <dune/grid/yaspgrid.hh>
#if HAVE_QUADMATH
#include <dune/common/quadmath.hh>
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

#ifndef GRIDTYPE // default to yasp grid if not provided by CMake
#define GRIDTYPE Dune::YaspGrid<2>
#endif

namespace Dumux {
// forward declarations
template<class TypeTag> class OnePTestProblem;

namespace Properties {
// Create new type tags
namespace TTag {
struct OnePIncompressible { using InheritsFrom = std::tuple<OneP>; };
struct OnePIncompressibleTpfa { using InheritsFrom = std::tuple<OnePIncompressible, CCTpfaModel>; };
struct OnePIncompressibleMpfa { using InheritsFrom = std::tuple<OnePIncompressible, CCMpfaModel>; };
struct OnePIncompressibleBox { using InheritsFrom = std::tuple<OnePIncompressible, BoxModel>; };
} // end namespace TTag

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::OnePIncompressible> { using type = GRIDTYPE; };

// Set the problem type
template<class TypeTag>
struct Problem<TypeTag, TTag::OnePIncompressible> { using type = OnePTestProblem<TypeTag>; };

// set the spatial params
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::OnePIncompressible>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = OnePTestSpatialParams<GridGeometry, Scalar>;
};

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
struct EnableGridGeometryCache<TypeTag, TTag::OnePIncompressible> { static constexpr bool value = false; };

// define a TypeTag for a quad precision test
#if HAVE_QUADMATH
namespace TTag {
struct OnePIncompressibleTpfaQuad { using InheritsFrom = std::tuple<OnePIncompressibleTpfa>; };
} // end namespace TTag
template<class TypeTag>
struct Scalar<TypeTag, TTag::OnePIncompressibleTpfaQuad> { using type = Dune::Float128; };
#endif
} // end namespace Properties

/*!
 * \ingroup OnePTests
 * \brief  Test problem for the incompressible one-phase model
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
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    static constexpr int dimWorld = GridView::dimensionworld;
    using GlobalPosition = Dune::FieldVector<Scalar,dimWorld>;

public:
    OnePTestProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry), velocity_(0.0)
    {
        extrusionFactor_ = getParam<Scalar>("Problem.ExtrusionFactor");
        Scalar permeability = getParam<Scalar>("SpatialParams.Permeability");
        dp_dy_ = -1.0e+5;

        const bool checkIsConstantVelocity = getParam<bool>("Problem.CheckIsConstantVelocity", false);
        if(checkIsConstantVelocity)
        {
            velocity_[dimWorld-1] = -permeability * dp_dy_;
            velocity_[dimWorld-1] /= FluidSystem::viscosity(temperature(), Scalar(1.0e5));
        }
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
        values[0] = 1.0e+5 + dp_dy_*(globalPos[dimWorld-1] - this->gridGeometry().bBoxMax()[dimWorld-1]);

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

    /*!
     * \brief Returns how much the domain is extruded at a given position.
     *
     * This means the factor by which a lower-dimensional
     * entity needs to be expanded to get a full-dimensional cell.
     */
    Scalar extrusionFactorAtPos(const GlobalPosition &globalPos) const
    {
        return extrusionFactor_;
    }

    /*!
     * \brief Returns the velocity
     *
     * The velocity is given for the case of a linear pressure solution
     * with constant permeablity and without gravity.
     */
    const GlobalPosition velocity() const
    {
        return velocity_;
    }

private:
    Scalar extrusionFactor_;
    Scalar dp_dy_;
    GlobalPosition velocity_;
};

} // end namespace Dumux

#endif
