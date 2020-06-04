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
/**
 * \file
 * \ingroup OnePNCTests
 * \brief Definition of a problem involving salt
 *        water intrusion into a fresh water aquifer.
 */

#ifndef DUMUX_SALTWATERINTRUSION_TEST_PROBLEM_HH
#define DUMUX_SALTWATERINTRUSION_TEST_PROBLEM_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/common/boundarytypes.hh>

#include <dumux/discretization/box.hh>
#include <dumux/porousmediumflow/1pnc/model.hh>
#include <dumux/porousmediumflow/problem.hh>

#include <dumux/material/fluidsystems/brine.hh>
#include "../../spatialparams.hh"

namespace Dumux {

template <class TypeTag>
class SaltWaterIntrusionTestProblem;

namespace Properties {
// Create new type tags
namespace TTag {
struct SaltWaterIntrusionTest { using InheritsFrom = std::tuple<OnePNC, BoxModel>; };
} // end namespace TTag

// Use a structured yasp grid
template<class TypeTag>
struct Grid<TypeTag, TTag::SaltWaterIntrusionTest> { using type = Dune::YaspGrid<2>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::SaltWaterIntrusionTest> { using type = SaltWaterIntrusionTestProblem<TypeTag>; };

// Set fluid configuration
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::SaltWaterIntrusionTest>
{ using type = FluidSystems::Brine< GetPropType<TypeTag, Properties::Scalar> >; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::SaltWaterIntrusionTest>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = OnePNCTestSpatialParams<GridGeometry, Scalar>;
};

// Use mass fractions to set salinity conveniently
template<class TypeTag>
struct UseMoles<TypeTag, TTag::SaltWaterIntrusionTest> { static constexpr bool value = false; };

} // end namespace Properties

/*!
 * \ingroup OnePNCTests
 * \brief Definition of a problem involving salt water intrusion into a
 *        fresh water aquifer.
 *
 * The aquifer is in contact with salt water with a salinity of 0.03 on the
 * right boundary.
 * \note To run the simulation execute the following line in shell:
 * <tt>./test_saltwaterintrusion -parameterFile saltwaterintrusion.input</tt> or
 */
template <class TypeTag>
class SaltWaterIntrusionTestProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;

    // copy pressure index for convenience
    enum { pressureIdx = Indices::pressureIdx };

    //! The test is defined using mass fractions
    static_assert(!getPropValue<TypeTag, Properties::UseMoles>(), "This test uses mass fractions!");

public:
    SaltWaterIntrusionTestProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        //initialize fluid system
        FluidSystem::init();
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief Returns the temperature within the domain [K].
     * This problem assumes a temperature of 20 degrees Celsius.
     */
    Scalar temperature() const
    { return 273.15 + 20; } // in [K]

    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param globalPos The position for which the bc type should be evaluated
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        values.setAllNeumann();

        // use Dirichlet BCs on the left and right boundary
        if(globalPos[0] < eps_ || globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_)
            values.setAllDirichlet();

        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet boundary segment.
     *
     * \param globalPos The position for which the bc type should be evaluated
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    {
        PrimaryVariables values = initialAtPos(globalPos);

        // salt water is in contact on the right boundary
        if (globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_)
            values[FluidSystem::NaClIdx] = 0.035; // 3.5% salinity (sea water)

        return values;
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluates the initial value for a control volume.
     *
     * \param globalPos The position for which the initial condition should be evaluated
     */
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    {
        const Scalar depth = this->gridGeometry().bBoxMax()[1] - globalPos[1];
        PrimaryVariables priVars;
        priVars[pressureIdx] = 1e5 + depth*9.81*1000; // hydrostatic pressure (assume rho_water = 1000.0)
        priVars[FluidSystem::NaClIdx] = 0.0;          // initially only fresh water is present
        return priVars;
    }

    // \}

private:
    static constexpr Scalar eps_ = 1e-6;
};

} // end namespace Dumux

#endif
