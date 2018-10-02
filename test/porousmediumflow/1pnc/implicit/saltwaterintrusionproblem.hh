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
/**
 * \file
 * \ingroup OnePNCTests
 * \brief Definition of a problem involving salt
 *        water intrusion into a fresh water aquifer.
 */
#ifndef DUMUX_SALTWATERINTRUSION_TEST_PROBLEM_HH
#define DUMUX_SALTWATERINTRUSION_TEST_PROBLEM_HH

#include <dumux/discretization/box/properties.hh>
#include <dumux/porousmediumflow/1pnc/model.hh>
#include <dumux/porousmediumflow/problem.hh>

#include <dumux/material/fluidsystems/brine.hh>
#include "1pnctestspatialparams.hh"

namespace Dumux {

template <class TypeTag>
class SaltWaterIntrusionTestProblem;

namespace Properties {
NEW_TYPE_TAG(SaltWaterIntrusionTestTypeTag, INHERITS_FROM(BoxModel, OnePNC));

// Use a structured yasp grid
SET_TYPE_PROP(SaltWaterIntrusionTestTypeTag, Grid, Dune::YaspGrid<2>);

// Set the problem property
SET_TYPE_PROP(SaltWaterIntrusionTestTypeTag, Problem, SaltWaterIntrusionTestProblem<TypeTag>);

// Set fluid configuration
SET_TYPE_PROP(SaltWaterIntrusionTestTypeTag,
              FluidSystem,
              FluidSystems::Brine< typename GET_PROP_TYPE(TypeTag, Scalar) >);

// Set the spatial parameters
SET_PROP(SaltWaterIntrusionTestTypeTag, SpatialParams)
{
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using type = OnePNCTestSpatialParams<FVGridGeometry, Scalar>;
};

// Use mass fractions to set salinity conveniently
SET_BOOL_PROP(SaltWaterIntrusionTestTypeTag, UseMoles, false);

} // end namespace Properties

/*!
 * \ingroup OnePNCTests
 * \brief Definition of a problem involving salt water intrusion into a
 *        fresh water aquifer. The aquifer is in contact with salt water
 *        with a salinity of 0.03 on the right boundary.
 * \note To run the simulation execute the following line in shell:
 * <tt>./test_saltwaterintrusion -parameterFile saltwaterintrusion.input</tt> or
 */
template <class TypeTag>
class SaltWaterIntrusionTestProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Indices = typename GET_PROP_TYPE(TypeTag, ModelTraits)::Indices;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);

    // copy pressure index for convenience
    enum { pressureIdx = Indices::pressureIdx };

    //! The test is defined using mass fractions
    static_assert(!GET_PROP_VALUE(TypeTag, UseMoles), "This test uses mass fractions!");

public:
    SaltWaterIntrusionTestProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
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
     * \param globalPos The position for which the bc type should be evaluated
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        values.setAllNeumann();

        // use Dirichlet BCs on the left and right boundary
        if(globalPos[0] < eps_ || globalPos[0] > this->fvGridGeometry().bBoxMax()[0] - eps_)
            values.setAllDirichlet();

        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet boundary segment.
     * \param globalPos The position for which the bc type should be evaluated
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    {
        PrimaryVariables values = initialAtPos(globalPos);

        // salt water is in contact on the right boundary
        if (globalPos[0] > this->fvGridGeometry().bBoxMax()[0] - eps_)
            values[FluidSystem::NaClIdx] = 0.03;

        return values;
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluate the initial value for a control volume.
     * \param globalPos The position for which the initial condition should be evaluated
     */
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    {
        const Scalar depth = this->fvGridGeometry().bBoxMax()[1] - globalPos[1];
        PrimaryVariables priVars;
        priVars[pressureIdx] = 1e5 + depth*9.81*1000; // hydrostatic pressure (assume rho_water = 1000.0)
        priVars[FluidSystem::NaClIdx] = 0.0;          // initially only fresh water is present
        return priVars;
    }

    // \}

private:
    static constexpr Scalar eps_ = 1e-6;
};

} //end namespace Dumux

#endif
