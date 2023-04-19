// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/**
 * \file
 * \ingroup OnePNCTests
 * \brief Definition of a problem involving salt
 *        water intrusion into a fresh water aquifer.
 */

#ifndef DUMUX_SALTWATERINTRUSION_TEST_PROBLEM_HH
#define DUMUX_SALTWATERINTRUSION_TEST_PROBLEM_HH

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>

#include <dumux/common/boundarytypes.hh>
#include <dumux/porousmediumflow/problem.hh>

namespace Dumux {

/*!
 * \ingroup OnePNCTests
 * \brief Definition of a problem involving salt water intrusion into a
 *        fresh water aquifer.
 *
 * The aquifer is in contact with salt water with a specific salinity on the
 * right boundary.
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
