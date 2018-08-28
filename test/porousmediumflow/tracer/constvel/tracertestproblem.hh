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
 * \ingroup TracerTests
 * \brief Definition of a problem, for the tracer problem:
 * A rotating velocity field mixes a tracer band in a porous groundwater reservoir.
 */
#ifndef DUMUX_TRACER_TEST_PROBLEM_HH
#define DUMUX_TRACER_TEST_PROBLEM_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/box/properties.hh>
#include <dumux/discretization/cellcentered/tpfa/properties.hh>
#include <dumux/discretization/cellcentered/mpfa/properties.hh>
#include <dumux/porousmediumflow/tracer/model.hh>
#include <dumux/porousmediumflow/problem.hh>
#include <dumux/material/fluidsystems/base.hh>

#include "tracertestspatialparams.hh"

#ifndef USEMOLES // default to true if not set through CMake
#define USEMOLES true
#endif

namespace Dumux {
/**
 * \ingroup TracerTests
 * \brief Definition of a problem, for the tracer problem:
 * A rotating velocity field mixes a tracer band in a porous groundwater reservoir.
 */
template <class TypeTag>
class TracerTest;

namespace Properties {
NEW_TYPE_TAG(TracerTest, INHERITS_FROM(Tracer));
NEW_TYPE_TAG(TracerTestTpfa, INHERITS_FROM(CCTpfaModel, TracerTest));
NEW_TYPE_TAG(TracerTestMpfa, INHERITS_FROM(CCMpfaModel, TracerTest));
NEW_TYPE_TAG(TracerTestBox, INHERITS_FROM(BoxModel, TracerTest));

// enable caching
SET_BOOL_PROP(TracerTest, EnableGridVolumeVariablesCache, true);
SET_BOOL_PROP(TracerTest, EnableGridFluxVariablesCache, true);
SET_BOOL_PROP(TracerTest, EnableFVGridGeometryCache, true);

// Set the grid type
SET_TYPE_PROP(TracerTest, Grid, Dune::YaspGrid<2>);

// Set the problem property
SET_TYPE_PROP(TracerTest, Problem, TracerTest<TypeTag>);

// Set the spatial parameters
SET_PROP(TracerTest, SpatialParams)
{
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using type = TracerTestSpatialParams<FVGridGeometry, Scalar>;
};

// Define whether mole(true) or mass (false) fractions are used
SET_BOOL_PROP(TracerTest, UseMoles, USEMOLES);

//! A simple fluid system with one tracer component
template<class TypeTag>
class TracerFluidSystem : public FluidSystems::Base<typename GET_PROP_TYPE(TypeTag, Scalar),
                                                               TracerFluidSystem<TypeTag>>
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;

public:
    static constexpr bool isTracerFluidSystem()
    { return true; }

    //! None of the components are the main component of the phase
    static constexpr int getMainComponent(int phaseIdx)
    { return -1; }

    //! The number of components
    static constexpr int numComponents = 2;
    static constexpr int numPhases = 1;

    //! Human readable component name (index compIdx) (for vtk output)
    static std::string componentName(int compIdx)
    { return "tracer_" + std::to_string(compIdx); }

    //! Human readable phase name (index phaseIdx) (for velocity vtk output)
    static std::string phaseName(int phaseIdx = 0)
    { return "Groundwater"; }

    //! Molar mass in kg/mol of the component with index compIdx
    static Scalar molarMass(unsigned int compIdx)
    { return 0.300; }

    //! binary diffusion coefficient
    //! (might depend on spatial parameters like pressure / temperature)
    static Scalar binaryDiffusionCoefficient(unsigned int compIdx,
                                             const Problem& problem,
                                             const Element& element,
                                             const SubControlVolume& scv)
    {
        static const Scalar D = getParam<Scalar>("Problem.D");
        static const Scalar D2 = getParam<Scalar>("Problem.D2");
        if (compIdx == 0)
            return D;
        else
            return D2;
    }

    /*!
     * \copydoc Base::isCompressible
     */
    static constexpr bool isCompressible(int phaseIdx)
    { return false; }

     /*!
     * \copydoc Base::viscosityIsConstant
     */
    static constexpr bool viscosityIsConstant(int phaseIdx)
    { return true; }
};

SET_TYPE_PROP(TracerTest, FluidSystem, TracerFluidSystem<TypeTag>);

} // end namespace Properties


/*!
 * \ingroup TracerModel
 * \ingroup ImplicitTestProblems
 *
 * \brief Definition of a problem, for the tracer problem:
 * A lens of contaminant tracer is diluted by diffusion and a base groundwater flow
 *
 * This problem uses the \ref TracerModel model.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_boxtracer -ParameterFile ./test_boxtracer.input</tt> or
 * <tt>./test_cctracer -ParameterFile ./test_cctracer.input</tt>
 */
template <class TypeTag>
class TracerTest : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Indices = typename GET_PROP_TYPE(TypeTag, ModelTraits)::Indices;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using SpatialParams = typename GET_PROP_TYPE(TypeTag, SpatialParams);

    //! property that defines whether mole or mass fractions are used
    static constexpr bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);
    using Element = typename FVGridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    TracerTest(std::shared_ptr<const FVGridGeometry> fvGridGeom)
    : ParentType(fvGridGeom)
    {
        // stating in the console whether mole or mass fractions are used
        if(useMoles)
            std::cout<<"problem uses mole fractions" << '\n';
        else
            std::cout<<"problem uses mass fractions" << '\n';
    }

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
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;
        values.setAllNeumann(); // no-flow
        return values;
    }
    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * \param globalPos The position for which the initial condition should be evaluated
     *
     * For this method, the \a values parameter stores primary
     * variables.
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables initialValues(0.0);
        if (globalPos[1] > 0.4 - eps_ && globalPos[1] < 0.6 + eps_)
        {
            if (useMoles)
                initialValues = 1e-9;
            else
                initialValues = 1e-9*FluidSystem::molarMass(0)/this->spatialParams().fluidMolarMass(globalPos);
        }
        return initialValues; }

    // \}

private:
    static constexpr Scalar eps_ = 1e-6;
};

} //end namespace Dumux

#endif
