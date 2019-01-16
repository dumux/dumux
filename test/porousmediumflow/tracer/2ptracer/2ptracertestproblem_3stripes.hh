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
#ifndef DUMUX_TWOP_TRACER_TEST_PROBLEM_HH
#define DUMUX_TWOP_TRACER_TEST_PROBLEM_HH

#include <dumux/discretization/cellcentered/tpfa/properties.hh>
#include <dumux/porousmediumflow/tracer/model.hh>
#include <dumux/porousmediumflow/problem.hh>
#include <dumux/material/fluidsystems/base.hh>

#include "2ptracertestspatialparams.hh"

namespace Dumux
{
/**
 * \ingroup TracerTests
 * \brief Definition of a problem, for the tracer problem:
 * A rotating velocity field mixes a tracer band in a porous groundwater reservoir.
 */
template <class TypeTag>
class TwoPTracerTestProblem;

namespace Properties
{
NEW_TYPE_TAG(TwoPTracerTestTypeTag, INHERITS_FROM(Tracer));
NEW_TYPE_TAG(TwoPTracerTestCCTypeTag, INHERITS_FROM(CCTpfaModel, TwoPTracerTestTypeTag));

// enable caching
SET_BOOL_PROP(TwoPTracerTestTypeTag, EnableGridVolumeVariablesCache, false);
SET_BOOL_PROP(TwoPTracerTestTypeTag, EnableGridFluxVariablesCache, false);
SET_BOOL_PROP(TwoPTracerTestTypeTag, EnableFVGridGeometryCache, false);

// Set the grid type
SET_TYPE_PROP(TwoPTracerTestTypeTag, Grid, Dune::YaspGrid<2>);

// Set the problem property
SET_TYPE_PROP(TwoPTracerTestTypeTag, Problem, TwoPTracerTestProblem<TypeTag>);

// Set the spatial parameters
SET_PROP(TwoPTracerTestTypeTag, SpatialParams)
{
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using type = TwoPTracerTestSpatialParams<FVGridGeometry, Scalar>;
};
// Define whether mole(true) or mass (false) fractions are used
SET_BOOL_PROP(TwoPTracerTestTypeTag, UseMoles, false);
SET_BOOL_PROP(TwoPTracerTestCCTypeTag, SolutionDependentMolecularDiffusion, false);

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
    //! If the fluid system only contains tracer components
    static constexpr bool isTracerFluidSystem()
    { return true; }

    //! No component is the main component
    static constexpr int getMainComponent(int phaseIdx)
    { return -1; }

    //! The number of components
    static constexpr int numComponents = 1;

    //! Human readable component name (index compIdx) (for vtk output)
    static std::string componentName(int compIdx)
    { return "tracer_" + std::to_string(compIdx); }

    //! Human readable phase name (index phaseIdx) (for velocity vtk output)
    //TODO adapt for multiphase and accordingly for the vtkOutput
    static std::string phaseName(int phaseIdx = 0)
    { if (phaseIdx == 0)
        return "Water";
      else
        return "NotWater";
    }

    //! Molar mass in kg/mol of the component with index compIdx
    static Scalar molarMass(unsigned int compIdx)
    { return 0.300; }

    //! binary diffusion coefficient
    //! (might depend on spatial parameters like pressure / temperature)
    static Scalar binaryDiffusionCoefficient(unsigned int compIdx,
                                             const Problem& problem,
                                             const Element& element,
                                             const SubControlVolume& scv)
    { return 0.0; }
};

SET_TYPE_PROP(TwoPTracerTestTypeTag, FluidSystem, TracerFluidSystem<TypeTag>);

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
class TwoPTracerTestProblem : public PorousMediumFlowProblem<TypeTag>
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

    static const int dimWorld = GridView::dimensionworld;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

public:
    TwoPTracerTestProblem(std::shared_ptr<const FVGridGeometry> fvGridGeom)
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
        values.setAllNeumann();
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

        // VERSION "3stripes"
        if (onUpperBoundary_(globalPos) ||  onStripe1_(globalPos) || onStripe2_(globalPos) || onStripe3_(globalPos))
        {
            if (useMoles)
                initialValues = 1e-9;
            else
                initialValues = 1e-9*FluidSystem::molarMass(0)/this->spatialParams().fluidMolarMass(globalPos);
        }
        return initialValues;
    }


private:


    static constexpr Scalar eps_ = 1e-6;

    Scalar yMax_ = this->fvGridGeometry().bBoxMax()[1];
    Scalar xMax_ = this->fvGridGeometry().bBoxMax()[0];

// TODO allgemeiner Abruf der Anzahl von Zellen in Y-Richtung/ X-Richtung
    Scalar yNumCells_ = 32;
    Scalar xNumCells_ = 48;
    Scalar cellHeight_ = yMax_/yNumCells_;
    Scalar cellWidth_ = xMax_/xNumCells_;
    Scalar width_ = xMax_ - this->fvGridGeometry().bBoxMin()[0];

    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[1] > yMax_ - 0.1 - eps_;
    }

    // VERSION "3stripes"
    bool onStripe1_(const GlobalPosition &globalPos) const
    {
       return  (
           ( (yMax_ /4.0 - cellHeight_*0.5) <= globalPos[1] ) &&
           ( (yMax_/4.0 + cellHeight_*0.5) > globalPos[1] )
       );
    }

    bool onStripe2_(const GlobalPosition &globalPos) const
    {
        return  (
            ( (2.0 * yMax_ /4.0 - cellHeight_*0.5) <= globalPos[1] ) &&
            ( (2.0 * yMax_/4.0 + cellHeight_*0.5) > globalPos[1] )
        );
    }

    bool onStripe3_(const GlobalPosition &globalPos) const
    {
        return  (
            ( (3.0 * yMax_ /4.0 - cellHeight_*0.5) <= globalPos[1] ) &&
            ( (3.0 * yMax_/4.0 + cellHeight_*0.5) > globalPos[1] )
        );
    }

};

} //end namespace Dumux

#endif
