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
 * \ingroup TracerTest
 * \brief Definition of a problem, for the MaxwellStefan problem:
 * A rotating velocity field mixes a MaxwellStefan band in a porous groundwater reservoir.
 */
#ifndef DUMUX_MAXWELL_STEFAN_TEST_PROBLEM_HH
#define DUMUX_MAXWELL_STEFAN_TEST_PROBLEM_HH

#include <dumux/discretization/elementsolution.hh>
#include <dumux/discretization/box/properties.hh>
#include <dumux/discretization/cellcentered/tpfa/properties.hh>

#include <dumux/porousmediumflow/tracer/model.hh>
#include <dumux/porousmediumflow/problem.hh>

#include "maxwellstefantestspatialparams.hh"
#include <dumux/discretization/maxwellstefanslaw.hh>

#include <dumux/io/gnuplotinterface.hh>
#include <dumux/material/fluidsystems/base.hh>

namespace Dumux {
template <class TypeTag>
class MaxwellStefanTestProblem;

namespace Properties {
NEW_TYPE_TAG(MaxwellStefanTestTypeTag, INHERITS_FROM(Tracer));
NEW_TYPE_TAG(MaxwellStefanTestCCTypeTag, INHERITS_FROM(CCTpfaModel, MaxwellStefanTestTypeTag));
NEW_TYPE_TAG(MaxwellStefanTestBoxTypeTag, INHERITS_FROM(BoxModel, MaxwellStefanTestTypeTag));

// Set the grid type
SET_TYPE_PROP(MaxwellStefanTestTypeTag, Grid, Dune::YaspGrid<2>);

// Set the problem property
SET_TYPE_PROP(MaxwellStefanTestTypeTag, Problem, MaxwellStefanTestProblem<TypeTag>);

// Set the spatial parameters
SET_TYPE_PROP(MaxwellStefanTestTypeTag, SpatialParams, MaxwellStefanTestSpatialParams<TypeTag>);

// Define whether mole(true) or mass (false) fractions are used
SET_BOOL_PROP(MaxwellStefanTestTypeTag, UseMoles, true);

//! Here we set FicksLaw or MaxwellStefansLaw
SET_TYPE_PROP(MaxwellStefanTestTypeTag, MolecularDiffusionType, MaxwellStefansLaw<TypeTag>);

//! A simple fluid system with one MaxwellStefan component
template<class TypeTag>
class H2N2CO2FluidSystem
: public FluidSystems::BaseFluidSystem<typename GET_PROP_TYPE(TypeTag, Scalar), H2N2CO2FluidSystem<TypeTag>>

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
    //! The number of components
    static constexpr int numComponents = 3;

    static constexpr int H2Idx = 0;//first major component
    static constexpr int N2Idx = 1;//second major component
    static constexpr int CO2Idx = 2;//secondary component

    //! Human readable component name (index compIdx) (for vtk output)
    static std::string componentName(int compIdx)
    {
        switch (compIdx)
        {
        case H2Idx: return "H2";
        case N2Idx: return "N2";
        case CO2Idx:return "CO2";
        }
        DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << compIdx);
    }

    //! Human readable phase name (index phaseIdx) (for velocity vtk output)
    static std::string phaseName(int phaseIdx = 0)
    { return "Gas"; }

    //! Molar mass in kg/mol of the component with index compIdx
    static Scalar molarMass(unsigned int compIdx)
    {
        switch (compIdx)
        {
        case H2Idx: return 0.002;
        case N2Idx: return 0.028;
        case CO2Idx:return 0.044;
        }
        DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << compIdx);
    }

    //! binary diffusion coefficient
    //! (might depend on spatial parameters like pressure / temperature)
    static Scalar binaryDiffusionCoefficient(unsigned int compIdx,
                                             const Problem& problem,
                                             const Element& element,
                                             const SubControlVolume& scv)
    {
      if (compIdx == H2Idx)
          return 0;
      if (compIdx == N2Idx)
          return 83.3e-6;
      if (compIdx == CO2Idx)
          return 68.0e-6;
       DUNE_THROW(Dune::InvalidStateException,
                       "Binary diffusion coefficient of component "
                       << compIdx <<" is undefined!\n");
    }

    //! binary diffusion coefficient
    //! (might depend on spatial parameters like pressure / temperature)
    static Scalar binaryDiffusionCoefficient(unsigned int compIIdx,
                                             unsigned int compJIdx,
                                             const Problem& problem,
                                             const Element& element,
                                             const SubControlVolume& scv)
    {
        if (compIIdx > compJIdx)
        {
            using std::swap;
            swap(compIIdx, compJIdx);
        }

        if (compIIdx == H2Idx && compJIdx == N2Idx)
            return 83.3e-6;
        if (compIIdx == H2Idx && compJIdx == CO2Idx)
            return 68.0e-6;
        if (compIIdx == N2Idx && compJIdx == CO2Idx)
            return 16.8e-6;
        DUNE_THROW(Dune::InvalidStateException,
                       "Binary diffusion coefficient of components "
                       << compIIdx << " and " << compJIdx << " is undefined!\n");
    }
};

SET_TYPE_PROP(MaxwellStefanTestTypeTag, FluidSystem, H2N2CO2FluidSystem<TypeTag>);

} // end namespace Properties

/*!
 * \ingroup TracerTest
 * \brief Definition of a problem, for the MaxwellStefan problem *
 * This problem uses the \ref MaxwellStefanModel model.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_boxmaxwellstefan -ParameterFile ./test_boxmaxwellstefan.input</tt> or
 * <tt>./test_ccmaxwellstefan -ParameterFile ./test_ccMaxwellstefan.input</tt>
 */
template <class TypeTag>
class MaxwellStefanTestProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Indices = typename GET_PROP_TYPE(TypeTag, ModelTraits)::Indices;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using NumEqVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using SpatialParams = typename GET_PROP_TYPE(TypeTag, SpatialParams);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);

    //! property that defines whether mole or mass fractions are used
    static constexpr bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);

    using Element = typename FVGridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    MaxwellStefanTestProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
    {
        name_ = getParam<std::string>("Problem.Name");

        // stating in the console whether mole or mass fractions are used
        if(useMoles)
            std::cout<<"problem uses mole fractions" << '\n';
        else
            std::cout<<"problem uses mass fractions" << '\n';

        plotOutput_ = false;
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief The problem name.
     * This is used as a prefix for files generated by the simulation.
     */
    const std::string& name() const
    { return name_; }

    //! Called after every time step
    void postTimeStep(const SolutionVector& curSol, Scalar time)
    {
        if (plotOutput_)
        {
            Scalar x_co2_left = 0.0;
            Scalar x_n2_left = 0.0;
            Scalar x_co2_right = 0.0;
            Scalar x_n2_right = 0.0;
            Scalar x_h2_left = 0.0;
            Scalar x_h2_right = 0.0;
            Scalar i = 0.0;
            Scalar j = 0.0;
            if (!(time < 0.0))
            {
                for (const auto& element : elements(this->fvGridGeometry().gridView()))
                {
                    auto fvGeometry = localView(this->fvGridGeometry());
                    fvGeometry.bindElement(element);

                    const auto elemSol = elementSolution(element, curSol, this->fvGridGeometry());
                    for (auto&& scv : scvs(fvGeometry))
                    {
                        const auto& globalPos = scv.dofPosition();
                        VolumeVariables volVars;
                        volVars.update(elemSol, *this, element, scv);

                        if (globalPos[0] < 0.5)
                        {
                            x_co2_left += volVars.moleFraction(0,2);

                            x_n2_left += volVars.moleFraction(0,1);
                            x_h2_left += volVars.moleFraction(0,0);
                            i +=1;
                        }
                        else
                        {
                            x_co2_right += volVars.moleFraction(0,2);
                            x_n2_right += volVars.moleFraction(0,1);
                            x_h2_right += volVars.moleFraction(0,0);
                            j +=1;
                        }

                    }
                }
                x_co2_left /= i;
                x_n2_left /= i;
                x_h2_left /= i;
                x_co2_right /= j;
                x_n2_right /= j;
                x_h2_right /= j;

                //do a gnuplot
                x_.push_back(time); // in seconds
                y_.push_back(x_n2_left);
                y2_.push_back(x_n2_right);
                y3_.push_back(x_co2_left);
                y4_.push_back(x_co2_right);
                y5_.push_back(x_h2_left);
                y6_.push_back(x_h2_right);

                gnuplot_.resetPlot();
                gnuplot_.setXRange(0, std::min(time, 72000.0));
                gnuplot_.setYRange(0.4, 0.6);
                gnuplot_.setXlabel("time [s]");
                gnuplot_.setYlabel("mole fraction mol/mol");
                gnuplot_.addDataSetToPlot(x_, y_, "N2_left.dat", "w l t 'N_2 left'");
                gnuplot_.addDataSetToPlot(x_, y2_, "N2_right.dat", "w l t 'N_2 right'");
                gnuplot_.plot("mole_fraction_N2");

                gnuplot2_.resetPlot();
                gnuplot2_.setXRange(0, std::min(time, 72000.0));
                gnuplot2_.setYRange(0.0, 0.6);
                gnuplot2_.setXlabel("time [s]");
                gnuplot2_.setYlabel("mole fraction mol/mol");
                gnuplot2_.addDataSetToPlot(x_, y3_, "CO2_left.dat", "w l t 'CO_2 left'");
                gnuplot2_.addDataSetToPlot(x_, y4_, "C02_right.dat", "w l t CO_2 right");
                gnuplot2_.plot("mole_fraction_C02");

                gnuplot3_.resetPlot();
                gnuplot3_.setXRange(0, std::min(time, 72000.0));
                gnuplot3_.setYRange(0.0, 0.6);
                gnuplot3_.setXlabel("time [s]");
                gnuplot3_.setYlabel("mole fraction mol/mol");
                gnuplot3_.addDataSetToPlot(x_, y5_, "H2_left.dat", "w l t 'H_2 left'");
                gnuplot3_.addDataSetToPlot(x_, y6_, "H2_right.dat", "w l t 'H_2 right'");
                gnuplot3_.plot("mole_fraction_H2");
           }
        }
    }

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
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;
        values.setAllNeumann();
        return values;
    }

    /*!
     * \brief Evaluate the boundary conditions for a Neumann
     *        boundary segment.
     *
     * \param globalPos The position for which the bc type should be evaluated
     * The units must be according to either using mole or mass fractions. (mole/(m^2*s) or kg/(m^2*s))
     */
    NumEqVector neumannAtPos(const GlobalPosition& globalPos) const
    { return NumEqVector(0.0); }

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
        if (globalPos[0] < 0.5)
        {
           initialValues[FluidSystem::H2Idx] = 0.0;
           initialValues[FluidSystem::N2Idx] = 0.50086;
           initialValues[FluidSystem::CO2Idx] = 0.49914;
        }
        else
        {
           initialValues[FluidSystem::H2Idx] = 0.50121;
           initialValues[FluidSystem::N2Idx] = 0.49879;
           initialValues[FluidSystem::CO2Idx] = 0.0;
        }
        return initialValues;
    }

    // \}

private:
    static constexpr Scalar eps_ = 1e-6;
    std::string name_;

    Dumux::GnuplotInterface<double> gnuplot_;
    Dumux::GnuplotInterface<double> gnuplot2_;
    Dumux::GnuplotInterface<double> gnuplot3_;

    std::vector<double> x_;
    std::vector<double> y_;
    std::vector<double> y2_;
    std::vector<double> y3_;
    std::vector<double> y4_;
    std::vector<double> y5_;
    std::vector<double> y6_;

    bool plotOutput_;
};

} //end namespace Dumux

#endif
