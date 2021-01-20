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
 * \brief Definition of a problem for a 1p3c problem:
 *        Component transport of N2, CO2 and H2 using the Maxwell-Stefan diffusion law.
 */

#ifndef DUMUX_1P3C_TEST_PROBLEM_HH
#define DUMUX_1P3C_TEST_PROBLEM_HH

#if HAVE_UG
#include <dune/grid/uggrid.hh>
#endif
#include <dune/grid/yaspgrid.hh>

#include <dumux/common/boundarytypes.hh>

#include <dumux/discretization/cctpfa.hh>
#include <dumux/discretization/box.hh>
#include <dumux/discretization/evalsolution.hh>
#include <dumux/discretization/evalgradients.hh>
#include <dumux/porousmediumflow/1pnc/model.hh>
#include <dumux/porousmediumflow/problem.hh>

#include <dumux/material/idealgas.hh>
#include <dumux/material/fluidsystems/base.hh>

#include "../1p2c/spatialparams.hh"

#include <dumux/io/gnuplotinterface.hh>
#include <dumux/flux/maxwellstefanslaw.hh>

namespace Dumux {

template <class TypeTag>
class MaxwellStefanOnePThreeCTestProblem;

namespace Properties {

// Create new type tags
namespace TTag {
struct MaxwellStefanOnePThreeCTest { using InheritsFrom = std::tuple<OnePNC>; };
struct MaxwellStefanOnePThreeCTestBox { using InheritsFrom = std::tuple<MaxwellStefanOnePThreeCTest, BoxModel>; };
struct MaxwellStefanOnePThreeCTestCCTpfa { using InheritsFrom = std::tuple<MaxwellStefanOnePThreeCTest, CCTpfaModel>; };
} // end namespace TTag

// Set the grid type
#if HAVE_UG
template<class TypeTag>
struct Grid<TypeTag, TTag::MaxwellStefanOnePThreeCTest> { using type = Dune::UGGrid<2>; };
#else
template<class TypeTag>
struct Grid<TypeTag, TTag::MaxwellStefanOnePThreeCTest> { using type = Dune::YaspGrid<2>; };
#endif

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::MaxwellStefanOnePThreeCTest> { using type = MaxwellStefanOnePThreeCTestProblem<TypeTag>; };

/*!
 * \ingroup OnePNCTests
 * \brief  A simple fluid system with three components for testing the  multi-component diffusion with the Maxwell-Stefan formulation.
 */
template<class TypeTag>
class H2N2CO2FluidSystem
: public FluidSystems::Base<GetPropType<TypeTag, Properties::Scalar>, H2N2CO2FluidSystem<TypeTag>>

{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ThisType = H2N2CO2FluidSystem<TypeTag>;
    using Base = FluidSystems::Base<Scalar, ThisType>;
    using IdealGas = Dumux::IdealGas<Scalar>;

public:
    //! The number of phases
    static constexpr int numPhases = 1;
    static constexpr int numComponents = 3;

    static constexpr int H2Idx = 0; //first major component
    static constexpr int N2Idx = 1; //second major component
    static constexpr int CO2Idx = 2; //secondary component

    //! Human readable component name (index compIdx) (for vtk output)
    static std::string componentName(int compIdx)
    { return "MaxwellStefan_" + std::to_string(compIdx); }

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
        DUNE_THROW(Dune::InvalidStateException, "Invalid component index " << compIdx);;
    }

    using Base::binaryDiffusionCoefficient;
    /*!
     * \brief Given a phase's composition, temperature and pressure,
     *        returns the binary diffusion coefficient \f$\mathrm{[m^2/s]}\f$ for components
     *        \f$i\f$ and \f$j\f$ in this phase.
     *
     * \param fluidState An arbitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     * \param compIIdx The index of the first component to consider
     * \param compJIdx The index of the second component to consider
     */
    template <class FluidState>
    static Scalar binaryDiffusionCoefficient(const FluidState &fluidState,
                                             int phaseIdx,
                                             int compIIdx,
                                             int compJIdx)
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
    using Base::density;
    /*!
     * \brief Given a phase's composition, temperature, pressure, and
     *        the partial pressures of all components, returns its
     *        density \f$\mathrm{[kg/m^3]}\f$.
     *
     * \param phaseIdx index of the phase
     * \param fluidState the fluid state
     */
    template <class FluidState>
    static Scalar density(const FluidState &fluidState,
                          const int phaseIdx)
    {
        Scalar T = fluidState.temperature(phaseIdx);
        Scalar p = fluidState.pressure(phaseIdx);
        return IdealGas::molarDensity(T, p) * fluidState.averageMolarMass(0);
    }

    using Base::viscosity;
    /*!
     * \brief Calculates the dynamic viscosity of a fluid phase \f$\mathrm{[Pa*s]}\f$
     *
     * \param fluidState An arbitrary fluid state
     * \param phaseIdx The index of the fluid phase to consider
     */
    template <class FluidState>
    static Scalar viscosity(const FluidState &fluidState,
                            int phaseIdx)
    {
        return 1e-6;
    }

    using Base::molarDensity;
    /*!
     * \brief The molar density \f$\rho_{mol,\alpha}\f$
     *   of a fluid phase \f$\alpha\f$ in \f$\mathrm{[mol/m^3]}\f$
     *
     * The molar density for the simple relation is defined by the
     * mass density \f$\rho_\alpha\f$ and the molar mass of the main component \f$M_\kappa\f$:
     *
     * \f[\rho_{mol,\alpha} = \frac{\rho_\alpha}{M_\kappa} \;.\f]
     */
    template <class FluidState>
    static Scalar molarDensity(const FluidState &fluidState, int phaseIdx)
    {
        Scalar T = fluidState.temperature(phaseIdx);
        Scalar p = fluidState.pressure(phaseIdx);
        return IdealGas::molarDensity(T,p);
    }
};

// Set fluid configuration
template<class TypeTag>
struct FluidSystem<TypeTag, TTag::MaxwellStefanOnePThreeCTest>
{using type = H2N2CO2FluidSystem<TypeTag>; };

// Set the spatial parameters
template<class TypeTag>
struct SpatialParams<TypeTag, TTag::MaxwellStefanOnePThreeCTest>
{
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using type = OnePNCTestSpatialParams<GridGeometry, Scalar>;
};

// Define whether mole(true) or mass (false) fractions are used
template<class TypeTag>
struct UseMoles<TypeTag, TTag::MaxwellStefanOnePThreeCTest> { static constexpr bool value = true; };

//! Here we set FicksLaw or MaxwellStefansLaw
template<class TypeTag>
struct MolecularDiffusionType<TypeTag, TTag::MaxwellStefanOnePThreeCTest> { using type = MaxwellStefansLaw<TypeTag>; };
}

/*!
 * \ingroup OnePNCTests
 * \brief Definition of a problem for a 1p3c problem:
 *        Component transport of N2, CO2 and H2.

 * The domain is closed on all sides. H2 constitutes the bulk gas phase.
 * Initially, there is N2 and CO2 in the left half of the domain,
 * while only N2 is present in the right half of the domain.
 * Over time, the concentrations will equilibrate.
 *
 * This problem uses the \ref OnePNCModel model and the Maxwell-Stefan diffusion law.
 *
 */
template <class TypeTag>
class MaxwellStefanOnePThreeCTestProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;

    //! property that defines whether mole or mass fractions are used
    static constexpr bool useMoles = getPropValue<TypeTag, Properties::UseMoles>();

    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    MaxwellStefanOnePThreeCTestProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        name_ = getParam<std::string>("Problem.Name");

        // stating in the terminal whether mole or mass fractions are used
        if (useMoles)
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


    /*!
     * \brief Returns the temperature within the domain [K].
     *
     * This problem assumes a temperature of 20 degrees Celsius.
     */
    Scalar temperature() const
    { return 273.15 + 20; } // in [K]

    //! Called after every time step
    void plotComponentsOverTime(const SolutionVector& curSol, Scalar time)
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
                for (const auto& element : elements(this->gridGeometry().gridView()))
                {
                    auto fvGeometry = localView(this->gridGeometry());
                    fvGeometry.bindElement(element);

                    const auto elemSol = elementSolution(element, curSol, this->gridGeometry());
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

                // do a gnuplot
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
                gnuplot2_.addDataSetToPlot(x_, y4_, "CO2_right.dat", "w l t 'CO_2 right");
                gnuplot2_.plot("mole_fraction_CO2");

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
     * \brief Evaluates the boundary conditions for a Neumann boundary segment.
     *
     * \param globalPos The position for which the bc type should be evaluated
     */
    NumEqVector neumannAtPos(const GlobalPosition& globalPos) const
    { return NumEqVector(0.0); }

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
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables initialValues(0.0);
        if (globalPos[0] < 0.5)
        {
           initialValues[Indices::pressureIdx] = 1e5;
           initialValues[FluidSystem::N2Idx] = 0.50086;
           initialValues[FluidSystem::CO2Idx] = 0.49914;
        }
        else
        {
           initialValues[Indices::pressureIdx] = 1e5;
           initialValues[FluidSystem::N2Idx] = 0.49879;
           initialValues[FluidSystem::CO2Idx] = 0.0;
        }
        return initialValues;
    }

    // \}

private:
    std::string name_;

    Dumux::GnuplotInterface<Scalar> gnuplot_;
    Dumux::GnuplotInterface<Scalar> gnuplot2_;
    Dumux::GnuplotInterface<Scalar> gnuplot3_;

    std::vector<Scalar> x_;
    std::vector<Scalar> y_;
    std::vector<Scalar> y2_;
    std::vector<Scalar> y3_;
    std::vector<Scalar> y4_;
    std::vector<Scalar> y5_;
    std::vector<Scalar> y6_;

    bool plotOutput_;
};

} // end namespace Dumux

#endif
