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
 * \ingroup NavierStokesNCTests
 * \brief Channel flow test for the multi-component staggered grid (Navier-)Stokes model.
 */

#ifndef DUMUX_CHANNEL_MAXWELL_STEFAN_TEST_PROBLEM_HH
#define DUMUX_CHANNEL_MAXWELL_STEFAN_TEST_PROBLEM_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/staggered/freeflow/properties.hh>

#include <dumux/flux/maxwellstefanslaw.hh>

#include <dumux/freeflow/compositional/navierstokesncmodel.hh>
#include <dumux/freeflow/navierstokes/boundarytypes.hh>
#include <dumux/freeflow/navierstokes/problem.hh>

#include <dumux/io/gnuplotinterface.hh>

#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/h2oair.hh>


namespace Dumux {
template <class TypeTag>
class MaxwellStefanNCTestProblem;

namespace Properties {

// Create new type tags
namespace TTag {
struct MaxwellStefanNCTest { using InheritsFrom = std::tuple<NavierStokesNC, StaggeredFreeFlowModel>; };
} // end namespace TTag

template<class TypeTag>
struct ReplaceCompEqIdx<TypeTag, TTag::MaxwellStefanNCTest> { static constexpr int value = 0; };

// Set the grid type
template<class TypeTag>
struct Grid<TypeTag, TTag::MaxwellStefanNCTest> { using type = Dune::YaspGrid<2>; };

// Set the problem property
template<class TypeTag>
struct Problem<TypeTag, TTag::MaxwellStefanNCTest> { using type = Dumux::MaxwellStefanNCTestProblem<TypeTag> ; };

template<class TypeTag>
struct EnableGridGeometryCache<TypeTag, TTag::MaxwellStefanNCTest> { static constexpr bool value = true; };

template<class TypeTag>
struct EnableGridFluxVariablesCache<TypeTag, TTag::MaxwellStefanNCTest> { static constexpr bool value = true; };
template<class TypeTag>
struct EnableGridVolumeVariablesCache<TypeTag, TTag::MaxwellStefanNCTest> { static constexpr bool value = true; };

template<class TypeTag>
struct UseMoles<TypeTag, TTag::MaxwellStefanNCTest> { static constexpr bool value = true; };


//! Here we set FicksLaw or MaxwellStefansLaw
template<class TypeTag>
struct MolecularDiffusionType<TypeTag, TTag::MaxwellStefanNCTest> { using type = MaxwellStefansLaw<TypeTag>; };


/*!
 * \ingroup NavierStokesNCTests
 * \brief  A simple fluid system with three components for testing the  multi-component diffusion with the Maxwell-Stefan formulation.
 */
template<class TypeTag>
class MaxwellStefanFluidSystem
: public FluidSystems::Base<GetPropType<TypeTag, Properties::Scalar>, MaxwellStefanFluidSystem<TypeTag>>

{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ThisType = MaxwellStefanFluidSystem<TypeTag>;
    using Base = FluidSystems::Base<Scalar, ThisType>;

public:
    //! The number of phases
    static constexpr int numPhases = 1;
    static constexpr int numComponents = 3;

    static constexpr int H2Idx = 0;//first major component
    static constexpr int N2Idx = 1;//second major component
    static constexpr int CO2Idx = 2;//secondary component

    //! Human readable component name (index compIdx) (for vtk output)
    static std::string componentName(int compIdx)
    { return "MaxwellStefan_" + std::to_string(compIdx); }

    //! Human readable phase name (index phaseIdx) (for velocity vtk output)
    static std::string phaseName(int phaseIdx = 0)
    { return "Gas"; }

    //! Molar mass in kg/mol of the component with index compIdx
    static Scalar molarMass(unsigned int compIdx)
    { return 0.02896; }


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
     *
     */
    template <class FluidState>
    static Scalar density(const FluidState &fluidState,
                          const int phaseIdx)
    {
      return 1;
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
        return density(fluidState, phaseIdx)/molarMass(0);
    }
};

template<class TypeTag>
struct FluidSystem<TypeTag, TTag::MaxwellStefanNCTest> { using type = MaxwellStefanFluidSystem<TypeTag>; };

} // end namespace Properties

/*!
 * \ingroup NavierStokesNCTests
 * \brief Test problem for the Maxwell-Stefan model
 */
template <class TypeTag>
class MaxwellStefanNCTestProblem : public NavierStokesProblem<TypeTag>
{
    using ParentType = NavierStokesProblem<TypeTag>;

    using BoundaryTypes = Dumux::NavierStokesBoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using NumEqVector = GetPropType<TypeTag, Properties::NumEqVector>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using TimeLoopPtr = std::shared_ptr<CheckPointTimeLoop<Scalar>>;

    static constexpr auto compOneIdx =  Indices::conti0EqIdx;
    static constexpr auto compTwoIdx =  Indices::conti0EqIdx + FluidSystem::N2Idx;
    static constexpr auto compThreeIdx = Indices::conti0EqIdx + FluidSystem::CO2Idx;

public:
    MaxwellStefanNCTestProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        plotOutput_ = getParam<bool>("Problem.PlotOutput", false);
    }

   /*!
     * \name Problem parameters
     */
    // \{

   /*!
     * \brief Writes out the diffusion rates from left to right
     *
     * Called after every time step
     *
     * \param curSol Vector containing the current solution
     * \param gridVariables The grid variables
     * \param time The time
     */
     template<class SolutionVector, class GridVariables>
     void plotComponentsOverTime(const SolutionVector& curSol,
                                 const GridVariables& gridVariables,
                                 const Scalar time)
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
            for (const auto& element : elements(this->gridGeometry().gridView()))
            {
                auto fvGeometry = localView(this->gridGeometry());
                fvGeometry.bindElement(element);

                auto elemVolVars = localView(gridVariables.curGridVolVars());
                elemVolVars.bind(element, fvGeometry, curSol);
                for (auto&& scv : scvs(fvGeometry))
                {
                    const auto globalPos = scv.dofPosition();

                    if (globalPos[0] < 0.5)
                    {
                        x_co2_left += elemVolVars[scv].moleFraction(FluidSystem::CO2Idx);
                        x_n2_left += elemVolVars[scv].moleFraction(FluidSystem::N2Idx);
                        x_h2_left += elemVolVars[scv].moleFraction(FluidSystem::H2Idx);
                        i +=1;
                    }
                    else
                    {
                        x_co2_right += elemVolVars[scv].moleFraction(FluidSystem::CO2Idx);
                        x_n2_right += elemVolVars[scv].moleFraction(FluidSystem::N2Idx);
                        x_h2_right += elemVolVars[scv].moleFraction(FluidSystem::H2Idx);
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
            y1_.push_back(x_n2_left);
            y2_.push_back(x_n2_right);
            y3_.push_back(x_co2_left);
            y4_.push_back(x_co2_right);
            y5_.push_back(x_h2_left);
            y6_.push_back(x_h2_right);

            gnuplot_.resetPlot();
            gnuplot_.setXRange(0, std::max(time, 72000.0));
            gnuplot_.setYRange(0.4, 0.6);
            gnuplot_.setXlabel("time [s]");
            gnuplot_.setYlabel("mole fraction mol/mol");
            gnuplot_.addDataSetToPlot(x_, y1_, "N2_left.dat", "w l t 'N_2 left'");
            gnuplot_.addDataSetToPlot(x_, y2_, "N2_right.dat", "w l t 'N_2 right'");
            gnuplot_.plot("mole_fraction_N2");

            gnuplot2_.resetPlot();
            gnuplot2_.setXRange(0, std::max(time, 72000.0));
            gnuplot2_.setYRange(0.0, 0.6);
            gnuplot2_.setXlabel("time [s]");
            gnuplot2_.setYlabel("mole fraction mol/mol");
            gnuplot2_.addDataSetToPlot(x_, y3_, "CO2_left.dat", "w l t 'CO_2 left'");
            gnuplot2_.addDataSetToPlot(x_, y4_, "CO2_right.dat", "w l t 'CO_2 right'");
            gnuplot2_.plot("mole_fraction_C02");

            gnuplot3_.resetPlot();
            gnuplot3_.setXRange(0, std::max(time, 72000.0));
            gnuplot3_.setYRange(0.0, 0.6);
            gnuplot3_.setXlabel("time [s]");
            gnuplot3_.setYlabel("mole fraction mol/mol");
            gnuplot3_.addDataSetToPlot(x_, y5_, "H2_left.dat", "w l t 'H_2 left'");
            gnuplot3_.addDataSetToPlot(x_, y6_, "H2_right.dat", "w l t 'H_2 right'");
            gnuplot3_.plot("mole_fraction_H2");
        }
    }

   /*!
     * \brief Returns the temperature within the domain in [K].
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperature() const
    { return 273.15 + 10; } // 10C

    // \}
   /*!
     * \name Boundary conditions
     */
    // \{

   /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary control volume.
     *
     * \param globalPos The position of the center of the finite volume
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;
        // Set no-flow/no-slip everywhere. Neumann defaults to zero.
        values.setDirichlet(Indices::velocityXIdx);
        values.setDirichlet(Indices::velocityYIdx);
        values.setNeumann(compOneIdx);
        values.setNeumann(compTwoIdx);
        values.setNeumann(compThreeIdx);
        return values;
    }

   /*!
     * \brief Evaluates the boundary conditions for a Dirichlet control volume.
     *
     * \param globalPos The center of the finite volume which ought to be set.
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        return initialAtPos(globalPos);
    }

   /*!
     * \name Volume terms
     */
    // \{

   /*!
     * \brief Evaluates the initial value for a control volume.
     *
     * \param globalPos The global position
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables initialValues(0.0);
        if (globalPos[0] < 0.5)
        {
           initialValues[compTwoIdx] = 0.50086;
           initialValues[compThreeIdx] = 0.49914;
        }
        else
        {
           initialValues[compTwoIdx] = 0.49879;
           initialValues[compThreeIdx] = 0.0;
        }

        initialValues[Indices::pressureIdx] = 1.1e+5;
        initialValues[Indices::velocityXIdx] = 0.0;
        initialValues[Indices::velocityYIdx] = 0.0;

        return initialValues;
    }

private:

    bool plotOutput_;

    const Scalar eps_{1e-6};

    Dumux::GnuplotInterface<Scalar> gnuplot_;
    Dumux::GnuplotInterface<Scalar> gnuplot2_;
    Dumux::GnuplotInterface<Scalar> gnuplot3_;

    std::vector<Scalar> x_;
    std::vector<Scalar> y1_;
    std::vector<Scalar> y2_;
    std::vector<Scalar> y3_;
    std::vector<Scalar> y4_;
    std::vector<Scalar> y5_;
    std::vector<Scalar> y6_;
};
} // end namespace Dumux

#endif
