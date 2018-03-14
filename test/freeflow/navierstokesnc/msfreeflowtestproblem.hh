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
 * \ingroup NavierStokesNCTests
 * \brief Channel flow test for the multi-component staggered grid (Navier-)Stokes model
 */
#ifndef DUMUX_CHANNEL_MAXWELL_STEFAN_TEST_PROBLEM_HH
#define DUMUX_CHANNEL_MAXWELL_STEFAN_TEST_PROBLEM_HH

#include <dune/grid/yaspgrid.hh>

#include <dumux/freeflow/compositional/navierstokesncmodel.hh>
#include <dumux/freeflow/navierstokes/problem.hh>

#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/h2oair.hh>

#include <dumux/discretization/maxwellstefanslaw.hh>
#include <dumux/discretization/staggered/freeflow/properties.hh>

#include <dumux/io/gnuplotinterface.hh>

namespace Dumux
{
template <class TypeTag>
class MaxwellStefanNCTestProblem;

namespace Properties {

NEW_TYPE_TAG(MaxwellStefanNCTestTypeTag, INHERITS_FROM(StaggeredFreeFlowModel, NavierStokesNC));

NEW_PROP_TAG(FluidSystem);

SET_INT_PROP(MaxwellStefanNCTestTypeTag, ReplaceCompEqIdx, 0);

// Set the grid type
SET_TYPE_PROP(MaxwellStefanNCTestTypeTag, Grid, Dune::YaspGrid<2>);

// Set the problem property
SET_TYPE_PROP(MaxwellStefanNCTestTypeTag, Problem, Dumux::MaxwellStefanNCTestProblem<TypeTag> );

SET_BOOL_PROP(MaxwellStefanNCTestTypeTag, EnableFVGridGeometryCache, true);

SET_BOOL_PROP(MaxwellStefanNCTestTypeTag, EnableGridFluxVariablesCache, true);
SET_BOOL_PROP(MaxwellStefanNCTestTypeTag, EnableGridVolumeVariablesCache, true);

SET_BOOL_PROP(MaxwellStefanNCTestTypeTag, UseMoles, true);

// #if ENABLE_NAVIERSTOKES
SET_BOOL_PROP(MaxwellStefanNCTestTypeTag, EnableInertiaTerms, false);

//! Here we set FicksLaw or MaxwellStefansLaw
SET_TYPE_PROP(MaxwellStefanNCTestTypeTag, MolecularDiffusionType, MaxwellStefansLaw<TypeTag>);


/*!
 * \ingroup NavierStokesNCTests
 * \brief  A simple fluid system with one MaxwellStefan component.
 * \todo doc me!
 */
template<class TypeTag>
class MaxwellStefanFluidSystem: public FluidSystems::BaseFluidSystem<typename GET_PROP_TYPE(TypeTag, Scalar),MaxwellStefanFluidSystem<TypeTag>>

{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using ThisType = MaxwellStefanFluidSystem<TypeTag>;
    using Base = FluidSystems::BaseFluidSystem<Scalar, ThisType>;

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
     *        return the binary diffusion coefficient \f$\mathrm{[m^2/s]}\f$ for components
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
                                             unsigned int compIIdx,
                                             unsigned int compJIdx)
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
     *        the partial pressures of all components, return its
     *        density \f$\mathrm{[kg/m^3]}\f$.
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
     * \brief Calculate the dynamic viscosity of a fluid phase \f$\mathrm{[Pa*s]}\f$
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

SET_TYPE_PROP(MaxwellStefanNCTestTypeTag, FluidSystem, MaxwellStefanFluidSystem<TypeTag>);

} //end namespace Property
/*!
 * \brief  Test problem for the maxwell stefan model
 * \todo doc me!
 */
template <class TypeTag>
class MaxwellStefanNCTestProblem : public NavierStokesProblem<TypeTag>
{
    using ParentType = NavierStokesProblem<TypeTag>;

    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);
    using Indices = typename GET_PROP_TYPE(TypeTag, ModelTraits)::Indices;
    using NumEqVector = typename GET_PROP_TYPE(TypeTag, NumEqVector);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);

    using Element = typename FVGridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using TimeLoopPtr = std::shared_ptr<CheckPointTimeLoop<Scalar>>;

    enum {
        compTwoIdx =  Indices::conti0EqIdx + FluidSystem::N2Idx,
        compThreeIdx = Indices::conti0EqIdx + FluidSystem::CO2Idx,
    };

public:
    MaxwellStefanNCTestProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
    {
        plotOutput_ = getParam<bool>("Problem.PlotOutput", false);
    }

   /*!
     * \name Problem parameters
     */
    // \{

   /*!
     * \brief Output the diffusion rates from left to right
     *
     * Called after every time step
     *
     * \param curSol Vector containing the current solution
     * \param gridVariables The grid variables
     * \param time The time
     */
     template<class SolutionVector, class GridVariables>
     void postTimeStep(const SolutionVector& curSol,
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
            for (const auto& element : elements(this->fvGridGeometry().gridView()))
            {
                auto fvGeometry = localView(this->fvGridGeometry());
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
     * \brief Return the temperature within the domain in [K].
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperature() const
    { return 273.15 + 10; } // 10C

   /*!
     * \brief Return the sources within the domain.
     *
     * \param globalPos The global position
     */
    NumEqVector sourceAtPos(const GlobalPosition &globalPos) const
    {
        return NumEqVector(0.0);
    }

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
        // set Dirichlet values for the velocity everywhere
        values.setDirichlet(Indices::velocityXIdx);
        values.setDirichlet(Indices::velocityYIdx);
        values.setOutflow(compTwoIdx);
        values.setOutflow(compThreeIdx);
        return values;
    }

   /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume.
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
     * \brief Evaluate the initial value for a control volume.
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
} //end namespace

#endif
