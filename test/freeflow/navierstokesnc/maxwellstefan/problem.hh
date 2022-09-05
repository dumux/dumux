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

#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/numeqvector.hh>

#include <dumux/freeflow/navierstokes/boundarytypes.hh>
#include <dumux/freeflow/navierstokes/staggered/problem.hh>

#include <dumux/io/gnuplotinterface.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesNCTests
 * \brief Test problem for the Maxwell-Stefan model
 */
template <class TypeTag>
class MaxwellStefanNCTestProblem : public NavierStokesStaggeredProblem<TypeTag>
{
    using ParentType = NavierStokesStaggeredProblem<TypeTag>;

    using BoundaryTypes = Dumux::NavierStokesBoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

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
            auto fvGeometry = localView(this->gridGeometry());
            auto elemVolVars = localView(gridVariables.curGridVolVars());
            for (const auto& element : elements(this->gridGeometry().gridView()))
            {
                fvGeometry.bindElement(element);
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
