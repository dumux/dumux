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
 *   along with this program. If not, see <http://www.gnu.org/licenses/>.    *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief A one-dimensional infiltration problem with a smoth, given solution.
 *
 * The source term is calculated analytically. Thus, this example can be used
 * to calclate the L2 error and to show convergence for grid and time-step
 * refinement.
 */
#ifndef DUMUX_RICHARDS_ANALYTICALPROBLEM_HH
#define DUMUX_RICHARDS_ANALYTICALPROBLEM_HH

#include <cmath>
#include <dune/geometry/type.hh>
#include <dune/geometry/quadraturerules.hh>
#include <dune/grid/io/file/dgfparser.hh>

#include <dumux/porousmediumflow/richards/implicit/model.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/liquidphase.hh>

#include "richardsanalyticalspatialparams.hh"

namespace Dumux
{

template <class TypeTag>
class RichardsAnalyticalProblem;

//////////
// Specify the properties for the analytical problem
//////////
namespace Properties
{
NEW_TYPE_TAG(RichardsAnalyticalProblem, INHERITS_FROM(Richards, RichardsAnalyticalSpatialParams));
NEW_TYPE_TAG(RichardsAnalyticalBoxProblem, INHERITS_FROM(BoxModel, RichardsAnalyticalProblem));
NEW_TYPE_TAG(RichardsAnalyticalCCProblem, INHERITS_FROM(CCModel, RichardsAnalyticalProblem));

// Use 2d YaspGrid
SET_TYPE_PROP(RichardsAnalyticalProblem, Grid, Dune::YaspGrid<2>);

// Set the physical problem to be solved
SET_PROP(RichardsAnalyticalProblem, Problem)
{ typedef Dumux::RichardsAnalyticalProblem<TypeTag> type; };

// Set the wetting phase
SET_PROP(RichardsAnalyticalProblem, WettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Dumux::FluidSystems::LiquidPhase<Scalar, Dumux::SimpleH2O<Scalar> > type;
};

// Enable gravity
SET_BOOL_PROP(RichardsAnalyticalProblem, ProblemEnableGravity, true);

// Enable partial reassembly of the Jacobian matrix
SET_BOOL_PROP(RichardsAnalyticalProblem, ImplicitEnablePartialReassemble, true);

// Enable re-use of the Jacobian matrix for the first iteration of a time step
SET_BOOL_PROP(RichardsAnalyticalProblem, ImplicitEnableJacobianRecycling, true);

// Use forward differences to approximate the Jacobian matrix
SET_INT_PROP(RichardsAnalyticalProblem, ImplicitNumericDifferenceMethod, +1);

// Set the maximum number of newton iterations of a time step
SET_INT_PROP(RichardsAnalyticalProblem, NewtonMaxSteps, 28);

// Set the "desireable" number of newton iterations of a time step
SET_INT_PROP(RichardsAnalyticalProblem, NewtonTargetSteps, 18);

// Do not write the intermediate results of the newton method
SET_BOOL_PROP(RichardsAnalyticalProblem, NewtonWriteConvergence, false);
}

/*!
 * \ingroup RichardsModel
 * \ingroup ImplicitTestProblems
 *
 *\brief A water infiltration problem using Richards model and comparing
 *        to an analytical solution. Implemented by using the source term
 *        defined as the analytical solution.
 *
 * The domain is box shaped. Top and bottom boundaries are Dirichlet
 * boundaries with fixed water pressure (fixed Saturation \f$S_w = 0\f$),
 * left and right boundary are closed (Neumann 0 boundary).
 * This problem uses the \ref RichardsModel
 *
 * The L2 error is decreasing with decreasing time and space discretization.
 */
template <class TypeTag>
class RichardsAnalyticalProblem : public RichardsProblem<TypeTag>
{
    typedef RichardsProblem<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GridView::template Codim<0>::Entity::Geometry Geometry;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        // copy some indices for convenience
        pwIdx = Indices::pwIdx,
        contiEqIdx = Indices::contiEqIdx
    };
    enum {
        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:
    /*!
     * \brief Constructor
     *
     * \param timeManager The Dumux TimeManager for simulation management.
     * \param gridView The grid view on the spatial domain of the problem
     */
    RichardsAnalyticalProblem(TimeManager &timeManager,
                        const GridView &gridView)
        : ParentType(timeManager, gridView)
    {
        eps_ = 3e-6;
        pnRef_ = 1e5; // air pressure
        name_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Problem, Name);
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief The problem name
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const std::string name() const
    { return name_; }

    /*!
     * \brief Returns the temperature [K] within a finite volume
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperature() const
    { return 273.15 + 10; } // -> 10°C

    /*!
     * \brief Returns the reference pressure [Pa] of the non-wetting
     *        fluid phase within a finite volume
     *
     * This problem assumes a constant reference pressure of 1 bar.
     *
     * \param element The DUNE Codim<0> entity which intersects with
     *                the finite volume in question
     * \param fvGeometry The finite volume geometry of the element
     * \param scvIdx The sub control volume index inside the finite
     *               volume geometry
     */
    Scalar referencePressure(const Element &element,
                             const FVElementGeometry &fvGeometry,
                             const int scvIdx) const
    { return pnRef_; }

   /*!
     * \brief Evaluate the source values for a control volume.
     *
     * For this method, the \a values parameter stores primary
     * variables. For this test case, the analytical solution is
     * used to calculate the source term. See the Matlab script
     * Richards.m which uses Matlab's Symbolic Toolbox to calclate
     * the source term.
     *
     * \param values Storage for all primary variables of the source term
     * \param globalPos The position for which the source term is set
     */
    void sourceAtPos(PrimaryVariables &values,
                const GlobalPosition &globalPos) const
    {
        const Scalar time = this->timeManager().time() + this->timeManager().timeStepSize();
        const Scalar pwTop = 98942.8;
        const Scalar pwBottom = 95641.1;

        // linear model with complex solution
        // calcluated with Matlab script "Richards.m"
        values = (std::pow(std::tanh(globalPos[1]*5.0+time*(1.0/1.0E1)-1.5E1),2.0)*(1.0/1.0E1)
            -1.0/1.0E1)*(pwBottom*(1.0/2.0)-pwTop*(1.0/2.0))*4.0E-8-((std::pow(std::tanh(globalPos[1]
            *5.0+time*(1.0/1.0E1)-1.5E1),2.0)*5.0-5.0)*(pwBottom*(1.0/2.0)-pwTop*(1.0/2.0))-1.0E3)
            *(std::pow(std::tanh(globalPos[1]*5.0+time*(1.0/1.0E1)-1.5E1),2.0)*5.0-5.0)*(pwBottom
            *(1.0/2.0)-pwTop*(1.0/2.0))*5.0E-16+std::tanh(globalPos[1]*5.0+time*(1.0/1.0E1)-1.5E1)
            *(std::pow(std::tanh(globalPos[1]*5.0+time*(1.0/1.0E1)-1.5E1),2.0)*5.0-5.0)*(pwBottom
            *(1.0/2.0)-pwTop*(1.0/2.0))*(pwBottom*5.0E-16-(std::tanh(globalPos[1]*5.0+time*(1.0/1.0E1)
            -1.5E1)+1.0)*(pwBottom*(1.0/2.0)-pwTop*(1.0/2.0))*5.0E-16+4.99995E-6)*1.0E1;
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
     * \param values The boundary types for the conservation equations
     * \param globalPos The position for which the boundary type is set
     */
    void boundaryTypesAtPos(BoundaryTypes &values,
                       const GlobalPosition &globalPos) const
    {
        if (onLowerBoundary_(globalPos) ||
            onUpperBoundary_(globalPos))
        {
            values.setAllDirichlet();
        }
        else
            values.setAllNeumann();
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        boundary segment.
     *
     * \param values The dirichlet values for the primary variables
     * \param globalPos The position for which the Dirichlet value is set
     *
     * For this method, the \a values parameter stores primary variables.
     */
    void dirichletAtPos(PrimaryVariables &values,
                   const GlobalPosition &globalPos) const
    {
        const Scalar time = this->timeManager().time() + this->timeManager().timeStepSize();
        const Scalar pwTop = 98942.8;
        const Scalar pwBottom = 95641.1;
        Scalar pw = pwBottom
          + 0.5 * (std::tanh( (5.0 * globalPos[1]) - 15.0 + time/10.0) + 1.0) * (pwTop - pwBottom);

        values[pwIdx] = pw;
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * For this method, the \a values parameter stores the mass flux
     * in normal direction of each phase. Negative values mean influx.
     *
     * \param values The neumann values for the conservation equations
     * \param globalPos The position for which the Neumann value is set
     */
    void neumannAtPos(PrimaryVariables &values,
                 const GlobalPosition &globalPos) const
    {
        values = 0.0;
    }

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluate the initial values for a control volume.
     *
     * For this method, the \a values parameter stores primary
     * variables.
     *
     * \param values Storage for all primary variables of the initial condition
     * \param globalPos The position for which the boundary type is set
     */
    void initialAtPos(PrimaryVariables &values,
                 const GlobalPosition &globalPos) const
    {
        const Scalar time = this->timeManager().time();
        analyticalSolution(values, time, globalPos);
    }

    // \}

    /*!
     * \brief Evaluate the analytical solution.
     *
     * \param values The dirichlet values for the primary variables
     * \param time The time at wich the solution should be evaluated
     * \param globalPos The position for which the Dirichlet value is set
     *
     * For this method, the \a values parameter stores primary variables.
     */
    void analyticalSolution(PrimaryVariables &values,
                   const Scalar time,
                   const GlobalPosition &globalPos) const
    {
        const Scalar pwTop = 98942.8;
        const Scalar pwBottom = 95641.1;
        Scalar pw = pwBottom
          + 0.5 * (std::tanh( (5.0 * globalPos[1]) - 15.0 + time/10.0) + 1.0) * (pwTop - pwBottom);

        values[pwIdx] = pw;
    }

    /*!
     * \brief Calculate the L2 error between the solution given by
     *        dirichletAtPos and the numerical approximation.
     * \note Works for cell-centered FV only because the numerical
     *       approximation is only evaluated in the cell center (once).
     *       To extend this function to the box method the evaluation
     *       has to be exted to box' subvolumes.
     */
    Scalar calculateL2Error()
    {
        const unsigned int qOrder = 4;
        Scalar l2error = 0.0;
        Scalar l2analytic = 0.0;
        const Scalar time = this->timeManager().time() + this->timeManager().timeStepSize();

        for (const auto& element : elements(this->gridView()))
        {
            // value from numerical approximation
            Scalar numericalSolution =
                this->model().curSol()[this->model().dofMapper().subIndex(element, 0, 0)];

            // integrate over element using a quadrature rule
            Geometry geometry = element.geometry();
            Dune::GeometryType gt = geometry.type();
            Dune::QuadratureRule<Scalar, dim> rule =
                Dune::QuadratureRules<Scalar, dim>::rule(gt, qOrder);

            for (auto qIt = rule.begin(); qIt != rule.end(); ++qIt)
            {
                // evaluate analytical solution
                Dune::FieldVector<Scalar, dim> globalPos = geometry.global(qIt->position());
                Dune::FieldVector<Scalar, 1> values(0.0);
                analyticalSolution(values, time, globalPos);
                // add contributino of current quadrature point
                l2error += (numericalSolution - values[0]) * (numericalSolution - values[0]) *
                    qIt->weight() * geometry.integrationElement(qIt->position());
                l2analytic += values[0] * values[0] *
                    qIt->weight() * geometry.integrationElement(qIt->position());
            }
        }
        return std::sqrt(l2error/l2analytic);
    }

    /*!
     * \brief Write the relevant secondary variables of the current
     *        solution into an VTK output file.
     */
    void writeOutput(const bool verbose = true)
    {
        ParentType::writeOutput(verbose);

        Scalar l2error = calculateL2Error();

        // compute L2 error if analytical solution is available
        std::cout.precision(8);
        std::cout << "L2 error for "
                  << std::setw(6) << this->gridView().size(0)
                  << " elements: "
                  << std::scientific
                  << l2error
                  << std::endl;
    }

    /*!
     * \brief If we should write output
     */
    bool shouldWriteOutput()
    {
        return this->timeManager().willBeFinished();
    }

    /*!
     * \brief If we should write output
     */
    bool shouldWriteRestartFile()
    {
        return false;
    }

private:

    // evalutates if global position is at lower boundary
    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[1] < this->bBoxMin()[1] + eps_;
    }

    // evalutates if global position is at upper boundary
    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[1] > this->bBoxMax()[1] - eps_;
    }

    Scalar eps_;
    Scalar pnRef_;
    std::string name_;
};
} //end namespace

#endif
