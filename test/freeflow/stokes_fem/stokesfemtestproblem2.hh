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
 * @file
 * @brief  Definition of a simple Stokes problem
 */
#ifndef DUMUX_STOKESFEMTESTPROBLEM_HH
#define DUMUX_STOKESFEMTESTPROBLEM_HH

// fluid -- TODO remove superfluous files
#include <dumux/material/fluidsystems/h2on2.hh>
#include <dumux/material/fluidsystems/gasphase.hh>
#include <dumux/material/fluidsystems/liquidphase.hh>
#include <dumux/material/components/constant.hh>
#include <dumux/material/components/simpleh2o.hh>

// fem related
#include <dumux/implicit/fem/properties.hh>
#include <dumux/implicit/fem/problem.hh>
// adapted from Stokes (box)
#include <dumux/freeflow/stokes_fem/model.hh>

#include <dumux/linear/seqsolverbackend.hh>
#include "elasticspatialparams.hh" // TODO check which methods are necessary for Stokes

namespace Dumux
{

template <class TypeTag>
class StokesFemTestProblem2;

//////////
// Specify the properties for the stokes problem
//////////
namespace Properties
{//valid to inherit from both?
NEW_TYPE_TAG(StokesFemTestProblem2, INHERITS_FROM(BoxStokes, FemModel, ElSpatialParams));

// Set the grid type
SET_TYPE_PROP(StokesFemTestProblem2, Grid, Dune::YaspGrid<2>);

// Set the problem property
SET_TYPE_PROP(StokesFemTestProblem2, Problem, StokesFemTestProblem2<TypeTag>);

////// Use nitrogen as gas phase
//SET_TYPE_PROP(StokesFemTestProblem, Fluid,
//              FluidSystems::GasPhase<typename GET_PROP_TYPE(TypeTag, Scalar),
//                                            N2<typename GET_PROP_TYPE(TypeTag, Scalar)> >);

////// Use nitrogen as gas phase
//SET_TYPE_PROP(StokesFemTestProblem2, Fluid,
//              FluidSystems::LiquidPhase<typename GET_PROP_TYPE(TypeTag, Scalar),
//                                        SimpleH2O<typename GET_PROP_TYPE(TypeTag, Scalar)> >);

// Use nitrogen as gas phase
//SET_TYPE_PROP(StokesFemTestProblem2, Fluid,
//              FluidSystems::GasPhase<typename GET_PROP_TYPE(TypeTag, Scalar),
//                                            N2<typename GET_PROP_TYPE(TypeTag, Scalar)> >);

// Use nitrogen as gas phase
SET_TYPE_PROP(StokesFemTestProblem2, Fluid,
              FluidSystems::LiquidPhase<typename GET_PROP_TYPE(TypeTag, Scalar),
                                            Constant<TypeTag ,typename GET_PROP_TYPE(TypeTag, Scalar)> >);


//added from elastic
// Quadrature order
SET_INT_PROP(StokesFemTestProblem2, FemBasisOrder, 1);

SET_TYPE_PROP(StokesFemTestProblem2, LinearSolver, ILUnBiCGSTABBackend<TypeTag>);
}

/*!
 * \ingroup BoxStokesModel
 * \ingroup ImplicitTestProblems
 * \brief Stokes flow problem
 *
 * This problem uses the \ref BoxStokesModel.
 * To run the simulation execute the following line in shell:
 * <tt>./test_stokes</tt>
 */
template <class TypeTag>
class StokesFemTestProblem2 : public ImplicitFemProblem<TypeTag>
{
    using ParentType = ImplicitFemProblem<TypeTag>;

    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using TimeManager = typename GET_PROP_TYPE(TypeTag, TimeManager);
    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);
    using IpData = typename GET_PROP_TYPE(TypeTag, FemIntegrationPointData);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    //needed?
    using SecondaryVariables = typename GET_PROP_TYPE(TypeTag, SecondaryVariables);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);


    enum {
        // Number of equations and grid dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,

        // copy some indices for convenience
        massBalanceIdx = Indices::massBalanceIdx, //!< Index of the mass balance
        momentumXIdx = Indices::momentumXIdx, //!< Index of the x-component of the momentum balance
        momentumYIdx = Indices::momentumYIdx //!< Index of the y-component of the momentum balance
    };
    enum { // indices of the primary variables
        velocityXIdx = Indices::velocityXIdx, //!< Index of the x-velocity
        velocityYIdx = Indices::velocityYIdx, //!< Index of the y-velocity
        pressureIdx = Indices::pressureIdx //!< Index of the pressure
    };

    using Element = typename GridView::template Codim<0>::Entity;
    using CoordScalar = typename GridView::ctype;
    using Intersection = typename GridView::Intersection;

    //typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    using Fluid = typename GET_PROP_TYPE(TypeTag, Fluid);
    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;


public:
    StokesFemTestProblem2(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView)
    {
        eps_ = 1e-6;
        inletVelocity_ = GET_RUNTIME_PARAM(TypeTag, Scalar, Problem.InletVelocity);
        kinematicViscosity_ = GET_RUNTIME_PARAM(TypeTag, Scalar, Component.LiquidKinematicViscosity);
    }

    /*!
     * \brief The problem name.
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const char *name() const
    { return "stokes_fem2"; }

    /*!
     * \brief Returns the temperature within the domain.
     *
     * This problem assumes a constant temperature of 10 degrees Celsius.
     */
    Scalar temperatureAtPos(const GlobalPosition &globalPos) const
    {
        return 273.15 + 10; // -> 10C
    }

    Scalar temperature() const
        {  return 273.15 + 10; }


//______________________________________________________________________________________________________________________________________________________
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        values.setAllDirichlet();

        // set pressure at one point
        const Scalar middle = (this->bBoxMax()[1] - this->bBoxMin()[1])/2;
        // if (onRightBoundary_(globalPos) && globalPos[1] > middle - eps_ && globalPos[1] < middle + eps_)
        //     values.setDirichlet(pressureIdx);

        // if (onRightBoundary_(globalPos))
        //     values.setNeumann(Indices::velocityXIdx);

        return values;
    }


   PrimaryVariables neumannAtPos(const GlobalPosition &globalPos) const
   {
       PrimaryVariables values(0.0);
       const Scalar time = this->timeManager().time() + this->timeManager().timeStepSize();
       const Scalar velocityVariation = 0.2;

       if (onRightBoundary_(globalPos))
       values[Indices::velocityXIdx] = 0;

       return values;
   }

    // old STOKES TEST PROBLEM
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values(0.0);
        const Scalar time = this->timeManager().time() + this->timeManager().timeStepSize();
        const Scalar velocityVariation = 0.2;

        values = initialAtPos(globalPos);

        // sinusoidal variation of the maximum velocity in time
        const Scalar v0 = inletVelocity_;

        // parabolic velocity profile
        values[velocityXIdx] =  v0*(globalPos[1] - this->bBoxMin()[1])*(this->bBoxMax()[1] - globalPos[1])
                               / (0.25*(this->bBoxMax()[1] - this->bBoxMin()[1])*(this->bBoxMax()[1] - this->bBoxMin()[1]));
        return values;
    }


    /*!
     * \brief Return the analytical solution of the problem at a given position
     *
     * \param globalPos The global position
     */
    PrimaryVariables analyticalSolution(const GlobalPosition& globalPos) const
    {
        Scalar x = globalPos[0];
        Scalar y = globalPos[1];

        const Scalar v0 = inletVelocity_;


        PrimaryVariables values(0.0);


        values[pressureIdx] = x * (1.0-x); // p(x,y) = x(1-x) [Donea2003]
        values[velocityXIdx] = v0*(globalPos[1] - this->bBoxMin()[1])*(this->bBoxMax()[1] - globalPos[1])
                                  / (0.25*(this->bBoxMax()[1] - this->bBoxMin()[1])*(this->bBoxMax()[1] - this->bBoxMin()[1]));
        values[velocityYIdx] = 0;

        return values;
    }


    // old STOKES TEST PROBLEM
    PrimaryVariables sourceAtPos(const GlobalPosition &globalPos) const
    {
        return PrimaryVariables(0.0);
    }


    // old STOKES TEST PROBLEM
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values(0.0);
        values[pressureIdx] = 1e5;
        values[velocityXIdx] =  0;
        values[velocityYIdx] =  0;
        return values;
    }


//______________________________________________________________________________________________________________________________________________________


private:
    bool onLeftBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] < this->bBoxMin()[0] + eps_; }

    bool onRightBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] > this->bBoxMax()[0] - eps_; }

    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] < this->bBoxMin()[1] + eps_; }

    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] > this->bBoxMax()[1] - eps_; }

    Scalar eps_, inletVelocity_, kinematicViscosity_;

};

} //end namespace

#endif
