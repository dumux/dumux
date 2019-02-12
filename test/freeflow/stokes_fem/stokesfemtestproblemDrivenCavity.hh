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
class StokesFemTestProblemDrivenCavity;

//////////
// Specify the properties for the stokes problem
//////////
namespace Properties
{//valid to inherit from both?
NEW_TYPE_TAG(StokesFemTestProblemDrivenCavity, INHERITS_FROM(BoxStokes, FemModel, ElSpatialParams));

// Set the grid type
SET_TYPE_PROP(StokesFemTestProblemDrivenCavity, Grid, Dune::YaspGrid<2>);

// Set the problem property
SET_TYPE_PROP(StokesFemTestProblemDrivenCavity, Problem, StokesFemTestProblemDrivenCavity<TypeTag>);

////// Use nitrogen as gas phase
//SET_TYPE_PROP(StokesFemTestProblem, Fluid,
//              FluidSystems::GasPhase<typename GET_PROP_TYPE(TypeTag, Scalar),
//                                            N2<typename GET_PROP_TYPE(TypeTag, Scalar)> >);

// Use nitrogen as gas phase
//SET_TYPE_PROP(StokesFemTestProblem2, Fluid,
//              FluidSystems::GasPhase<typename GET_PROP_TYPE(TypeTag, Scalar),
//                                            N2<typename GET_PROP_TYPE(TypeTag, Scalar)> >);

// Use nitrogen as gas phase
SET_TYPE_PROP(StokesFemTestProblemDrivenCavity, Fluid,
              FluidSystems::LiquidPhase<typename GET_PROP_TYPE(TypeTag, Scalar),
                                            Constant<TypeTag ,typename GET_PROP_TYPE(TypeTag, Scalar)> >);


////// Use simple h2o
//SET_TYPE_PROP(StokesFemTestProblemDrivenCavity, Fluid,
//              FluidSystems::LiquidPhase<typename GET_PROP_TYPE(TypeTag, Scalar),
//                                        SimpleH2O<typename GET_PROP_TYPE(TypeTag, Scalar)> >);

//added from elastic
// Quadrature order
SET_INT_PROP(StokesFemTestProblemDrivenCavity, FemBasisOrder, 1);

SET_TYPE_PROP(StokesFemTestProblemDrivenCavity, LinearSolver, ILUnBiCGSTABBackend<TypeTag>);
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
class StokesFemTestProblemDrivenCavity : public ImplicitFemProblem<TypeTag>
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

 //   typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    using DimVector = Dune::FieldVector<Scalar, dim>;

public:
    StokesFemTestProblemDrivenCavity(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView)
    {
        eps_ = 1e-6;
        lidVelocity_ = GET_RUNTIME_PARAM(TypeTag, Scalar, Problem.LidVelocity);
        cellSizeX_ = GET_RUNTIME_PARAM(TypeTag, DimVector, Grid.UpperRight)[0]/GET_RUNTIME_PARAM(TypeTag, DimVector, Grid.Cells)[0];
        inletVelocity_ = 10.0;
        kinematicViscosity_ = GET_RUNTIME_PARAM(TypeTag, Scalar, Component.LiquidKinematicViscosity);
    }

    /*!
     * \brief The problem name.
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const char *name() const
    { return "stokes_femDrivenCavity"; }

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

/////////////////////////////////////////////////////////////////////////////
    /*!
        * \brief Return the sources within the domain.
        *
        * \param values Stores the source values, acts as return value
        * \param globalPos The global position
        */
       PrimaryVariables sourceAtPos(const GlobalPosition &globalPos) const
       {
           return PrimaryVariables(0.0);
       }


      /*!
        * \brief Specifies which kind of boundary condition should be
        *        used for which equation on a given boundary control volume.
        *
        * \param values The boundary types for the conservation equations
        * \param globalPos The position of the center of the finite volume
        */
       BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
       {
           BoundaryTypes values;

           // set Dirichlet values for the velocity everywhere
           values.setDirichlet(Indices::velocityXIdx);
           values.setDirichlet(Indices::velocityYIdx);

           // set a fixed pressure in one cell
           if (isLowerLeftCell_(globalPos))
               values.setDirichlet(Indices::pressureIdx);
//               values.setDirichletCell(Indices::pressureIdx);

           return values;
       }

      /*!
        * \brief Return dirichlet boundary values at a given position
        *
        * \param globalPos The global position
        */
       PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
       {
           PrimaryVariables values;
           values[Indices::pressureIdx] = 1.1e+5;
           values[Indices::velocityXIdx] = 0.0;
           values[Indices::velocityYIdx] = 0.0;

//           if(globalPos[1] > this->bBoxMax()[1] - eps_)
           if(onUpperBoundary_(globalPos))
               values[Indices::velocityXIdx] = lidVelocity_;

  //         values = analyticalSolution(globalPos);

           return values;
       }

      /*!
        * \brief Evaluate the initial value for a control volume.
        *
        * \param globalPos The global position
        */
       PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
       {
           PrimaryVariables values;
           values[Indices::pressureIdx] = 1.0e+5;
           values[Indices::velocityXIdx] = 0.0;
           values[Indices::velocityYIdx] = 0.0;

           return values;
       }



        /*!
         * \brief Return the analytical solution of the problem at a given position
         *
         * \param globalPos The global position
         */
        PrimaryVariables analyticalSolution(const GlobalPosition& globalPos) const
        {
            //Scalar kinematicViscosity_ = 0.000001;
            Scalar reynoldsNumber = 1.0 / kinematicViscosity_;

            Scalar v0 = 1.0;
            reynoldsNumber = v0 / kinematicViscosity_;

            Scalar lambda_ = 0.5 * reynoldsNumber
                             - std::sqrt(reynoldsNumber * reynoldsNumber * 0.25 + 4.0 * M_PI * M_PI);


            Scalar x = globalPos[0];
            Scalar y = globalPos[1];

            PrimaryVariables values;

            //analytical solution for Stokes problem
            //in x-direction
//            values[Indices::pressureIdx] = 0.5 * (1.0 - std::exp(2.0 * lambda_ * x));
//            values[Indices::velocityXIdx] = 1.0 - std::exp(lambda_ * x) * std::cos(2.0 * M_PI * y);
//            values[Indices::velocityYIdx] = 0.5 * lambda_ / M_PI * std::exp(lambda_ * x) * std::sin(2.0 * M_PI * y);

            //in y-direction
//            values[Indices::pressureIdx] = 0.5 * (1.0 - std::exp(2.0 * lambda_ * y));
//            values[Indices::velocityXIdx] = 0.5 * lambda_ / M_PI * std::exp(lambda_ * y) * std::sin(2.0 * M_PI * x);
//            values[Indices::velocityYIdx] = 1.0 - std::exp(lambda_ * y) * std::cos(2.0 * M_PI * x);

//            values[Indices::pressureIdx] = 0.5 * (1.0 - std::exp(2.0 * lambda_ * y));
//            values[Indices::velocityXIdx] = 1.0-y;
//            values[Indices::velocityYIdx] = 1.0-x;



            values[Indices::pressureIdx] = 0.5 * (1.0 - std::exp(2.0 * lambda_ * x));
            values[Indices::velocityXIdx] = 8*(x*x*x*x -2*x*x*x+x*x)*(4*y*y*y-2*y);
            values[Indices::velocityYIdx] = -8*(4*x*x*x-6*x*x+2*x)*(y*y*y*y-y*y);
            values *= lidVelocity_;

            return values;
        }



private:
    bool onLeftBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] < this->bBoxMin()[0] + eps_; }

    bool onRightBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] > this->bBoxMax()[0] - eps_; }

    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] < this->bBoxMin()[1] + eps_; }

    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] > this->bBoxMax()[1] - eps_; }

    bool isLowerLeftCell_(const GlobalPosition& globalPos) const
    {
        return globalPos[0] < (0.5*cellSizeX_ + eps_) && globalPos[1] < eps_;
    }

    bool isInlet(const GlobalPosition& globalPos) const
    {
        return globalPos[0] < eps_;
    }

    bool isOutlet(const GlobalPosition& globalPos) const
    {
        return globalPos[0] > this->bBoxMax()[0] - eps_;
    }

    bool isWall(const GlobalPosition& globalPos) const
    {
        return globalPos[0] > eps_ || globalPos[0] < this->bBoxMax()[0] - eps_;
    }


    Scalar eps_;
    Scalar lidVelocity_;
    Scalar cellSizeX_;
    Scalar inletVelocity_;
    Scalar kinematicViscosity_;
};
} //end namespace

#endif
