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

#include <dumux/material/fluidsystems/h2on2.hh>
#include <dumux/material/fluidsystems/gasphase.hh>
#include <dumux/material/fluidsystems/liquidphase.hh>
#include <dumux/material/components/constant.hh>
#include <dumux/material/components/simpleh2o.hh>


//added for running purposes, added to resemble elastic
#include <dumux/implicit/fem/properties.hh>
#include <dumux/implicit/fem/problem.hh>
#include <dumux/geomechanics/elastic/model.hh>
#include <dumux/linear/seqsolverbackend.hh>
#include "elasticspatialparams.hh"
//?
#include <dumux/geomechanics/elastic/model.hh>

//variated to stokes_fem from stokes model
#include <dumux/freeflow/stokes_fem/model.hh>

namespace Dumux
{

template <class TypeTag>
class StokesFemTestProblem;

//////////
// Specify the properties for the stokes problem
//////////
namespace Properties
{//valid to inherit from both?
NEW_TYPE_TAG(StokesFemTestProblem, INHERITS_FROM(BoxStokes, FemModel, ElSpatialParams));

// Set the grid type
SET_TYPE_PROP(StokesFemTestProblem, Grid, Dune::YaspGrid<2>);

// Set the problem property
SET_TYPE_PROP(StokesFemTestProblem, Problem, StokesFemTestProblem<TypeTag>);

////// Use nitrogen as gas phase
//SET_TYPE_PROP(StokesFemTestProblem, Fluid,
//              FluidSystems::GasPhase<typename GET_PROP_TYPE(TypeTag, Scalar),
//                                            N2<typename GET_PROP_TYPE(TypeTag, Scalar)> >);

//// Use nitrogen as gas phase
SET_TYPE_PROP(StokesFemTestProblem, Fluid,
              FluidSystems::LiquidPhase<typename GET_PROP_TYPE(TypeTag, Scalar),
                                        SimpleH2O<typename GET_PROP_TYPE(TypeTag, Scalar)> >);

// Use nitrogen as gas phase
//SET_TYPE_PROP(StokesFemTestProblem, Fluid,
//              FluidSystems::LiquidPhase<typename GET_PROP_TYPE(TypeTag, Scalar),
//                                            Constant<TypeTag ,typename GET_PROP_TYPE(TypeTag, Scalar)> >);


//added from elastic
// Quadrature order
SET_INT_PROP(StokesFemTestProblem, FemBasisOrder, 1);

SET_TYPE_PROP(StokesFemTestProblem, LinearSolver, ILUnBiCGSTABBackend<TypeTag>);
}

/*!
 * \ingroup BoxStokesModel
 * \ingroup ImplicitTestProblems
 * \brief Stokes flow problem with nitrogen (N2) flowing
 *        from the left to the right.
 *
 * The domain is sized 1m times 1m. The boundary conditions for the momentum balances
 * are set to Dirichlet with outflow on the right boundary. The mass balance has
 * outflow bcs, which are replaced in the localresidual by the sum
 * of the momentum balance equations in case of Dirichlet bcs for the momentum balance.
 * In the middle of the right boundary, one vertex receives Dirichlet bcs to set the pressure level.
 * The flow velocity starts with 0 m/s. A flow field evolves with a maximum velocity, which is
 * varied time-dependently using a sinus function and a period of 3000s.
 *
 * This problem uses the \ref BoxStokesModel.
 * To run the simulation execute the following line in shell:
 * <tt>./test_stokes</tt>
 */
template <class TypeTag>
class StokesFemTestProblem : public ImplicitFemProblem<TypeTag>
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
    using Vertex = typename GridView::template Codim<dim>::Entity;
    using CoordScalar = typename GridView::ctype;
    using Intersection = typename GridView::Intersection;

    //typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    using Fluid = typename GET_PROP_TYPE(TypeTag, Fluid);
    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;

//added for running purposes
    //from stokes
 //   typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;

public:
    StokesFemTestProblem(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView)
    {
        eps_ = 1e-6;
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief The problem name.
     *
     * This is used as a prefix for files generated by the simulation.
     */
    const char *name() const
    { return "stokes_fem"; }

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
    // \}

    /*!
     * \name Boundary conditions
     */
    // \{


    //DONEA_TEST_PROBLEM
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
        {
            BoundaryTypes values;

            // set Dirichlet values for the velocity and pressure everywhere
            //values.setAllDirichlet();

            values.setAllDirichlet();

            return values;
        }


        //from ELASTIC
        PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
        {
            return analyticalSolution(globalPos);
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

//std::cout << "Koord. x _ y: " << x << "  " << y << std::endl;


            PrimaryVariables values;
            values[pressureIdx] = x * (1.0-x); // p(x,y) = x(1-x) [Donea2003]
            values[velocityXIdx] = x*x * (1.0 - x)*(1.0 - x) * (2.0*y - 6.0*y*y + 4.0*y*y*y);
            values[velocityYIdx] = -1.0*y*y * (1.0 - y)*(1.0 - y) * (2.0*x - 6.0*x*x + 4.0*x*x*x);

//printvector(std::cout, values, "valuesAnalytical","");


            return values;
        }




    //from ELASTIC
//        PrimaryVariables source(const Element& element,
//                                const IpData& ipData,
//                                const SecondaryVariables& secVars) const
//        {
//            const auto ipGlobal = ipData.ipGlobal();
//            const auto x = ipGlobal[0];
//            const auto y = ipGlobal[1];

//            PrimaryVariables source(0.0);
//      /*      source[Indices::momentumXIdx] = (12.0-24.0*y) * x*x*x*x + (-24.0 + 48.0*y)* x*x*x
//                                          + (-48.0*y + 72.0*y*y - 48.0*y*y*y + 12.0)* x*x
//                                          + (-2.0 + 24.0*y - 72.0*y*y + 48.0*y*y*y)*x
//                                          + 1.0 - 4.0*y + 12.0*y*y - 8.0*y*y*y;
//            source[Indices::momentumYIdx] = (8.0 - 48.0*y + 48.0*y*y)*x*x*x + (-12.0 + 72.0*y - 72.0*y*y)*x*x
//                                          + (4.0 - 24.0*y + 48.0*y*y - 48.0*y*y*y + 24.0*y*y*y*y)*x - 12.0*y*y
//                                          + 24.0*y*y*y - 12.0*y*y*y*y;
//    */
//            return source;
//        }

    //from donea
    PrimaryVariables sourceAtPos(const GlobalPosition &globalPos) const
        {
            PrimaryVariables source(0.0);
            Scalar x = globalPos[0];
            Scalar y = globalPos[1];

//std::cout << "Koord. x _ y: " << x << "  " << y << std::endl;


            source[momentumXIdx] = (12.0-24.0*y) * x*x*x*x + (-24.0 + 48.0*y)* x*x*x
                                          + (-48.0*y + 72.0*y*y - 48.0*y*y*y + 12.0)* x*x
                                          + (-2.0 + 24.0*y - 72.0*y*y + 48.0*y*y*y)*x
                                          + 1.0 - 4.0*y + 12.0*y*y - 8.0*y*y*y;
            source[momentumYIdx] = (8.0 - 48.0*y + 48.0*y*y)*x*x*x + (-12.0 + 72.0*y - 72.0*y*y)*x*x
                                          + (4.0 - 24.0*y + 48.0*y*y - 48.0*y*y*y + 24.0*y*y*y*y)*x - 12.0*y*y
                                          + 24.0*y*y*y - 12.0*y*y*y*y;

//printvector(std::cout, source, "sourceTest","");
            return source;
        }


        //from ELASTIC
        PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
        {
            PrimaryVariables values;
            values[velocityXIdx] = 0.0;
            values[velocityYIdx] = 0.0;
            values[pressureIdx] = 0.0;
          //  std::cout << "globalPos=" << globalPos << std::endl;
            return values;
        }


    PrimaryVariables neumannAtPos(const GlobalPosition& globalPos) const
            {
                PrimaryVariables neumannPrim(0.0);
                Scalar x = globalPos[0];
                Scalar y = globalPos[1];

                if(onLeftBoundary_(globalPos)){

            }


                if(onRightBoundary_(globalPos)){
                    neumannPrim[velocityXIdx] = 0.05;
                    neumannPrim[velocityYIdx] = 0;
                    neumannPrim[pressureIdx] = -0.5;
                }


                if(onUpperBoundary_(globalPos)){
                    neumannPrim[pressureIdx] = 0;
                }


                if(onLowerBoundary_(globalPos)){
                    neumannPrim[pressureIdx] = 0;
                }

                //return analyticalSolution(globalPos);
                return neumannPrim;
            }

















    /* CHANNEL_TEST_PROBLEM
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;

        // set Dirichlet values for the velocity and pressure everywhere
        //values.setAllDirichlet();

        values.setAllNeumann();

        if(onLeftBoundary_(globalPos)){
            values.setDirichlet(velocityXIdx,momentumXIdx);
            values.setDirichlet(velocityYIdx,momentumYIdx);
        }


        if(onRightBoundary_(globalPos)){
            values.setNeumann(momentumXIdx);
            values.setNeumann(momentumYIdx);
        }


        if(onUpperBoundary_(globalPos)){
            values.setDirichlet(velocityXIdx,momentumXIdx);
            values.setDirichlet(velocityYIdx,momentumYIdx);
            //p pet setAllNeumann zu 0 gesetzt
        }


        if(onLowerBoundary_(globalPos)){
            values.setDirichlet(velocityXIdx,momentumXIdx);
            values.setDirichlet(velocityYIdx,momentumYIdx);
            //p pet setAllNeumann zu 0 gesetzt
        }


        return values;
    }


    //from ELASTIC
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    {
        PrimaryVariables dirichletPrim(0.0);
        Scalar x = globalPos[0];
        Scalar y = globalPos[1];

        if(onLeftBoundary_(globalPos)){
            //  dirichletPrim[velocityXIdx] = 4/boxY * y - 4/boxY² * y²
            dirichletPrim[velocityXIdx] = 4/3 * y - 4/9 * y*y;
            dirichletPrim[velocityYIdx] = 0;
            dirichletPrim[pressureIdx]  = 1e5;
        }


        if(onRightBoundary_(globalPos)){
        }


        if(onUpperBoundary_(globalPos)){
            dirichletPrim[velocityXIdx] = 0;
            dirichletPrim[velocityYIdx] = 0;
            //p pet setAllNeumann zu 0 gesetzt
        }


        if(onLowerBoundary_(globalPos)){
            dirichletPrim[velocityXIdx] = 0;
            dirichletPrim[velocityYIdx] = 0;
            //p pet setAllNeumann zu 0 gesetzt
        }


            //return analyticalSolution(globalPos);
    }



    PrimaryVariables neumannAtPos(const GlobalPosition& globalPos) const
    {
        PrimaryVariables neumannPrim(0.0);
        Scalar x = globalPos[0];
        Scalar y = globalPos[1];

        if(onLeftBoundary_(globalPos)){

    }


        if(onRightBoundary_(globalPos)){
            neumannPrim[velocityXIdx] = 0.05;
            neumannPrim[velocityYIdx] = 0;
            neumannPrim[pressureIdx] = -0.5;
        }


        if(onUpperBoundary_(globalPos)){
            neumannPrim[pressureIdx] = 0;
        }


        if(onLowerBoundary_(globalPos)){
            neumannPrim[pressureIdx] = 0;
        }

        //return analyticalSolution(globalPos);
    }



    PrimaryVariables analyticalSolution(const GlobalPosition& globalPos) const
    {
        Scalar x = globalPos[0];
        Scalar y = globalPos[1];

        PrimaryVariables values;
        values[pressureIdx] = x * (1.0-x); // p(x,y) = x(1-x) [Donea2003]
        values[velocityXIdx] = x*x * (1.0 - x)*(1.0 - x) * (2.0*y - 6.0*y*y + 4.0*y*y*y);
        values[velocityYIdx] = -1.0*y*y * (1.0 - y)*(1.0 - y) * (2.0*x - 6.0*x*x + 4.0*x*x*x);

        return values;
    }


 //   PrimaryVariables neumann(const Element& element,
 //                            const Intersection& intersection,
 //                            const ElementSolutionVector& elemSol,
 //                            const IpData& ipData) const
 //   { return PrimaryVariables(0.0); }




//from ELASTIC
    PrimaryVariables source(const Element& element,
                            const IpData& ipData,
                            const SecondaryVariables& secVars) const
    {
        const auto ipGlobal = ipData.ipGlobal();
        const auto x = ipGlobal[0];
        const auto y = ipGlobal[1];

        PrimaryVariables source(0.0);
//       source[Indices::momentumXIdx] = (12.0-24.0*y) * x*x*x*x + (-24.0 + 48.0*y)* x*x*x
//                                      + (-48.0*y + 72.0*y*y - 48.0*y*y*y + 12.0)* x*x
//                                      + (-2.0 + 24.0*y - 72.0*y*y + 48.0*y*y*y)*x
//                                      + 1.0 - 4.0*y + 12.0*y*y - 8.0*y*y*y;
//        source[Indices::momentumYIdx] = (8.0 - 48.0*y + 48.0*y*y)*x*x*x + (-12.0 + 72.0*y - 72.0*y*y)*x*x
//                                      + (4.0 - 24.0*y + 48.0*y*y - 48.0*y*y*y + 24.0*y*y*y*y)*x - 12.0*y*y
//                                      + 24.0*y*y*y - 12.0*y*y*y*y;

        return source;
    }


    //from ELASTIC
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    {
        PrimaryVariables values;
        values[pressureIdx] = 1e5;
        values[velocityXIdx] = 1.0;
        values[velocityYIdx] = 0.0;
      //  std::cout << "globalPos=" << globalPos << std::endl;
        return values;
    }
*/

private:

    bool onLeftBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] < this->bBoxMin()[0] + eps_; }

    bool onRightBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] > this->bBoxMax()[0] - eps_; }

    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] < this->bBoxMin()[1] + eps_; }

    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] > this->bBoxMax()[1] - eps_; }

    Scalar eps_;

};

} //end namespace

#endif
