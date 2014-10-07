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
 *
 * \brief test problem for the decoupled one-phase model.
 */
#ifndef DUMUX_TEST_1P_PROBLEM_HH
#define DUMUX_TEST_1P_PROBLEM_HH

#if HAVE_UG
#include <dune/grid/uggrid.hh>
#endif

#include <dune/grid/yaspgrid.hh>
#include <dumux/io/cubegridcreator.hh>

#include <dumux/material/fluidsystems/liquidphase.hh>
#include <dumux/material/components/unit.hh>

#include <dumux/decoupled/1p/diffusion/fv/fvpressureproperties1p.hh>
#include <dumux/decoupled/1p/diffusion/diffusionproblem1p.hh>
#include <dumux/decoupled/common/fv/fvvelocity.hh>

#include "test_1pspatialparams.hh"

namespace Dumux
{

template<class TypeTag>
class TestProblemOneP;

//////////
// Specify the properties
//////////
namespace Properties
{
NEW_TYPE_TAG(TestProblemOneP, INHERITS_FROM(FVPressureOneP));

// set the GridCreator property
SET_TYPE_PROP(TestProblemOneP, GridCreator, CubeGridCreator<TypeTag>);

// Set the grid type
SET_PROP(TestProblemOneP, Grid)
{
        typedef Dune::YaspGrid<2> type;
};

// Set the wetting phase
SET_PROP(TestProblemOneP, Fluid)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::Unit<Scalar> > type;
};

// Set the spatial parameters
SET_TYPE_PROP(TestProblemOneP, SpatialParams, Dumux::TestOnePSpatialParams<TypeTag>);

// Enable gravity
SET_BOOL_PROP(TestProblemOneP, ProblemEnableGravity, false);

//Set the problem
SET_TYPE_PROP(TestProblemOneP, Problem, Dumux::TestProblemOneP<TypeTag>);


SET_INT_PROP(TestProblemOneP, LinearSolverVerbosity, 1);
}

/*!
 * \ingroup IMPETtests
 *
 * \brief test problem for the decoupled one-phase model.
 */
template<class TypeTag>
class TestProblemOneP: public DiffusionProblem1P<TypeTag >
{
    typedef DiffusionProblem1P<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, Fluid) Fluid;

    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;

    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::Intersection Intersection;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldVector<Scalar, dim> LocalPosition;
    typedef typename GET_PROP(TypeTag, ParameterTree) ParameterTree;
    typedef typename GET_PROP_TYPE(TypeTag, GridCreator) GridCreator;


public:
    TestProblemOneP(TimeManager &timeManager, const GridView &gridView) :
        ParentType(gridView), velocity_(*this)
    {
        delta_ = 1e-3 ;

        try
        {
            if (ParameterTree::tree().hasKey("Problem.Delta"))
                delta_       = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, Delta);
            int numRefine;
            numRefine = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, Grid, NumRefine);
            GridCreator::grid().globalRefine(numRefine);
        }
        catch (Dumux::ParameterException &e) {
            std::cerr << e << ". Abort!\n";
            exit(1) ;
        }
        catch (...) {
            std::cerr << "Unknown exception thrown!\n";
            exit(1);
        }

        this->spatialParams().initialize(delta_);
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
    {
        return "test_1p";
    }

    bool shouldWriteRestartFile() const
    { return false; }

    void addOutputVtkFields()
    {
        velocity_.calculateVelocity();
        velocity_.addOutputVtkFields(this->resultWriter());
    }

    /*!
    * \brief Returns the temperature within the domain.
    *
    * This problem assumes a temperature of 10 degrees Celsius.
    */
    Scalar temperatureAtPos(const GlobalPosition& globalPos) const
    {
        return 273.15 + 10; // -> 10°C
    }

    // \}

    //! Returns the reference pressure for evaluation of constitutive relations
    Scalar referencePressureAtPos(const GlobalPosition& globalPos) const
    {
        return 1e5; // -> 10°C
    }

    //!source term [kg/(m^3 s)]
    void source(PrimaryVariables &values, const Element& element) const
        {
        values = 0;

        values = integratedSource_(element, 4);
        }

    /*!
    * \brief Returns the type of boundary condition.
    *
    * BC can be dirichlet (pressure) or neumann (flux).
    */
    void boundaryTypes(BoundaryTypes &bcType,
            const Intersection& intersection) const
    {
        bcType.setAllDirichlet();
    }

    //! return dirichlet condition  (pressure, [Pa])
    void dirichletAtPos(PrimaryVariables &values,
                        const GlobalPosition &globalPos) const
    {
        values = exact(globalPos);
    }


    //! return neumann condition  (flux, [kg/(m^2 s)])
    void neumann(PrimaryVariables &values, const Intersection& intersection) const
        {
        values = 0;
        }

private:
    Scalar exact (const GlobalPosition& globalPos) const
    {
        double pi = 4.0*atan(1.0);

        return (sin(pi*globalPos[0])*sin(pi*globalPos[1]));
    }

    Dune::FieldVector<Scalar,dim> exactGrad (const GlobalPosition& globalPos) const
        {
        Dune::FieldVector<Scalar,dim> grad(0);
        double pi = 4.0*atan(1.0);
        grad[0] = pi*cos(pi*globalPos[0])*sin(pi*globalPos[1]);
        grad[1] = pi*cos(pi*globalPos[1])*sin(pi*globalPos[0]);

        return grad;
        }

    Scalar integratedSource_(const Element& element, int integrationPoints) const
    {
        Scalar source = 0.;
        LocalPosition localPos(0.0);
        GlobalPosition globalPos(0.0);
        Scalar halfInterval = 1.0/double(integrationPoints)/2.;
        for (int i = 1; i <= integrationPoints; i++)
        {
            for (int j = 1; j <= integrationPoints; j++)
            {
                localPos[0] = double(i)/double(integrationPoints) - halfInterval;
                localPos[1] = double(j)/double(integrationPoints) - halfInterval;
                globalPos = element.geometry().global(localPos);
                source += 1./(integrationPoints*integrationPoints) * evaluateSource_(globalPos);
            }
        }

        return source;
    }

    Scalar evaluateSource_(const GlobalPosition& globalPos) const
    {
        Scalar temp = temperatureAtPos(globalPos);
        Scalar referencePress = referencePressureAtPos(globalPos);

        Scalar pi = 4.0 * atan(1.0);
        Scalar x = globalPos[0];
        Scalar y = globalPos[1];

        Scalar dpdx = pi * cos(pi * x) * sin(pi * y);
        Scalar dpdy = pi * sin(pi * x) * cos(pi * y);
        Scalar dppdxx = -pi * pi * sin(pi * x) * sin(pi * y);
        Scalar dppdxy = pi * pi * cos(pi * x) * cos(pi * y);
        Scalar dppdyx = dppdxy;
        Scalar dppdyy = dppdxx;
        Scalar kxx = (delta_* x*x + y*y)/(x*x + y*y);
        Scalar kxy = -(1.0 - delta_) * x * y / (x*x + y*y);
        Scalar kyy = (x*x + delta_*y*y)/(x*x + y*y);
        Scalar dkxxdx = 2 * x * y*y * (delta_ - 1.0)/((x*x + y*y) * (x*x + y*y));
        Scalar dkyydy = 2 * x*x * y * (delta_ - 1.0)/((x*x + y*y) * (x*x + y*y));
        Scalar dkxydx = (1.0 - delta_) * y * (x*x - y*y) /((x*x + y*y) * (x*x + y*y));
        Scalar dkxydy = (1.0 - delta_) * x * (y*y - x*x) /((x*x + y*y) * (x*x + y*y));

        Scalar fx = dkxxdx * dpdx + kxx * dppdxx + dkxydx * dpdy + kxy * dppdyx;
        Scalar fy = dkxydy * dpdx + kxy * dppdxy + dkyydy * dpdy + kyy * dppdyy;

        return -(fx + fy) / Fluid::viscosity(temp, referencePress) * Fluid::density(temp, referencePress);
    }

    double delta_;
    Dumux::FVVelocity<TypeTag, typename GET_PROP_TYPE(TypeTag, Velocity) > velocity_;
};
} //end namespace

#endif
