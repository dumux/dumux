/*****************************************************************************
*   Copyright (C) 2007-2008 by Klaus Mosthaf                                *
*   Copyright (C) 2007-2008 by Bernd Flemisch                               *
*   Copyright (C) 2008-2009 by Andreas Lauser                               *
*   Institute of Hydraulic Engineering                                      *
*   University of Stuttgart, Germany                                        *
*   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
*                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
*****************************************************************************/
/*!
 * \file
 *
 * \brief test problem for diffusion models from the FVCA5 benchmark.
 */
#ifndef DUMUX_TEST_2P_PROBLEM_HH
#define DUMUX_TEST_2P_PROBLEM_HH

#if HAVE_UG
#include <dune/grid/uggrid.hh>
#endif

#include <dune/grid/yaspgrid.hh>
#include <dune/grid/sgrid.hh>

#include <dumux/material/components/unit.hh>

#include <dumux/decoupled/2p/diffusion/fvmpfa/mpfaproperties.hh>
#include <dumux/decoupled/2p/diffusion/diffusionproblem2p.hh>
#include <dumux/decoupled/2p/diffusion/fv/fvvelocity2p.hh>
#include <dumux/decoupled/2p/diffusion/fvmpfa/fvmpfaovelocity2p.hh>
#include <dumux/decoupled/2p/diffusion/mimetic/mimeticpressure2p.hh>

#include "test_diffusion_spatialparams.hh"

namespace Dumux
{
/*!
 * \ingroup IMPETtests
 */
template<class TypeTag>
class TestDiffusionProblem;

//////////
// Specify the properties
//////////
namespace Properties
{
NEW_TYPE_TAG(DiffusionTestProblem, INHERITS_FROM(DecoupledTwoP, MPFAProperties, TestDiffusionSpatialParams));

// Set the grid type
SET_PROP(DiffusionTestProblem, Grid)
{//    typedef Dune::YaspGrid<2> type;
    typedef Dune::SGrid<2, 2> type;
};

SET_INT_PROP(DiffusionTestProblem, Formulation,
        DecoupledTwoPCommonIndices::pGlobalSw);

// Set the wetting phase
SET_PROP(DiffusionTestProblem, WettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::Unit<Scalar> > type;
};

// Set the non-wetting phase
SET_PROP(DiffusionTestProblem, NonwettingPhase)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
public:
    typedef Dumux::LiquidPhase<Scalar, Dumux::Unit<Scalar> > type;
};

// Enable gravity
SET_BOOL_PROP(DiffusionTestProblem, EnableGravity, false);

// set the types for the 2PFA FV method
NEW_TYPE_TAG(FVVelocity2PTestProblem, INHERITS_FROM(DiffusionTestProblem));
SET_TYPE_PROP(FVVelocity2PTestProblem, Model, Dumux::FVVelocity2P<TTAG(FVVelocity2PTestProblem)>);
SET_TYPE_PROP(FVVelocity2PTestProblem, Problem, Dumux::TestDiffusionProblem<TTAG(FVVelocity2PTestProblem)>);

// set the types for the MPFA-O FV method
NEW_TYPE_TAG(FVMPFAOVelocity2PTestProblem, INHERITS_FROM(DiffusionTestProblem));
SET_TYPE_PROP(FVMPFAOVelocity2PTestProblem, LinearSolver, Dumux::ILUnBiCGSTABBackend<TypeTag>);
SET_INT_PROP(FVMPFAOVelocity2PTestProblem, PreconditionerIterations, 1);
SET_TYPE_PROP(FVMPFAOVelocity2PTestProblem, Model, Dumux::FVMPFAOVelocity2P<TypeTag>);
SET_TYPE_PROP(FVMPFAOVelocity2PTestProblem, Problem, Dumux::TestDiffusionProblem<TTAG(FVMPFAOVelocity2PTestProblem)>);

// set the types for the mimetic FD method
NEW_TYPE_TAG(MimeticPressure2PTestProblem, INHERITS_FROM(DiffusionTestProblem, Mimetic));
SET_TYPE_PROP(MimeticPressure2PTestProblem, Model, Dumux::MimeticPressure2P<TTAG(MimeticPressure2PTestProblem)>);
SET_TYPE_PROP(MimeticPressure2PTestProblem, Problem, Dumux::TestDiffusionProblem<TTAG(MimeticPressure2PTestProblem)>);
}

/*!
 * \ingroup DecoupledProblems
 *
 * \brief test problem for diffusion models from the FVCA5 benchmark.
 */
template<class TypeTag = TTAG(DiffusionTestProblem)>
class TestDiffusionProblem: public DiffusionProblem2P<TypeTag>
{
    typedef DiffusionProblem2P<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPIndices)) Indices;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidSystem)) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FluidState)) FluidState;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(BoundaryTypes)) BoundaryTypes;

    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };

    enum
    {
        wPhaseIdx = Indices::wPhaseIdx,
        nPhaseIdx = Indices::nPhaseIdx,
        pGlobalIdx = Indices::pGlobalIdx,
        SwIdx = Indices::SwIdx,
        pressEqIdx = Indices::pressEqIdx,
    };

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;

    typedef typename GridView::Traits::template Codim<0>::Entity Element;
    typedef typename GridView::Intersection Intersection;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:
    typedef typename GET_PROP(TypeTag, PTAG(SolutionTypes))::PrimaryVariables PrimaryVariables;

    TestDiffusionProblem(const GridView &gridView, const double delta = 1.0) :
        ParentType(gridView), delta_(delta)
    {
        this->variables().saturation() = 1.0;
        this->spatialParameters().setDelta(delta_);
    }

    /*!
    * \name Problem parameters
    */
    // \{

    bool shouldWriteRestartFile() const
    { return false; }

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

    void sourceAtPos(PrimaryVariables &values,const GlobalPosition& globalPos) const
    {
        double pi = 4.0*atan(1.0);
        double rt = globalPos[0]*globalPos[0]+globalPos[1]*globalPos[1];
        double ux = pi*cos(pi*globalPos[0])*sin(pi*globalPos[1]);
        double uy = pi*cos(pi*globalPos[1])*sin(pi*globalPos[0]);
        double kxx = (delta_*globalPos[0]*globalPos[0] + globalPos[1]*globalPos[1])/rt;
        double kxy = -(1.0 - delta_)*globalPos[0]*globalPos[1]/rt;
        double kyy = (globalPos[0]*globalPos[0] + delta_*globalPos[1]*globalPos[1])/rt;
        double f0 = sin(pi*globalPos[0])*sin(pi*globalPos[1])*pi*pi*(1.0 + delta_)*(globalPos[0]*globalPos[0] + globalPos[1]*globalPos[1])
        + cos(pi*globalPos[0])*sin(pi*globalPos[1])*pi*(1.0 - 3.0*delta_)*globalPos[0]
                                                                                    + cos(pi*globalPos[1])*sin(pi*globalPos[0])*pi*(1.0 - 3.0*delta_)*globalPos[1] + cos(pi*globalPos[1])*cos(pi*globalPos[0])*2.0*pi*pi*(1.0 - delta_)*globalPos[0]*globalPos[1];

        values[wPhaseIdx] = (f0 + 2.0*(globalPos[0]*(kxx*ux + kxy*uy) + globalPos[1]*(kxy*ux + kyy*uy)))/rt;
    }

    /*!
    * \brief Returns the type of boundary condition.
    *
    *
    * BC for saturation equation can be dirichlet (saturation), neumann (flux), or outflow.
    */
    void boundaryTypesAtPos(BoundaryTypes &bcTypes, const GlobalPosition& globalPos) const
    {

                bcTypes.setAllDirichlet();
    }

    //! set dirichlet condition  (saturation [-])
    void dirichletAtPos(PrimaryVariables &values, const GlobalPosition& globalPos) const
    {
            values[pGlobalIdx] = exact(globalPos);
            values[SwIdx] = 1.0;
    }

    //! set neumann condition for phases (flux, [kg/(m^2 s)])
    void neumannAtPos(PrimaryVariables &values, const GlobalPosition& globalPos) const
    {
        values = 0;
    }

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

private:
    double delta_;
};
} //end namespace

#endif
