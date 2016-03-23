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
 * \brief A test problem for the one-phase model:
 * water is flowing from bottom to top through and around a low permeable lens.
 */
#ifndef DUMUX_1PTEST_PROBLEM_HH
#define DUMUX_1PTEST_PROBLEM_HH

#include <dumux/porousmediumflow/1p/implicit/model.hh>
#include <dumux/porousmediumflow/implicit/problem.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/fluidsystems/liquidphase.hh>

#include <dumux/linear/amgbackend.hh>

#include "1ptestspatialparams.hh"

namespace Dumux
{
template <class TypeTag>
class OnePTestProblem;

namespace Properties
{
NEW_TYPE_TAG(OnePTestProblem, INHERITS_FROM(OneP));
NEW_TYPE_TAG(OnePTestBoxProblem, INHERITS_FROM(BoxModel, OnePTestProblem));
NEW_TYPE_TAG(OnePTestCCProblem, INHERITS_FROM(CCModel, OnePTestProblem));

SET_PROP(OnePTestProblem, Fluid)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Dumux::FluidSystems::LiquidPhase<Scalar, Dumux::SimpleH2O<Scalar> > type;
};

// Set the grid type
SET_TYPE_PROP(OnePTestProblem, Grid, Dune::YaspGrid<2>);
//SET_TYPE_PROP(OnePTestProblem, Grid, Dune::UGGrid<2>);
//SET_TYPE_PROP(OnePTestProblem, Grid, Dune::ALUGrid<2, 2, Dune::cube, Dune::nonconforming>);

// Set the problem property
SET_TYPE_PROP(OnePTestProblem, Problem, Dumux::OnePTestProblem<TypeTag> );

// Set the spatial parameters
SET_TYPE_PROP(OnePTestProblem, SpatialParams, Dumux::OnePTestSpatialParams<TypeTag> );

// Linear solver settings
SET_TYPE_PROP(OnePTestProblem, LinearSolver, Dumux::ILU0BiCGSTABBackend<TypeTag> );

NEW_TYPE_TAG(OnePTestBoxProblemWithAMG, INHERITS_FROM(OnePTestBoxProblem));
NEW_TYPE_TAG(OnePTestCCProblemWithAMG, INHERITS_FROM(OnePTestCCProblem));
// Solver settings for the tests using AMG
SET_TYPE_PROP(OnePTestBoxProblemWithAMG, LinearSolver, Dumux::AMGBackend<TypeTag> );
SET_TYPE_PROP(OnePTestCCProblemWithAMG, LinearSolver, Dumux::AMGBackend<TypeTag> );

// if FoamGrid is available, test for dim < dimWorld
#if HAVE_DUNE_FOAMGRID
NEW_TYPE_TAG(OnePOneDThreeDTestProblem, INHERITS_FROM(OnePTestProblem));
NEW_TYPE_TAG(OnePOneDThreeDTestBoxProblem, INHERITS_FROM(BoxModel, OnePOneDThreeDTestProblem));
NEW_TYPE_TAG(OnePOneDThreeDTestCCProblem, INHERITS_FROM(CCModel, OnePOneDThreeDTestProblem));
SET_TYPE_PROP(OnePOneDThreeDTestProblem, Grid, Dune::FoamGrid<1, 3>);

NEW_TYPE_TAG(OnePTwoDThreeDTestProblem, INHERITS_FROM(OnePTestProblem));
NEW_TYPE_TAG(OnePTwoDThreeDTestBoxProblem, INHERITS_FROM(BoxModel, OnePTwoDThreeDTestProblem));
NEW_TYPE_TAG(OnePTwoDThreeDTestCCProblem, INHERITS_FROM(CCModel, OnePTwoDThreeDTestProblem));
SET_TYPE_PROP(OnePTwoDThreeDTestProblem, Grid, Dune::FoamGrid<2, 3>);
#endif

// Enable gravity
SET_BOOL_PROP(OnePTestProblem, ProblemEnableGravity, true);
}

/*!
 * \ingroup OnePModel
 * \ingroup ImplicitTestProblems
 * \brief  Test problem for the one-phase model:
 * water is flowing from bottom to top through and around a low permeable lens.
 *
 * The domain is box shaped. All sides are closed (Neumann 0 boundary)
 * except the top and bottom boundaries (Dirichlet), where water is
 * flowing from bottom to top.
 *
 * In the middle of the domain, a lens with low permeability (\f$K=10e-12\f$)
 * compared to the surrounding material (\f$ K=10e-10\f$) is defined.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_box1p -parameterFile test_box1p.input</tt> or
 * <tt>./test_cc1p -parameterFile test_cc1p.input</tt>
 *
 * The same parameter file can be also used for 3d simulation but you need to change line
 * <tt>typedef Dune::YaspGrid<2> type;</tt> to
 * <tt>typedef Dune::YaspGrid<3> type;</tt> in the problem file
 * and use <tt>test_1p_3d.dgf</tt> in the parameter file.
 */
template <class TypeTag>
class OnePTestProblem : public ImplicitPorousMediaProblem<TypeTag>
{
    typedef ImplicitPorousMediaProblem<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    // copy some indices for convenience
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    enum {
        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };
    enum {
        // indices of the primary variables
        conti0EqIdx = Indices::conti0EqIdx,
        pressureIdx = Indices::pressureIdx
    };

    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::Intersection Intersection;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:
    OnePTestProblem(TimeManager &timeManager, const GridView &gridView)
    : ParentType(timeManager, gridView)
    {
        name_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag,
                                             std::string,
                                             Problem,
                                             Name);
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
        return name_.c_str();
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
     * \param values Stores the source values, acts as return value
     * \param globalPos The global position
     */
    void sourceAtPos(PrimaryVariables &values,
                const GlobalPosition &globalPos) const
    {
        values = 0;
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
     * \param values The boundary types for the conservation equations
     * \param globalPos The position of the center of the finite volume
     */
    void boundaryTypesAtPos(BoundaryTypes &values,
            const GlobalPosition &globalPos) const
    {
        double eps = 1.0e-3;
        if (globalPos[dimWorld-1] < eps || globalPos[dimWorld-1] > this->bBoxMax()[dimWorld-1] - eps)
            values.setAllDirichlet();
        else
            values.setAllNeumann();
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume.
     *
     * \param values The dirichlet values for the primary variables
     * \param globalPos The center of the finite volume which ought to be set.
     *
     * For this method, the \a values parameter stores primary variables.
     */
    void dirichletAtPos(PrimaryVariables &values,
                        const GlobalPosition &globalPos) const
    {
        values[pressureIdx] = 1.0e+5*(2.0 - globalPos[dimWorld-1]);
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * For this method, the \a priVars parameter stores the mass flux
     * in normal direction of each component. Negative values mean
     * influx.
     */
    using ParentType::neumann;
    void neumann(PrimaryVariables &priVars,
                 const Element &element,
                 const FVElementGeometry &fvGeometry,
                 const Intersection &intersection,
                 const int scvIdx,
                 const int boundaryFaceIdx) const
    {
        priVars[conti0EqIdx] = 0;
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * For this method, the \a priVars parameter stores primary
     * variables.
     */
    void initial(PrimaryVariables &priVars,
                 const Element &element,
                 const FVElementGeometry &fvGeometry,
                 const int scvIdx) const
    {
        priVars[pressureIdx] = 1.0e+5;
    }

    // \}

private:
    std::string name_;
};
} //end namespace

#endif
