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
 * \file
 * \brief Definition of an isothermal ZeroEq problem.
 *
 * This problem implements a pipe flow experiment performed by John Laufer
 * (Laufer J., The structure of turbulence in fully developed pipe flow, NACA Report, 1954).
 */
#ifndef DUMUX_ZEROEQTESTPROBLEM_HH
#define DUMUX_ZEROEQTESTPROBLEM_HH

#include <dumux/material/fluidsystems/h2oair.hh>
#include <dumux/material/fluidsystems/gasphase.hh>

#include <dumux/freeflow/zeroeq/model.hh>
#include <dumux/freeflow/zeroeq/problem.hh>

namespace Dumux
{

template <class TypeTag>
class ZeroEqTestProblem;

namespace Properties
{
NEW_TYPE_TAG(ZeroEqTestProblem, INHERITS_FROM(BoxZeroEq));

// Set the grid type
SET_TYPE_PROP(ZeroEqTestProblem, Grid, Dune::YaspGrid<2>);

// Set the problem property
SET_TYPE_PROP(ZeroEqTestProblem, Problem, ZeroEqTestProblem<TypeTag>);

// Set the air as the gas phase
SET_TYPE_PROP(ZeroEqTestProblem, Fluid,
              FluidSystems::GasPhase<typename GET_PROP_TYPE(TypeTag, Scalar),
                                            Air<typename GET_PROP_TYPE(TypeTag, Scalar)> >);

// Disable gravity
SET_BOOL_PROP(ZeroEqTestProblem, ProblemEnableGravity, false);
}

/*!
 * \ingroup BoxZeroEqModel
 * \ingroup ImplicitTestProblems
 * \brief ZeroEq problem with air flowing from the left to the right.
 *
 * The domain is sized 10m times 0.2469m. The problem is taken from an
 * experimental setup by John Laufer (J. Laufer, The structure of turbulence in
 * fully developed pipe flow, NACA Report, 1954).
 * The boundary conditions for the momentum balances
 * are set to Dirichlet on the left (inflow) and outflow on the right boundary.
 * The mass balance has outflow boundary conditions, which are replaced in the
 * localresidual by the sum of the two momentum balances. On the right boundary,
 * the mass balance receives a Dirichlet value to set the pressure level.
 *
 * This problem uses the \ref ZeroEqModel with the Baldwin-Lomax turbulence model
 * (Baldwin, B. S. & Lomax, H. Thin Layer Approximation and Algebraic Model for
 * Seperated Turbulent Flows AIAA Journal, 1978).
 *
 * To run the simulation execute the following line in shell:<br>
 * <tt>./test_zeroeq -ParameterFile ./test_zeroeq.input</tt>
 */
template <class TypeTag>
class ZeroEqTestProblem : public ZeroEqProblem<TypeTag>
{
    typedef ZeroEqProblem<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    // Number of equations and grid dimension
    enum { dim = GridView::dimension };
    enum { // equation indices
        massBalanceIdx = Indices::massBalanceIdx,
        momentumXIdx = Indices::momentumXIdx, // Index of the x-component of the momentum balance
        momentumYIdx = Indices::momentumYIdx // Index of the y-component of the momentum balance
    };
    enum { // indices for primary variables
        velocityXIdx = Indices::velocityXIdx,
        velocityYIdx = Indices::velocityYIdx,
        pressureIdx = Indices::pressureIdx
    };

    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::ctype CoordScalar;
    typedef typename GridView::Intersection Intersection;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, Fluid) Fluid;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dune::FieldVector<CoordScalar, dim> GlobalPosition;

public:
    ZeroEqTestProblem(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView)
    {
        eps_ = 1e-6;

        injectionVelocity_ = GET_RUNTIME_PARAM(TypeTag, Scalar, Problem.InjectionVelocity);
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
    const std::string &name() const
    {
        return GET_RUNTIME_PARAM(TypeTag, std::string, Problem.OutputName);
    }

    /*!
     * \brief Returns the temperature within the domain.
     *
     * This problem assumes a constant temperature of 10 degrees Celsius.
     */
    Scalar temperatureAtPos(const GlobalPosition &globalPos) const
    {
        return 273.15 + 10; // -> 10C
    }

    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    //! \copydoc ImplicitProblem::boundaryTypesAtPos()
    void boundaryTypesAtPos(BoundaryTypes &values,
                            const GlobalPosition &globalPos) const
    {
        values.setAllDirichlet();

        // the mass balance has to be of type outflow
        values.setOutflow(massBalanceIdx);

        if (onRightBoundary_(globalPos)
            && globalPos[1] < this->bBoxMax()[1]-eps_ && globalPos[1] > this->bBoxMin()[1]+eps_)
            values.setAllOutflow();

        // set pressure on left boundary (at least at one point)
        if (onRightBoundary_(globalPos))
            values.setDirichlet(pressureIdx);
    }

    //! \copydoc ImplicitProblem::dirichletAtPos()
    void dirichletAtPos(PrimaryVariables &values,
                        const GlobalPosition &globalPos) const
    {
        initial_(values, globalPos);
    }

    /*!
     * \copydoc ImplicitProblem::neumann()
     * A neumann condition for the RANS momentum equation equation corresponds to:
     * \f[ - \left[ \mu + \mu_\textrm{t} \right] \nabla {\bf v} \cdot {\bf n} + p \cdot {\bf n} = q_N \f]
     */
    void neumann(PrimaryVariables &values,
                 const Element &element,
                 const FVElementGeometry &fvGeometry,
                 const Intersection &intersection,
                 int scvIdx,
                 int boundaryFaceIdx) const
    {
        values = 0.0;
    }
    // \}

    /*!
     * \name Volume terms
     */
    // \{

    //! \copydoc ImplicitProblem::sourceAtPos()
    void sourceAtPos(PrimaryVariables &values,
                     const GlobalPosition &globalPos) const
    {
        values = 0.0;
    }

    //! \copydoc ImplicitProblem::initialAtPos()
    void initialAtPos(PrimaryVariables &values,
                      const GlobalPosition &globalPos) const
    {
        initial_(values, globalPos);
    }
    // \}

private:
    // internal method for the initial condition (reused for the Dirichlet condition)
    void initial_(PrimaryVariables &values, const GlobalPosition &globalPos) const
    {
        values = 0.0;

        if (!(onLowerBoundary_(globalPos) || onUpperBoundary_(globalPos)))
        {
            values[velocityXIdx] = injectionVelocity_;
        }

        values[pressureIdx] = 1e5;
    }

    bool onLeftBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] < this->bBoxMin()[0] + eps_; }

    bool onRightBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] > this->bBoxMax()[0] - eps_; }

    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] < this->bBoxMin()[1] + eps_; }

    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] > this->bBoxMax()[1] - eps_; }

    Scalar eps_;
    Scalar injectionVelocity_;
};

} //end namespace

#endif // DUMUX_ZEROEQTESTPROBLEM_HH
