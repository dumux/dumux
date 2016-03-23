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
 * \brief  Definition of an isothermal compositional ZeroEqnc problem
 */
#ifndef DUMUX_ZEROEQTWOCTESTPROBLEM_HH
#define DUMUX_ZEROEQTWOCTESTPROBLEM_HH

#include <dumux/material/fluidsystems/h2oair.hh>
#include <dumux/material/fluidsystems/gasphase.hh>

#include <dumux/freeflow/zeroeqnc/model.hh>
#include <dumux/freeflow/zeroeq/problem.hh>

namespace Dumux
{

template <class TypeTag>
class ZeroEq2cTestProblem;

namespace Properties
{
NEW_TYPE_TAG(ZeroEq2cTestProblem, INHERITS_FROM(BoxZeroEqnc));

// Set the grid type
SET_TYPE_PROP(ZeroEq2cTestProblem, Grid, Dune::YaspGrid<2>);

// Set the problem property
SET_TYPE_PROP(ZeroEq2cTestProblem, Problem, Dumux::ZeroEq2cTestProblem<TypeTag>);

// Select the fluid system
SET_TYPE_PROP(ZeroEq2cTestProblem, FluidSystem,
              Dumux::FluidSystems::H2OAir<typename GET_PROP_TYPE(TypeTag, Scalar)>);

// Disable gravity
SET_BOOL_PROP(ZeroEq2cTestProblem, ProblemEnableGravity, false);
}

/*!
 * \ingroup BoxZeroEqnCModel
 * \ingroup ImplicitTestProblems
 * \brief isothermal compositional ZeroEqnc flow problem with injection from the left.
 *
 * The domain is sized 3m times 1m. Air is flowing in a pipe from left to right.
 * Thus the momentum balance has Dirichlet conditions (inflow) on the left and
 * no-slip on top and on bottom, on the right outflow conditions are applied. The total
 * mass balance has outflow boundary conditions everywhere, except at the right where
 * Dirichlet values are set. Water vapor is injected in the middle of the left border
 * and is mixed and transported in the air stream.
 *
 * This problem uses the \ref ZeroEqncModel with the Prandtl mixing length model for
 * the eddy viscosity and the model by Meier and Rotta for the eddy diffusivity.
 *
 * To run the simulation execute the following line in shell:<br>
 * <tt>./test_zeroeq2c -ParameterFile ./test_zeroeq2c.input</tt>
 */
template <class TypeTag>
class ZeroEq2cTestProblem : public ZeroEqProblem<TypeTag>
{
    typedef ZeroEqProblem<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    // Number of equations and grid dimension
    enum { dim = GridView::dimension };
    enum { // equation indices
        massBalanceIdx = Indices::massBalanceIdx,
        momentumXIdx = Indices::momentumXIdx, // Index of the x-component of the momentum balance
        momentumYIdx = Indices::momentumYIdx, // Index of the y-component of the momentum balance
        transportEqIdx = Indices::transportEqIdx  // Index of the transport equation
    };
    enum { // indices for primary variables
        velocityXIdx = Indices::velocityXIdx,
        velocityYIdx = Indices::velocityYIdx,
        pressureIdx = Indices::pressureIdx,
        massOrMoleFracIdx = Indices::massOrMoleFracIdx
    };

    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::ctype CoordScalar;
    typedef typename GridView::Intersection Intersection;
    typedef Dune::FieldVector<CoordScalar, dim> GlobalPosition;

public:
    ZeroEq2cTestProblem(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView)
    {
        eps_ = 1e-6;

        injectionVelocity_ = GET_RUNTIME_PARAM(TypeTag, Scalar, Problem.InjectionVelocity);
        injectionConcentration_ = GET_RUNTIME_PARAM(TypeTag, Scalar, Problem.InjectionConcentration);

        // initialize the tables of the fluid system
        FluidSystem::init();
    }

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
    Scalar temperature() const
    {
        return 273.15 + 10; // -> 10C
    }

    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \name Boundary conditions
     */
    // \{

    //! \copydoc ImplicitProblem::boundaryTypesAtPos()
    void boundaryTypesAtPos(BoundaryTypes &values,
                            const GlobalPosition &globalPos) const
    {
        values.setAllDirichlet();

        if (onRightBoundary_(globalPos)
            && globalPos[1] < this->bBoxMax()[1]-eps_ && globalPos[1] > this->bBoxMin()[1]+eps_)
            values.setAllOutflow();

        // the mass balance has to be of type outflow
        values.setOutflow(massBalanceIdx);

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
            values[velocityXIdx] = injectionVelocity_;

        Scalar middle = 0.5 * (this->bBoxMax()[1] - this->bBoxMin()[1]);
        if (onLeftBoundary_(globalPos) && globalPos[1] > middle - 0.2 && globalPos[1] < middle + 0.2)
            values[massOrMoleFracIdx] = injectionConcentration_;

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

    Scalar injectionVelocity_;
    Scalar injectionConcentration_;
    Scalar eps_;
};

} //end namespace

#endif // DUMUX_ZEROEQTWOCTESTPROBLEM_HH
