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
 * \brief  Definition of an isothermal ZeroEq channel flow problem.
 */
#ifndef DUMUX_ZEROEQCHANNELTESTPROBLEM_HH
#define DUMUX_ZEROEQCHANNELTESTPROBLEM_HH

#include <dumux/material/fluidsystems/h2oair.hh>
#include <dumux/material/fluidsystems/gasphase.hh>

#include <dumux/freeflow/zeroeq/model.hh>
#include <dumux/freeflow/zeroeq/problem.hh>

namespace Dumux
{

template <class TypeTag>
class ZeroEqChannelTestProblem;

namespace Properties
{
NEW_TYPE_TAG(ZeroEqChannelTestProblem, INHERITS_FROM(BoxZeroEq));

// Set the grid type
SET_TYPE_PROP(ZeroEqChannelTestProblem, Grid, Dune::YaspGrid<2>);

// Set the problem property
SET_TYPE_PROP(ZeroEqChannelTestProblem, Problem, Dumux::ZeroEqChannelTestProblem<TypeTag>);

// Set the air as the gas phase
SET_TYPE_PROP(ZeroEqChannelTestProblem, Fluid,
              Dumux::FluidSystems::GasPhase<typename GET_PROP_TYPE(TypeTag, Scalar),
                                            Dumux::Air<typename GET_PROP_TYPE(TypeTag, Scalar)> >);

// Disable gravity
SET_BOOL_PROP(ZeroEqChannelTestProblem, ProblemEnableGravity, false);

// Set only bottom as wall
SET_BOOL_PROP(ZeroEqChannelTestProblem, BBoxMaxIsWall, false);

#if HAVE_UMFPACK
// Use UMFPack as linear solver
SET_TYPE_PROP(ZeroEqChannelTestProblem, LinearSolver, UMFPackBackend<TypeTag>);
#endif
}

/*!
 * \ingroup BoxZeroEqModel
 * \ingroup ImplicitTestProblems
 * \brief ZeroEq problem with air flowing from the left to the right.
 *
 * The domain is 5.0m long and 1.5m high. Air is flow from left to right.
 * The momentum balance has Dirichlet conditions on the left (inflow) and
 * on bottom (no-slip), on the right it has outflow condition. On the top,
 * which corresponds to an open surface or the center line of a pipe Neumann
 * no-flow for tangential momentum are applied. The mass balance is outflow
 * except at the right side.
 *
 * This problem uses the \ref ZeroEqModel with the modified Van Driest turbulence model.
 *
 * To run the simulation execute the following line in shell:<br>
 * <tt>./test_zeroeq_channel -ParameterFile ./test_zeroeq_channel.input</tt>
 */
template <class TypeTag>
class ZeroEqChannelTestProblem : public ZeroEqProblem<TypeTag>
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
    ZeroEqChannelTestProblem(TimeManager &timeManager, const GridView &gridView)
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
        return 273.15 + 10; // -> 10 °C
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

        if (onUpperBoundary_(globalPos))
        {
            values.setNeumann(momentumXIdx);
            values.setDirichlet(velocityYIdx);
        }

        if (onRightBoundary_(globalPos))
            values.setAllOutflow();

        if (onLeftBoundary_(globalPos))
            values.setAllDirichlet();

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

        if (!onLowerBoundary_(globalPos))
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

#endif // DUMUX_ZEROEQCHANNELTESTPROBLEM_HH
