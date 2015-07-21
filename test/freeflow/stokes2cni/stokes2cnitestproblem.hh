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
#ifndef DUMUX_STOKES2CNITESTPROBLEM_HH
#define DUMUX_STOKES2CNITESTPROBLEM_HH

#include <dune/grid/io/file/dgfparser/dgfyasp.hh>

#if HAVE_PARDISO
#include <dumux/linear/pardisobackend.hh>
#endif
#include <dumux/linear/seqsolverbackend.hh>

#include <dumux/material/fluidsystems/h2oairfluidsystem.hh>
#include <dumux/freeflow/stokesncni/stokesncnimodel.hh>

namespace Dumux
{

template <class TypeTag>
class Stokes2cniTestProblem;

//////////
// Specify the properties for the stokes2cni problem
//////////
namespace Properties
{
NEW_TYPE_TAG(Stokes2cniTestProblem, INHERITS_FROM(BoxStokesncni));

// Set the grid type
SET_TYPE_PROP(Stokes2cniTestProblem, Grid, Dune::YaspGrid<2>);

// Set the problem property
SET_TYPE_PROP(Stokes2cniTestProblem, Problem, Stokes2cniTestProblem<TypeTag>);

// Select the fluid system
SET_TYPE_PROP(Stokes2cniTestProblem, FluidSystem,
              Dumux::FluidSystems::H2OAir<typename GET_PROP_TYPE(TypeTag, Scalar)>);

// Use Pardiso as linear solver, if available
#if HAVE_PARDISO
SET_TYPE_PROP(Stokes2cniTestProblem, LinearSolver, PardisoBackend<TypeTag>);
#elif HAVE_UMFPACK
SET_TYPE_PROP(Stokes2cniTestProblem, LinearSolver, UMFPackBackend<TypeTag>);
#endif
}

/*!
 * \ingroup BoxStokesncniModel
 * \ingroup ImplicitTestProblems
 * \brief Stokesncni problem with air flowing
 *        from the bottom to the top, blowing away a warm and dry square.
 *
 * The domain is sized 1m times 1m. An air flow from the bottom boundary blows a warm and dry square,
 * which is initially in the center of the domain out of the upper boundary.
 * The boundary conditions for the momentum balances
 * are all set to Dirichlet. The mass balance has outflow boundary conditions, which are
 * replaced in the localresidual by the sum of the two momentum balances equations in case of Dirichlet bcs for the momentum balance.
 * On the upper boundary a Dirichlet condition is set for the mass balance to fix the pressure.
 * Gravity is on in this example.
 *
 * This problem uses the \ref StokesncniModel.
 * To run the simulation execute the following line in shell:
 * <tt>./test_stokes2cni  -parameterFile ./test_stokes2cni.input</tt>
 */
template <class TypeTag>
class Stokes2cniTestProblem : public StokesProblem<TypeTag>
{
    typedef StokesProblem<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    enum { // grid dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };
    enum { // copy some indices for convenience
        massBalanceIdx = Indices::massBalanceIdx,
        momentumXIdx = Indices::momentumXIdx, //!< Index of the x-component of the momentum balance
        momentumYIdx = Indices::momentumYIdx, //!< Index of the y-component of the momentum balance
        momentumZIdx = Indices::momentumZIdx, //!< Index of the z-component of the momentum balance
        transportEqIdx = Indices::transportEqIdx, //!< Index of the transport equation
    };
    enum { // indices for primary variables
        velocityXIdx = Indices::velocityXIdx,
        velocityYIdx = Indices::velocityYIdx,
        pressureIdx = Indices::pressureIdx,
        massOrMoleFracIdx = Indices::massOrMoleFracIdx,
        temperatureIdx = Indices::temperatureIdx
    };
    enum { //component indices
        transportCompIdx = Indices::transportCompIdx,
        phaseCompIdx = Indices::phaseCompIdx
    };

    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::ctype CoordScalar;
    typedef typename GridView::Intersection Intersection;
    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;

    static const bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);

public:
    Stokes2cniTestProblem(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView)
    {
        eps_ = 1e-6;

        // initialize the tables of the fluid system
        FluidSystem::init();

        velocity_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, Velocity);
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
    { return GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Problem, Name); }

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

        if (onUpperBoundary_(globalPos)
            && !onLeftBoundary_(globalPos)
            && !onRightBoundary_(globalPos)
            && velocity_ > 1e-10 /*if not conduction case*/)
            values.setAllOutflow();

        // set pressure at the upper boundary
        if (onUpperBoundary_(globalPos) &&
                !onLeftBoundary_(globalPos) && !onRightBoundary_(globalPos))
            values.setDirichlet(pressureIdx);
    }

    //! \copydoc ImplicitProblem::dirichletAtPos()
    void dirichletAtPos(PrimaryVariables &values,
                        const GlobalPosition &globalPos) const
    {
        initial_(values, globalPos);
    }

    //! \copydoc ImplicitProblem::neumann()
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

    //! \copydoc ImplicitProblem::source()
    void source(PrimaryVariables &values,
                const Element &element,
                const FVElementGeometry &fvGeometry,
                int scvIdx) const
    {
        // ATTENTION: The source term of the mass balance has to be chosen as
        // div (q_momentum) in the problem file
        values = Scalar(0.0);
    }

    //! \copydoc ImplicitProblem::initialAtPos()
    void initialAtPos(PrimaryVariables &values,
                      const GlobalPosition &globalPos) const
    {
        initial_(values, globalPos);
    }
   // \}

private:
    // internal method for the initial condition (reused for the
    // dirichlet conditions!)
    void initial_(PrimaryVariables &values,
                  const GlobalPosition &globalPos) const
    {
        values = 0.0;
        values[velocityYIdx] = velocity_
                               * (globalPos[0] - this->bBoxMin()[0])*(this->bBoxMax()[0] - globalPos[0])
                               / (0.25*(this->bBoxMax()[0] - this->bBoxMin()[0])*(this->bBoxMax()[0] - this->bBoxMin()[0]));
        values[pressureIdx] = 1e5 + 1.189*this->gravity()[1]*globalPos[1];
        values[massOrMoleFracIdx] = 1e-4;
        values[temperatureIdx] = 283.15;
        if(globalPos[0] < 0.75 && globalPos[0] > 0.25 &&
                globalPos[1] < 0.75 && globalPos[1] > 0.25)
        {
            values[massOrMoleFracIdx] = 0.9e-4;
            values[temperatureIdx] = 284.15;
        }

        if(useMoles)
        {
            Scalar M1 = FluidSystem::molarMass(transportCompIdx);
            Scalar M2 = FluidSystem::molarMass(phaseCompIdx);
            Scalar X1 = values[massOrMoleFracIdx];
            Scalar X2 = 1.0 - X1;

            values[massOrMoleFracIdx] = (X1/M1)/(X1/M1 + X2/M2);
        }
    }
    // \}

private:
    bool onLeftBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] < this->bBoxMin()[0] + eps_; }

    bool onRightBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] > this->bBoxMax()[0] - eps_; }

    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] < this->bBoxMin()[1] + eps_; }

    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] > this->bBoxMax()[1] - eps_; }

    bool onBoundary_(const GlobalPosition &globalPos) const
    {
        return (onLeftBoundary_(globalPos) || onRightBoundary_(globalPos)
                || onLowerBoundary_(globalPos) || onUpperBoundary_(globalPos));
    }

    Scalar eps_;
    Scalar velocity_;
};
} //end namespace

#endif
