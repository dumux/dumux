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
 * @brief  Definition of a simple two-component Stokes problem for the n-component model
 */
#ifndef DUMUX_STOKES2CTESTPROBLEM_HH
#define DUMUX_STOKES2CTESTPROBLEM_HH

#include <dune/grid/io/file/dgfparser/dgfyasp.hh>

#include <dumux/material/fluidsystems/h2oairfluidsystem.hh>

#include <dumux/freeflow/stokesnc/stokesncmodel.hh>

namespace Dumux
{

template <class TypeTag>
class Stokes2cTestProblem;

//////////
// Specify the properties for the stokes2c problem
//////////
namespace Properties
{
NEW_TYPE_TAG(Stokes2cTestProblem, INHERITS_FROM(BoxStokesnc));

// Set the grid type
SET_TYPE_PROP(Stokes2cTestProblem, Grid, Dune::YaspGrid<2>);

// Set the problem property
SET_TYPE_PROP(Stokes2cTestProblem, Problem, Dumux::Stokes2cTestProblem<TypeTag>);

//! Select the fluid system
SET_PROP(BoxStokesnc, FluidSystem)
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dumux::FluidSystems::H2OAir<Scalar> type;
};

//! Scalar is set to type long double for higher accuracy
//SET_TYPE_PROP(BoxStokes, Scalar, long double);
}

/*!
 * \ingroup BoxStokesncModel
 * \ingroup ImplicitTestProblems
 * \brief Stokes transport problem with dryer air flowing
 *        from the top to the bottom.
 *
 * The domain is sized 1m times 1m. Dry air enters the domain from the top boundary
 * and is transported downwards. The boundary conditions for the momentum balance
 * and the component transport equation are all set to Dirichlet, except on the lower
 * boundary, where outflow conditions are set. The mass balance receives
 * outflow bcs everywhere, which are replaced in the localresidual by the sum
 * of the momentum balance equations in case of Dirichlet bcs for the momentum balance.
 * In the middle of the lower boundary one vertex receives Dirichlet bcs, to set the pressure level.
 *
 * This problem uses the \ref StokesncModel.
 * To run the simulation execute the following line in a shell:
 * <tt>./test_stokes2c -parameterFile ./test_stokes2c.input</tt>
 */
template <class TypeTag>
class Stokes2cTestProblem : public StokesProblem<TypeTag>
{
    typedef StokesProblem<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    // Number of equations and grid dimension
    enum { dim = GridView::dimension };
    enum { // equation indices
        massBalanceIdx = Indices::massBalanceIdx,
        momentumXIdx = Indices::momentumXIdx, //!< Index of the x-component of the momentum balance
        momentumYIdx = Indices::momentumYIdx, //!< Index of the y-component of the momentum balance
        transportEqIdx = Indices::transportEqIdx  //!< Index of the transport equation
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
    Stokes2cTestProblem(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView)
    {
        eps_ = 1e-6;

        // initialize the tables of the fluid system
        FluidSystem::init();
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

    /*!
     * \brief Returns the temperature within the domain.
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperatureAtPos(const GlobalPosition &globalPos) const
    {
        return 273.15 + 10; // -> 10
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

        if (onLowerBoundary_(globalPos)
                && !onLeftBoundary_(globalPos) && !onRightBoundary_(globalPos))
            values.setAllOutflow();

        // the mass balance has to be of type outflow
        values.setOutflow(massBalanceIdx);

        // set pressure at one point
        const Scalar middle = (this->bBoxMax()[0] - this->bBoxMin()[0])/2;
        if (onLowerBoundary_(globalPos) &&
                globalPos[0] > middle - eps_ && globalPos[0] < middle + eps_)
            values.setDirichlet(massBalanceIdx);
    }

    //! \copydoc ImplicitProblem::dirichletAtPos()
    void dirichletAtPos(PrimaryVariables &values,
                        const GlobalPosition &globalPos) const
    {
        initial_(values, globalPos);

        if (onUpperBoundary_(globalPos))
        {
            values[massOrMoleFracIdx] = 0.005;
        }
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

        values[pressureIdx] = 1e5;
        values[velocityXIdx] = 0.0;

        //parabolic profile
        values[velocityYIdx] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Problem, Velocity)
                               * (globalPos[0] - this->bBoxMin()[0])*(this->bBoxMax()[0] - globalPos[0])
                               / (0.25*(this->bBoxMax()[0] - this->bBoxMin()[0])*(this->bBoxMax()[0] - this->bBoxMin()[0]));

        if (onUpperBoundary_(globalPos))
            values[massOrMoleFracIdx] = 0.005;
        else
            values[massOrMoleFracIdx] = 0.007;
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
};
} //end namespace

#endif
