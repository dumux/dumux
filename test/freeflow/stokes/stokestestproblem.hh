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
#ifndef DUMUX_STOKESTESTPROBLEM_HH
#define DUMUX_STOKESTESTPROBLEM_HH

#include <dumux/material/fluidsystems/h2on2.hh>
#include <dumux/material/fluidsystems/gasphase.hh>

#include <dumux/freeflow/stokes/model.hh>

namespace Dumux
{

template <class TypeTag>
class StokesTestProblem;

//////////
// Specify the properties for the stokes problem
//////////
namespace Properties
{
NEW_TYPE_TAG(StokesTestProblem, INHERITS_FROM(BoxStokes));

// Set the grid type
SET_TYPE_PROP(StokesTestProblem, Grid, Dune::YaspGrid<2>);

// Set the problem property
SET_TYPE_PROP(StokesTestProblem, Problem, StokesTestProblem<TypeTag>);

// Use nitrogen as gas phase
SET_TYPE_PROP(StokesTestProblem, Fluid,
              FluidSystems::GasPhase<typename GET_PROP_TYPE(TypeTag, Scalar),
                                            N2<typename GET_PROP_TYPE(TypeTag, Scalar)> >);
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
class StokesTestProblem : public StokesProblem<TypeTag>
{
    typedef StokesProblem<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
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

    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::ctype CoordScalar;
    typedef typename GridView::Intersection Intersection;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, Fluid) Fluid;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;

public:
    StokesTestProblem(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView)
    {}

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
    { return "stokes"; }

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

        if(onRightBoundary_(globalPos) &&
                globalPos[1] < this->bBoxMax()[1]-eps_ && globalPos[1] > this->bBoxMin()[1]+eps_)
            values.setAllOutflow();

        // set pressure at one point
        const Scalar middle = (this->bBoxMax()[1] - this->bBoxMin()[1])/2;
        if (onRightBoundary_(globalPos) &&
                globalPos[1] > middle - eps_ && globalPos[1] < middle + eps_)
            values.setDirichlet(pressureIdx);
    }

    //! \copydoc ImplicitProblem::dirichletAtPos()
    void dirichletAtPos(PrimaryVariables &values,
                        const GlobalPosition &globalPos) const
    {
        const Scalar time = this->timeManager().time() + this->timeManager().timeStepSize();
        const Scalar velocityVariation = 0.2;

        initial_(values, globalPos);

        // sinusoidal variation of the maximum velocity in time
        using std::sin;
        const Scalar v0 = 1.0 + sin(2*M_PI*time/3000) * velocityVariation;

        // parabolic velocity profile
        values[velocityXIdx] =  v0*(globalPos[1] - this->bBoxMin()[1])*(this->bBoxMax()[1] - globalPos[1])
                               / (0.25*(this->bBoxMax()[1] - this->bBoxMin()[1])*(this->bBoxMax()[1] - this->bBoxMin()[1]));
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     *
     * \param values The neumann values for the conservation equations
     * \param element The finite element
     * \param fvGeometry The finite-volume geometry in the box scheme
     * \param intersection The intersection between element and boundary
     * \param scvIdx The local vertex index
     * \param boundaryFaceIdx The index of the boundary face
     *
     * For this method, the \a values parameter stores the mass flux
     * in normal direction of each phase. Negative values mean influx.
     *
     * A NEUMANN condition for the Stokes equation corresponds to:
     * \f[ -\mu \nabla {\bf v} \cdot {\bf n} + p \cdot {\bf n} = q_N \f]
     *
     *
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
    }

    bool onLeftBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] < this->bBoxMin()[0] + eps_; }

    bool onRightBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] > this->bBoxMax()[0] - eps_; }

    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] < this->bBoxMin()[1] + eps_; }

    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] > this->bBoxMax()[1] - eps_; }

    static constexpr Scalar eps_ = 1e-6;
};

} //end namespace

#endif
