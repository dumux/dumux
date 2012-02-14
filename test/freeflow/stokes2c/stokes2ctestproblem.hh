// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009-2010 by Klaus Mosthaf                                *
 *   Copyright (C) 2007-2008 by Bernd Flemisch                               *
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
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
/**
 * @file
 * @brief  Definition of a simple Stokes problem
 * @author Klaus Mosthaf, Andreas Lauser, Bernd Flemisch
 */
#ifndef DUMUX_STOKES2CTESTPROBLEM_HH
#define DUMUX_STOKES2CTESTPROBLEM_HH

#if HAVE_UG
#include <dune/grid/io/file/dgfparser/dgfug.hh>
#endif
#include <dune/grid/io/file/dgfparser/dgfs.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>

#include <dumux/material/fluidsystems/h2oairfluidsystem.hh>

#include <dumux/freeflow/stokes2c/stokes2cmodel.hh>

namespace Dumux
{

template <class TypeTag>
class Stokes2cTestProblem;

//////////
// Specify the properties for the stokes2c problem
//////////
namespace Properties
{
NEW_TYPE_TAG(Stokes2cTestProblem, INHERITS_FROM(BoxStokes2c));

// Set the grid type
SET_TYPE_PROP(Stokes2cTestProblem, Grid, Dune::SGrid<2,2>);

// Set the problem property
SET_TYPE_PROP(Stokes2cTestProblem, Problem, Dumux::Stokes2cTestProblem<TypeTag>);

//! Select the fluid system
SET_PROP(BoxStokes2c, FluidSystem)
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dumux::FluidSystems::H2OAir<Scalar> type;
};

//! Scalar is set to type long double for higher accuracy
//SET_TYPE_PROP(BoxStokes, Scalar, long double);

//! a stabilization factor. Set to zero for no stabilization
SET_SCALAR_PROP(BoxStokes2c, StabilizationAlpha, -1.0);

//! stabilization at the boundaries
SET_SCALAR_PROP(BoxStokes2c, StabilizationBeta, 0.0);

// Enable gravity
SET_BOOL_PROP(Stokes2cTestProblem, EnableGravity, false);
}

/*!
 * \ingroup BoxStokes2cModel
 * \ingroup BoxTestProblems
 * \brief Stokes transport problem with air flowing
 *        from the left to the right.
 *
 * The domain is sized 1m times 1m. The boundary conditions for the momentum balances
 * are all set to Dirichlet. The mass balance receives
 * outflow bcs, which are replaced in the localresidual by the sum
 * of the two momentum balances. In the middle of the right boundary,
 * one vertex receives Dirichlet bcs, to set the pressure level.
 *
 * This problem uses the \ref BoxStokes2cModel.
 *
 * This problem is non-stationary and can be simulated until \f$t_{\text{end}} =
 * 1e5\;s\f$ is reached. A good choice for the initial time step size
 * is \f$t_{\text{inital}} = 1\;s\f$.
 * To run the simulation execute the following line in shell:
 * <tt>./test_stokes2c -parameterFile ./test_stokes2c.input</tt>
 */
template <class TypeTag>
class Stokes2cTestProblem : public StokesProblem<TypeTag>
{
    typedef StokesProblem<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, Stokes2cIndices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

    // Number of equations and grid dimension
    enum { dim = GridView::dimension };
    enum {
        // copy some indices for convenience
        massBalanceIdx = Indices::massBalanceIdx,
        momentumXIdx = Indices::momentumXIdx, //!< Index of the x-component of the momentum balance
        momentumYIdx = Indices::momentumYIdx, //!< Index of the y-component of the momentum balance
        transportIdx = Indices::transportIdx  //!< Index of the transport equation (massfraction)
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
    const char *name() const
    { return "stokes2c"; }

    /*!
     * \brief Returns the temperature within the domain.
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar boxTemperature(const Element &element,
                       const FVElementGeometry &fvElemGeom,
                       int scvIdx) const
    {
        return 273.15 + 10; // -> 10
    };

    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param values The boundary types for the conservation equations
     * \param vertex The vertex on the boundary for which the
     *               conditions needs to be specified
     */
    void boundaryTypes(BoundaryTypes &values, const Vertex &vertex) const
    {
        const GlobalPosition globalPos = vertex.geometry().center();

        values.setAllDirichlet();

        if (onLowerBoundary_(globalPos)
                && !onLeftBoundary_(globalPos) && !onRightBoundary_(globalPos))
            values.setAllOutflow();

        // the mass balance has to be of type outflow
        values.setOutflow(massBalanceIdx);

        // set pressure at one point
        const Scalar middle = (this->bboxMax()[0] - this->bboxMin()[0])/2;
        if (onLowerBoundary_(globalPos) &&
                globalPos[0] > middle - eps_ && globalPos[1] < middle + eps_)
            values.setDirichlet(massBalanceIdx);
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume.
     *
     * \param values The dirichlet values for the primary variables
     * \param vertex The vertex representing the "half volume on the boundary"
     *
     * For this method, the \a values parameter stores primary variables.
     */
    void dirichlet(PrimaryVariables &values, const Vertex &vertex) const
    {
        const GlobalPosition globalPos = vertex.geometry().center();
        initial_(values, globalPos);

        if (onUpperBoundary_(globalPos))
        {
            values[transportIdx] = 0.005;
        }
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * For this method, the \a values parameter stores the mass flux
     * in normal direction of each phase. Negative values mean influx.
     */
    void neumann(PrimaryVariables &values,
                 const Element &element,
                 const FVElementGeometry &fvElemGeom,
                 const Intersection &is,
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

    /*!
     * \brief Evaluate the source term for all phases within a given
     *        sub-control-volume.
     *
     * For this method, the \a values parameter stores the rate mass
     * generated or annihilate per volume unit. Positive values mean
     * that mass is created, negative ones mean that it vanishes.
     */
    void source(PrimaryVariables &values,
                const Element &element,
                const FVElementGeometry &,
                int subControlVolumeIdx) const
    {
        // ATTENTION: The source term of the mass balance has to be chosen as
        // div (q_momentum) in the problem file
        values = Scalar(0.0);
    }

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * For this method, the \a values parameter stores primary
     * variables.
     */
    void initial(PrimaryVariables &values,
                 const Element &element,
                 const FVElementGeometry &fvElemGeom,
                 int scvIdx) const
    {
        const GlobalPosition &globalPos = element.geometry().corner(scvIdx);
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

       values[massBalanceIdx] = 1e5;
       values[momentumXIdx] = 0.0;

       //parabolic profile
       const Scalar v1 = 1.0;
       values[momentumYIdx] = -v1*(globalPos[0] - this->bboxMin()[0])*(this->bboxMax()[0] - globalPos[0])
                                    / (0.25*(this->bboxMax()[0] - this->bboxMin()[0])*(this->bboxMax()[0] - this->bboxMin()[0]));

       if (onUpperBoundary_(globalPos))
           values[transportIdx] = 0.005;
       else
           values[transportIdx] = 0.007;
   }
    bool onLeftBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] < this->bboxMin()[0] + eps_; }

    bool onRightBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[0] > this->bboxMax()[0] - eps_; }

    bool onLowerBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] < this->bboxMin()[1] + eps_; }

    bool onUpperBoundary_(const GlobalPosition &globalPos) const
    { return globalPos[1] > this->bboxMax()[1] - eps_; }

    Scalar eps_;
};
} //end namespace

#endif
