// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2009 by Katherina Baber
 *   Copyright (C) 2009 by Karin Erbertseder                                 *
 *   Copyright (C) 2009 by Andreas Lauser                                    *
 *   Copyright (C) 2008 by Bernd Flemisch                                    *
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
 * \file
 * \brief Definition of a problem, for the 1p2c box problem:
 * Component transport of nitrogen dissolved in the water phase.
 */
#ifndef DUMUX_1P2C_OUTFLOW_PROBLEM_HH
#define DUMUX_1P2C_OUTFLOW_PROBLEM_HH

#ifdef HAVE_UG
#include <dune/grid/io/file/dgfparser/dgfug.hh>
#endif
#include <dune/grid/io/file/dgfparser/dgfs.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>

#include <dumux/boxmodels/1p2c/1p2cmodel.hh>

#include <dumux/material/fluidsystems/h2on2liquidphasefluidsystem.hh>
#include "1p2coutflowspatialparameters.hh"

namespace Dumux
{

template <class TypeTag>
class OnePTwoCOutflowProblem;

namespace Properties
{
NEW_TYPE_TAG(OnePTwoCOutflowProblem, INHERITS_FROM(BoxOnePTwoC));

// Set the grid type
SET_PROP(OnePTwoCOutflowProblem, Grid)
{
#if HAVE_UG
    typedef Dune::UGGrid<2> type;
#else
    typedef Dune::SGrid<2, 2> type;
    //typedef Dune::YaspGrid<2> type;
#endif
};

// Set the problem property
SET_PROP(OnePTwoCOutflowProblem, Problem)
{
    typedef Dumux::OnePTwoCOutflowProblem<TTAG(OnePTwoCOutflowProblem)> type;
};

// Set fluid configuration
SET_PROP(OnePTwoCOutflowProblem, FluidSystem)
{ private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
    typedef Dumux::FluidSystems::H2ON2LiquidPhase<Scalar, false> type;
};

// Set the spatial parameters
SET_TYPE_PROP(OnePTwoCOutflowProblem,
              SpatialParameters,
              Dumux::OnePTwoCOutflowSpatialParameters<TypeTag>);

//Define whether mole(true) or mass (false) fractions are used
SET_BOOL_PROP(OnePTwoCOutflowProblem, UseMoles, false);

// Disable gravity
SET_BOOL_PROP(OnePTwoCOutflowProblem, EnableGravity, false);
}


/*!
 * \ingroup OnePTwoCBoxModel
 * \ingroup BoxTestProblems
 *
 * \brief Definition of a problem, for the 1p2c box problem:
 * Nitrogen is dissolved in the water phase and
 * is transported with the water flow from the left side to the right.
 *
 * The model domain is 1m times 1m with a discretization length of 0.05m
 * and homogeneous soil properties (\f$ \mathrm{K=10e-10, \Phi=0.4}\f$).
 * Initially the domain is filled with pure water.
 *
 * At the left side, a Dirichlet condition defines a nitrogen mole fraction
 * of 0.3 mol/mol.
 * The water phase flows from the left side to the right due to the applied pressure
 * gradient of 1e5Pa/m. The nitrogen is transported with the water flow
 * and leaves the domain at the right boundary
 * where an outflow boundary condition is applied.
 * This problem uses the \ref OnePTwoCBoxModel.
 *
 * To run the simulation execute the following line in shell:
 * <tt>./test_1p2c -parameterFile ./test_1p2c.input</tt>
 */
template <class TypeTag>
class OnePTwoCOutflowProblem : public OnePTwoCBoxProblem<TypeTag>
{
    typedef OnePTwoCBoxProblem<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

    // copy some indices for convenience
    typedef typename GET_PROP_TYPE(TypeTag, OnePTwoCIndices) Indices;
    enum {
        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld,

        // indices of the primary variables
        pressureIdx = Indices::pressureIdx,
        x1Idx = Indices::x1Idx,

        // indices of the equations
        contiEqIdx = Indices::contiEqIdx,
        transEqIdx = Indices::transEqIdx
    };


    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::Intersection Intersection;

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:
    OnePTwoCOutflowProblem(TimeManager &timeManager, const GridView &gridView)
        : ParentType(timeManager, gridView)
        , eps_(1e-6)
    {
        //initialize fluid system
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
    { return "outflow"; }

    /*!
     * \brief Returns the temperature within the domain.
     *
     * This problem assumes a temperature of 20 degrees Celsius.
     */
    Scalar temperature() const
    { return 273.15 + 20; }; // in [K]

    // \}

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     */
    void boundaryTypes(BoundaryTypes &values, const Vertex &vertex) const
    {
        const GlobalPosition globalPos = vertex.geometry().center();

        if(globalPos[0] < eps_ || globalPos[0] > this->bboxMax()[0] - eps_)
            values.setAllDirichlet();
        else
            values.setAllNeumann();

        //outflow condition for the transport equation at right boundary
        if(globalPos[0] > this->bboxMax()[0] - eps_)
            values.setOutflow(transEqIdx);
    }

    /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        boundary segment.
     *
     * For this method, the \a values parameter stores primary variables.
     */
    void dirichlet(PrimaryVariables &values, const Vertex &vertex) const
    {
        const GlobalPosition globalPos = vertex.geometry().center();

        initial_(values, globalPos);
        //condition for the trail molefraction at left boundary
        if(globalPos[0] < eps_)
            values[x1Idx] = 0.3;
    }

    /*!
     * \brief Evaluate the boundary conditions for a neumann
     *        boundary segment.
     *
     * For this method, the \a values parameter stores the mass flux
     * in normal direction of each component. Negative values mean
     * influx.
     */
    void neumann(PrimaryVariables &values,
                 const Element &element,
                 const FVElementGeometry &fvElemGeom,
                 const Intersection &is,
                 int scvIdx,
                 int boundaryFaceIdx) const
    {
        //const GlobalPosition &globalPos = element.geometry().corner(scvIdx);
        values = 0;
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
     * of a component is generated or annihilate per volume
     * unit. Positive values mean that mass is created, negative ones
     * mean that it vanishes.
     */
    void sourceAtPos(PrimaryVariables &values,
                     const GlobalPosition &globalPos) const
    {
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
        const GlobalPosition &globalPos
            = element.geometry().corner(scvIdx);

        initial_(values, globalPos);
    }

    // \}

private:
    // the internal method for the initial condition
    void initial_(PrimaryVariables &values,
                  const GlobalPosition &globalPos) const
    {
        values[pressureIdx] = 2e5 - 1e5*globalPos[0];//0.0; //initial condition for the pressure
        values[x1Idx] = 0.0; //initial condition for the trail molefraction
    }

    const Scalar eps_;
};

} //end namespace
#endif
