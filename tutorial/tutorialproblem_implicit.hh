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
 * \brief Tutorial problem for a fully coupled twophase box model.
 */
#ifndef DUMUX_TUTORIAL_PROBLEM_IMPLICIT_HH    // guardian macro /*@\label{tutorial-implicit:guardian1}@*/
#define DUMUX_TUTORIAL_PROBLEM_IMPLICIT_HH    // guardian macro /*@\label{tutorial-implicit:guardian2}@*/

// The numerical model
#include <dumux/porousmediumflow/2p/implicit/model.hh>

// The base porous media box problem
#include <dumux/porousmediumflow/implicit/problem.hh>

// The DUNE grid used
#if  HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#elif HAVE_UG
#include <dune/grid/uggrid.hh>
#else
#include <dune/grid/yaspgrid.hh>
#endif // HAVE_DUNE_ALUGRID, HAVE_UG

// Spatially dependent parameters
#include "tutorialspatialparams_implicit.hh"

// The components that are used
#include <dumux/material/components/h2o.hh>
#include <dumux/material/components/lnapl.hh>
#include <dumux/io/cubegridcreator.hh>

namespace Dumux{
// Forward declaration of the problem class
template <class TypeTag>
class TutorialProblemImplicit;

namespace Properties {
// Create a new type tag for the problem
NEW_TYPE_TAG(TutorialProblemImplicit, INHERITS_FROM(BoxTwoP, TutorialSpatialParamsImplicit)); /*@\label{tutorial-implicit:create-type-tag}@*/

// Set the "Problem" property
SET_PROP(TutorialProblemImplicit, Problem) /*@\label{tutorial-implicit:set-problem}@*/
{ typedef Dumux::TutorialProblemImplicit<TypeTag> type;};

// Set grid and the grid creator to be used
#if HAVE_DUNE_ALUGRID /*@\label{tutorial-implicit:set-grid}@*/
SET_TYPE_PROP(TutorialProblemImplicit, Grid, Dune::ALUGrid</*dim=*/2, 2, Dune::cube, Dune::nonconforming>); /*@\label{tutorial-implicit:set-grid-ALU}@*/
#elif HAVE_UG
SET_TYPE_PROP(TutorialProblemImplicit, Grid, Dune::UGGrid<2>);
#else
SET_TYPE_PROP(TutorialProblemImplicit, Grid, Dune::YaspGrid<2>);
#endif // HAVE_DUNE_ALUGRID
SET_TYPE_PROP(TutorialProblemImplicit, GridCreator, Dumux::CubeGridCreator<TypeTag>); /*@\label{tutorial-implicit:set-gridcreator}@*/

// Set the wetting phase
SET_PROP(TutorialProblemImplicit, WettingPhase) /*@\label{tutorial-implicit:2p-system-start}@*/
{
private: typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public: typedef Dumux::FluidSystems::LiquidPhase<Scalar, Dumux::H2O<Scalar> > type; /*@\label{tutorial-implicit:wettingPhase}@*/
};

// Set the non-wetting phase
SET_PROP(TutorialProblemImplicit, NonwettingPhase)
{
private: typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public: typedef Dumux::FluidSystems::LiquidPhase<Scalar, Dumux::LNAPL<Scalar> > type; /*@\label{tutorial-implicit:nonwettingPhase}@*/
}; /*@\label{tutorial-implicit:2p-system-end}@*/

SET_TYPE_PROP(TutorialProblemImplicit, FluidSystem, Dumux::TwoPImmiscibleFluidSystem<TypeTag>);/*@\label{tutorial-implicit:set-fluidsystem}@*/
// Disable gravity
SET_BOOL_PROP(TutorialProblemImplicit, ProblemEnableGravity, false); /*@\label{tutorial-implicit:gravity}@*/
}

/*!
 * \ingroup TwoPBoxModel
 *
 * \brief  Tutorial problem for a fully coupled twophase box model.
 */
template <class TypeTag>
class TutorialProblemImplicit : public ImplicitPorousMediaProblem<TypeTag> /*@\label{tutorial-implicit:def-problem}@*/
{
    typedef ImplicitPorousMediaProblem<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    // Grid dimension
    enum { dim = GridView::dimension,
           dimWorld = GridView::dimensionworld
    };

    // Types from DUNE-Grid
    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::template Codim<dim>::Entity Vertex;
    typedef typename GridView::Intersection Intersection;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

    // Dumux specific types
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;

public:
    TutorialProblemImplicit(TimeManager &timeManager,
                           const GridView &gridView)
        : ParentType(timeManager, gridView)
        , eps_(3e-6)
    {
#if !(HAVE_DUNE_ALUGRID || HAVE_UG)
      std::cout << "If you want to use simplices instead of cubes, install and use dune-ALUGrid or UGGrid." << std::endl;
#endif // !(HAVE_DUNE_ALUGRID || HAVE_UG)
    }

    //! Specifies the problem name. This is used as a prefix for files
    //! generated by the simulation.
    const char *name() const
    { return "tutorial_implicit"; }

    //! Returns true if a restart file should be written.
    bool shouldWriteRestartFile() const /*@\label{tutorial-implicit:restart}@*/
    { return false; }

    //! Returns true if the current solution should be written to disk
    //! as a VTK file
    bool shouldWriteOutput() const /*@\label{tutorial-implicit:output}@*/
    {
        return
            (this->timeManager().timeStepIndex() > 0
             && (this->timeManager().timeStepIndex() % 1 == 0));
    }

    //! Returns the temperature within a finite volume. We use constant
    //! 10 degrees Celsius.
    Scalar temperature() const
    { return 283.15; }

    //! Specifies which kind of boundary condition should be used for
    //! which equation for a finite volume on the boundary.
    void boundaryTypes(BoundaryTypes &bcTypes, const Vertex &vertex) const
    {
        const GlobalPosition &globalPos = vertex.geometry().center();
        if (globalPos[0] < eps_) // Dirichlet conditions on left boundary
           bcTypes.setAllDirichlet();
        else // neuman for the remaining boundaries
           bcTypes.setAllNeumann();

    }

    //! Evaluates the Dirichlet boundary conditions for a finite volume
    //! on the grid boundary. Here, the 'values' parameter stores
    //! primary variables.
    void dirichlet(PrimaryVariables &values, const Vertex &vertex) const
    {
        values[Indices::pwIdx] = 200.0e3; // 200 kPa = 2 bar
        values[Indices::snIdx] = 0.0; // 0 % oil saturation on left boundary
    }

    //! Evaluates the boundary conditions for a Neumann boundary
    //! segment. Here, the 'values' parameter stores the mass flux in
    //! [kg/(m^2 * s)] in normal direction of each phase. Negative
    //! values mean influx.
    void neumann(PrimaryVariables &values,
                 const Element &element,
                 const FVElementGeometry &fvGeometry,
                 const Intersection &intersection,
                 int scvIdx,
                 int boundaryFaceIdx) const
    {
        const GlobalPosition &globalPos =
            fvGeometry.boundaryFace[boundaryFaceIdx].ipGlobal;
        Scalar right = this->bBoxMax()[0];
        // extraction of oil on the right boundary for approx. 1.e6 seconds
        if (globalPos[0] > right - eps_) {
            // oil outflux of 30 g/(m * s) on the right boundary.
            values[Indices::contiWEqIdx] = 0;
            values[Indices::contiNEqIdx] = 3e-2;
        } else {
            // no-flow on the remaining Neumann-boundaries.
            values[Indices::contiWEqIdx] = 0;
            values[Indices::contiNEqIdx] = 0;
        }
    }

    //! Evaluates the initial value for a control volume. For this
    //! method, the 'values' parameter stores primary variables.
    void initial(PrimaryVariables &values,
                 const Element &element,
                 const FVElementGeometry &fvGeometry,
                 int scvIdx) const
    {
        values[Indices::pwIdx] = 200.0e3; // 200 kPa = 2 bar
        values[Indices::snIdx] = 1.0;
    }

    //! Evaluates the source term for all phases within a given
    //! sub-control-volume. In this case, the 'values' parameter
    //! stores the rate mass generated or annihilated per volume unit
    //! in [kg / (m^3 * s)]. Positive values mean that mass is created.
    void source(PrimaryVariables &values,
                const Element &element,
                const FVElementGeometry &fvGeometry,
                int scvIdx) const
    {
        values[Indices::contiWEqIdx] = 0.0;
        values[Indices::contiNEqIdx]= 0.0;
    }

private:
    // small epsilon value
    Scalar eps_;
};
}

#endif
