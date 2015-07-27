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
 * \brief Problem that uses a generalized Dirichlet boundary condition.
 */
#ifndef DUMUX_GENERALIZED_DIRICHLET_PROBLEM_HH
#define DUMUX_GENERALIZED_DIRICHLET_PROBLEM_HH

// The numerical model
#include <dumux/implicit/2p/2pmodel.hh>

// The base porous media box problem
#include <dumux/implicit/common/implicitporousmediaproblem.hh>

// The DUNE grid used
#include <dune/grid/yaspgrid.hh>

// Spatially dependent parameters
#include "generalizeddirichletspatialparams.hh"

// The components that are used
#include <dumux/material/components/h2o.hh>
#include <dumux/material/components/lnapl.hh>
#include <dumux/io/cubegridcreator.hh>
#include <dumux/linear/seqsolverbackend.hh>

namespace Dumux{
// Forward declaration of the problem class
template <class TypeTag>
class GeneralizedDirichletProblem;

namespace Properties {
// Create a new type tag for the problem
NEW_TYPE_TAG(GeneralizedDirichletProblem, INHERITS_FROM(BoxTwoP, GeneralizedDirichletSpatialParams));

// Set the "Problem" property
SET_PROP(GeneralizedDirichletProblem, Problem)
{ typedef Dumux::GeneralizedDirichletProblem<TypeTag> type;};

// Set grid and the grid creator to be used
SET_TYPE_PROP(GeneralizedDirichletProblem, Grid, Dune::YaspGrid<1>);
SET_TYPE_PROP(GeneralizedDirichletProblem, GridCreator, Dumux::CubeGridCreator<TypeTag>);

// Set the wetting phase
SET_PROP(GeneralizedDirichletProblem, WettingPhase)
{
private: typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public: typedef Dumux::LiquidPhase<Scalar, Dumux::H2O<Scalar> > type;
};

// Set the non-wetting phase
SET_PROP(GeneralizedDirichletProblem, NonwettingPhase)
{
private: typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public: typedef Dumux::LiquidPhase<Scalar, Dumux::LNAPL<Scalar> > type;
};

SET_INT_PROP(GeneralizedDirichletProblem, Formulation, TwoPFormulation::pnsw);

SET_BOOL_PROP(GeneralizedDirichletProblem, ProblemEnableGravity, false);
}

/*!
 * \ingroup TwoPBoxModel
 *
 * \brief  Problem that uses a generalized Dirichlet boundary condition.
 *
 * On the right boundary, a Dirichlet value for the nonwetting-phase pressure
 * should be set. At the same time, nonwetting phase should be allowed to leave
 * the domain by means of an outflow boundary condition. In order to achieve
 * this the general 'setDirichlet' method is used in the function
 * 'boundaryTypesAtPos' which allows to set the value of a primary variable by
 * replacing an arbitrary balance equation. In this case, the index of the
 * primary variable is Indices::pnIdx, while the index of the replaced equation
 * is Indices::contiWEqIdx. This allows to set the outflow condition for
 * equation index Indices::contiNEqIdx.
 */
template <class TypeTag>
class GeneralizedDirichletProblem : public ImplicitPorousMediaProblem<TypeTag>
{
    typedef ImplicitPorousMediaProblem<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;

    // Grid dimension
    enum { dim = GridView::dimension,
           dimWorld = GridView::dimensionworld
    };

    // Types from DUNE-Grid
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

    // Dumux specific types
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;
    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;
    typedef typename GET_PROP_TYPE(TypeTag, BoundaryTypes) BoundaryTypes;

public:
    GeneralizedDirichletProblem(TimeManager &timeManager,
                                const GridView &gridView)
    : ParentType(timeManager, gridView)
    , eps_(3e-6)
    {
        name_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Problem, Name);
    }

    //! Specifies the problem name. This is used as a prefix for files
    //! generated by the simulation.
    const char *name() const
    { return name_.c_str(); }

    //! Returns the temperature within a finite volume. We use constant
    //! 10 degrees Celsius.
    Scalar temperature() const
    { return 283.15; };

    bool shouldWriteRestartFile() const
    {
        return false;
    }

    //! Specifies which kind of boundary condition should be used for
    //! which equation for a finite volume on the boundary.
    void boundaryTypesAtPos(BoundaryTypes &bcTypes,
                            const GlobalPosition &globalPos) const
    {
        if (globalPos[0] > this->bBoxMax()[0] - eps_)
        {
           bcTypes.setDirichlet(Indices::pnIdx, Indices::contiWEqIdx);
           bcTypes.setOutflow(Indices::contiNEqIdx);
        }
        else // Neumann for the remaining boundaries
           bcTypes.setAllNeumann();

    }

    //! Evaluates the Dirichlet boundary conditions for a finite volume
    //! on the grid boundary. Here, the 'values' parameter stores
    //! primary variables.
    void dirichletAtPos(PrimaryVariables &values,
                        const GlobalPosition &globalPos) const
    {
        values = 0;
        if (globalPos[0] > this->bBoxMax()[0] - eps_)
            values[Indices::pnIdx] = 200e3; // 200 kPa = 2 bar oil pressure on right boundary
    }

    //! Evaluates the boundary conditions for a Neumann boundary
    //! segment. Here, the 'values' parameter stores the mass flux in
    //! [kg/(m^2 * s)] in normal direction of each phase. Negative
    //! values mean influx.
    void neumannAtPos(PrimaryVariables &values,
                      const GlobalPosition &globalPos) const
    {
        // water influx of 0.05 kg / (m.s) on the left boundary
        if (globalPos[0] < eps_) {
            values[Indices::contiWEqIdx] = -5e-2;
            values[Indices::contiNEqIdx] = 0;
        }
    }

    //! Evaluates the initial value for a control volume. For this
    //! method, the 'values' parameter stores primary variables.
    void initialAtPos(PrimaryVariables &values,
                      const GlobalPosition &globalPos) const
    {
        values[Indices::pnIdx] = 200e3; // 200 kPa = 2 bar
        values[Indices::swIdx] = 0;
    }

    //! Evaluates the source term for all phases within a given
    //! sub-control-volume. In this case, the 'values' parameter
    //! stores the rate mass generated or annihilated per volume unit
    //! in [kg / (m^3 * s)]. Positive values mean that mass is created.
    void sourceAtPos(PrimaryVariables &values,
                     const GlobalPosition &globalPos) const
    {
        values = 0;
    }

private:
    Scalar eps_;
    std::string name_;
};
}

#endif
