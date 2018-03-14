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
 * \ingroup RANSTests
 * \brief Pipe flow test for the staggered grid RANS model
 *
 * This test simulates is based on pipe flow experiments by
 * John Laufers experiments in 1954 \cite Laufer1954a.
 */
#ifndef DUMUX_PIPE_LAUFER_PROBLEM_HH
#define DUMUX_PIPE_LAUFER_PROBLEM_HH

#include <dumux/material/fluidsystems/gasphase.hh>
#include <dumux/material/components/air.hh>

#include <dumux/freeflow/rans/zeroeq/model.hh>
#include <dumux/freeflow/rans/problem.hh>
#include <dumux/discretization/staggered/freeflow/properties.hh>

namespace Dumux
{
template <class TypeTag>
class PipeLauferProblem;

namespace Properties
{
NEW_TYPE_TAG(PipeLauferProblem, INHERITS_FROM(StaggeredFreeFlowModel, ZeroEq));

// the fluid system
SET_PROP(PipeLauferProblem, FluidSystem)
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using type = FluidSystems::GasPhase<Scalar, Air<Scalar> >;
};

// Set the grid type
SET_TYPE_PROP(PipeLauferProblem, Grid,
              Dune::YaspGrid<2, Dune::TensorProductCoordinates<typename GET_PROP_TYPE(TypeTag, Scalar), 2> >);

// Set the problem property
SET_TYPE_PROP(PipeLauferProblem, Problem, Dumux::PipeLauferProblem<TypeTag> );

SET_BOOL_PROP(PipeLauferProblem, EnableFVGridGeometryCache, true);

SET_BOOL_PROP(PipeLauferProblem, EnableGridFluxVariablesCache, true);
SET_BOOL_PROP(PipeLauferProblem, EnableGridVolumeVariablesCache, true);
}

/*!
 * \ingroup NavierStokesTests
 * \brief  Test problem for the one-phase (Navier-) Stokes problem in a channel.
 *
 * This test simulates is based on pipe flow experiments by
 * John Laufers experiments in 1954 \cite Laufer1954a.
 */
template <class TypeTag>
class PipeLauferProblem : public RANSProblem<TypeTag>
{
    using ParentType = RANSProblem<TypeTag>;

    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);

    // copy some indices for convenience
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    enum { dimWorld = GridView::dimensionworld };
    enum {
        massBalanceIdx = Indices::massBalanceIdx,
        momentumBalanceIdx = Indices::momentumBalanceIdx,
        pressureIdx = Indices::pressureIdx,
        velocityXIdx = Indices::velocityXIdx,
        velocityYIdx = Indices::velocityYIdx
    };

    using BoundaryTypes = typename GET_PROP_TYPE(TypeTag, BoundaryTypes);

    using FVGridGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry);

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using SourceValues = typename GET_PROP_TYPE(TypeTag, NumEqVector);

    using TimeLoopPtr = std::shared_ptr<CheckPointTimeLoop<Scalar>>;

public:
    PipeLauferProblem(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry), eps_(1e-6)
    {
        inletVelocity_ = getParam<Scalar>("Problem.InletVelocity");
    }

   /*!
     * \name Problem parameters
     */
    // \{

    bool isOnWall(const GlobalPosition &globalPos) const
    {
        Scalar localEps_ = 1e-6; // cannot use the epsilon, because ParentType is initialized first
        return globalPos[1] < this->fvGridGeometry().bBoxMin()[1] + localEps_
              || globalPos[1] > this->fvGridGeometry().bBoxMax()[1] - localEps_;
    }

    bool shouldWriteRestartFile() const
    {
        return false;
    }

   /*!
     * \brief Return the temperature within the domain in [K].
     *
     * This problem assumes a temperature of 10 degrees Celsius.
     */
    Scalar temperature() const
    { return 273.15 + 10; } // 10C

   /*!
     * \brief Return the sources within the domain.
     *
     * \param globalPos The global position
     */
    SourceValues sourceAtPos(const GlobalPosition &globalPos) const
    {
        return SourceValues(0.0);
    }
    // \}
   /*!
     * \name Boundary conditions
     */
    // \{

   /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary control volume.
     *
     * \param globalPos The position of the center of the finite volume
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;

        // set Dirichlet values for the velocity everywhere
        values.setDirichlet(momentumBalanceIdx);

        // set a fixed pressure in one cell
        if (isOutlet(globalPos))
        {
            values.setDirichlet(massBalanceIdx);
            values.setOutflow(momentumBalanceIdx);
        }
        else
            values.setOutflow(massBalanceIdx);

        return values;
    }

   /*!
     * \brief Evaluate the boundary conditions for a dirichlet
     *        control volume.
     *
     * \param globalPos The center of the finite volume which ought to be set.
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        return initialAtPos(globalPos);
    }

    // \}

   /*!
     * \name Volume terms
     */
    // \{

   /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * \param globalPos The global position
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values(0.0);
        values[pressureIdx] = 1.0e+5;
        if(!isOnWall(globalPos))
        {
            values[velocityXIdx] = inletVelocity_;
        }

        return values;
    }

    // \}
    void setTimeLoop(TimeLoopPtr timeLoop)
    {
        timeLoop_ = timeLoop;
        if(inletVelocity_ > eps_)
            timeLoop_->setCheckPoint({200.0, 210.0});
    }

    Scalar time() const
    {
        return timeLoop_->time();
    }

private:

    bool isInlet(const GlobalPosition& globalPos) const
    {
        return globalPos[0] < eps_;
    }

    bool isOutlet(const GlobalPosition& globalPos) const
    {
        return globalPos[0] > this->fvGridGeometry().bBoxMax()[0] - eps_;
    }

    bool isWall(const GlobalPosition& globalPos) const
    {
        return globalPos[0] > eps_ || globalPos[0] < this->fvGridGeometry().bBoxMax()[0] - eps_;
    }

    Scalar eps_;
    Scalar inletVelocity_;
    TimeLoopPtr timeLoop_;
};
} //end namespace

#endif
