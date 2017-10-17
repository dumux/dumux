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
 * \brief The problem class for the coupling of an isothermal two-component Stokes
 *        and an isothermal two-phase two-component Darcy model.
 *
 * The problem class for the coupling of an isothermal two-component Stokes (stokes2c)
 * and an isothermal two-phase two-component Darcy model (2p2c).
 * The initial and boundary conditions of the submodels are specified in the two subproblems,
 * 2p2csubproblem.hh and stokes2csubproblem.hh, which are accessible via the coupled problem.
 */

#ifndef DUMUX_2CSTOKES2P2CPROBLEM_HH
#define DUMUX_2CSTOKES2P2CPROBLEM_HH

// TODO copied and adapted from stokesdarcy1p
#include "2p2csubproblem.hh"
#include "stokes2csubproblem.hh"
#include <dumux/material/fluidsystems/h2oair.hh>

#include <dumux/multidomain/problem.hh>
#include <dumux/multidomain/couplingmanager.hh>
#include <dumux/multidomain/map.hh>

#include <dumux/multidomain/staggered-ccfv/properties.hh>
#include <dumux/multidomain/staggered-ccfv/model.hh>
#include <dumux/multidomain/staggered-ccfv/assembler.hh>
#include <dumux/multidomain/staggered-ccfv/newtoncontroller.hh>
#include <dumux/multidomain/staggered-ccfv/subproblemlocaljacobian.hh>

#ifdef HAVE_PARDISO
#include <dumux/linear/pardisobackend.hh>
#endif // HAVE_PARDISO

namespace Dumux
{
template <class TypeTag>
class TwoCStokesTwoPTwoCTestProblem;

namespace Properties
{
NEW_TYPE_TAG(TwoCStokesTwoPTwoCTestProblem, INHERITS_FROM(MultiDomain, CouplingStokesStaggeredModel));

// Set the problem property
SET_TYPE_PROP(TwoCStokesTwoPTwoCTestProblem, Problem, Dumux::TwoCStokesTwoPTwoCTestProblem<TypeTag>);

// Set the coupling manager
SET_TYPE_PROP(TwoCStokesTwoPTwoCTestProblem, CouplingManager, Dumux::CouplingManagerStokesDarcy<TypeTag>);

//////////////////////////////////////////////////////////////////////////
// Set the two sub-problems of the global problem
SET_TYPE_PROP(TwoCStokesTwoPTwoCTestProblem, DarcyProblemTypeTag, TTAG(TwoPTwoCSubProblem));
SET_TYPE_PROP(TwoCStokesTwoPTwoCTestProblem, StokesProblemTypeTag, TTAG(Stokes2cSubProblem));
////////////////////////////////////////////////////////////////////////////

// publish this problem in the sub problems
SET_TYPE_PROP(TwoPTwoCSubProblem, GlobalProblemTypeTag, TTAG(TwoCStokesTwoPTwoCTestProblem));
SET_TYPE_PROP(Stokes2cSubProblem, GlobalProblemTypeTag, TTAG(TwoCStokesTwoPTwoCTestProblem));

// The subproblems inherit the parameter tree from this problem
SET_PROP(TwoPTwoCSubProblem, ParameterTree) : GET_PROP(TTAG(TwoCStokesTwoPTwoCTestProblem), ParameterTree) {};
SET_PROP(Stokes2cSubProblem, ParameterTree) : GET_PROP(TTAG(TwoCStokesTwoPTwoCTestProblem), ParameterTree) {};

NEW_PROP_TAG(DarcyToStokesMapValue); // TODO: make specialized map value class
SET_TYPE_PROP(TwoCStokesTwoPTwoCTestProblem, DarcyToStokesMapValue, Dumux::DarcyToStokesMap<TypeTag>);

// Set the fluid system to use simple relations (last argument)
SET_TYPE_PROP(TwoCStokesTwoPTwoCTestProblem, FluidSystem, H2OAirFluidSystem<TypeTag>);

// if you do not have PARDISO, the SuperLU solver is used:
#ifdef HAVE_PARDISO
SET_TYPE_PROP(TwoCStokesTwoPTwoCTestProblem, LinearSolver, PardisoBackend<TypeTag>);
#else
SET_TYPE_PROP(TwoCStokesTwoPTwoCTestProblem, LinearSolver, SuperLUBackend<TypeTag>);
#endif // HAVE_PARDISO
} // end namespace Properties

/*!
 * \ingroup ImplicitTestProblems
 * \ingroup TwoPTwoCStokesTwoCModel
 * \brief The problem class for the coupling of an isothermal two-component Stokes
 *        and an isothermal two-phase two-component Darcy model.
 *
 * The problem class for the coupling of an isothermal two-component Stokes (stokes2c)
 * and an isothermal two-phase two-component Darcy model (2p2c).
 * The initial and boundary conditions of the submodels are specified in the two subproblems,
 * 2p2csubproblem.hh and stokes2csubproblem.hh, which are accessible via the coupled problem.
 */
template <class TypeTag>
class TwoCStokesTwoPTwoCTestProblem : public MultiDomainProblem<TypeTag>
{
    using ParentType = MultiDomainProblem<TypeTag>;
    using TimeManager = typename GET_PROP_TYPE(TypeTag, TimeManager);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);

    // obtain the type tags of the sub problems
    using StokesProblemTypeTag = typename GET_PROP_TYPE(TypeTag, StokesProblemTypeTag);
    using DarcyProblemTypeTag = typename GET_PROP_TYPE(TypeTag, DarcyProblemTypeTag);

    // obtain types from the sub problem type tags
    using StokesProblem = typename GET_PROP_TYPE(StokesProblemTypeTag, Problem);
    using DarcyProblem = typename GET_PROP_TYPE(DarcyProblemTypeTag, Problem);

    using StokesGridView = typename GET_PROP_TYPE(StokesProblemTypeTag, GridView);
    using DarcyGridView = typename GET_PROP_TYPE(DarcyProblemTypeTag, GridView);

    // TODO not needed for stokesdarcy1p, necessary for 2cstokes2p2c?
//    using GlobalPosition = Dune::FieldVector<Scalar, dim>;
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);

public:
    /*!
     * \brief The problem for the coupling of Stokes and Darcy flow
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     */
    TwoCStokesTwoPTwoCTestProblem(TimeManager &timeManager, const StokesGridView &stokesGridView, const DarcyGridView &darcygridView)
    : ParentType(timeManager, stokesGridView, darcygridView)
    {
        // TODO from stokesdarcy1p
        verbose_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, bool, Problem, Verbose);

        // TODO from 2cstokes2p2c
        interfacePosY_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, DarcyGrid, InterfacePosY);
        noDarcyX_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, DarcyGrid, NoDarcyX);
        episodeLength_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, TimeManager, EpisodeLength);
        initializationTime_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, TimeManager, InitTime);

        // define output options
        freqRestart_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, Output, FreqRestart);
        freqOutput_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, Output, FreqOutput);

//        stokes2c_ = this->sdID1();
//        twoPtwoC_ = this->sdID2();

        // initialize the tables of the fluid system
        FluidSystem::init(/*tempMin=*/273.15, /*tempMax=*/323.15, /*numTemp=*/50,
                          /*pMin=*/5e4, /*pMax=*/1.5e5, /*numP=*/100);

        if (initializationTime_ > 0.0)
            this->timeManager().startNextEpisode(initializationTime_);
        else
            this->timeManager().startNextEpisode(episodeLength_);
    }

    /*!
     * \brief Called by the time manager after the time integration to
     *        do some post processing on the solution.
     */
    void postTimeStep()
    {
        // TODO!
//        // call the postTimeStep function of the subproblems
//        this->sdProblem1().postTimeStep();
//        this->sdProblem2().postTimeStep();
    }

    /*!
     * \brief Called when the end of an simulation episode is reached.
     *
     * Typically a new episode should be started in this method.
     */
    void episodeEnd()
    { this->timeManager().startNextEpisode(episodeLength_); }

    /*!
     * \brief Returns true if a restart file should be written to
     *        disk.
     *
     * The default behavior is to write one restart file every 5 time
     * steps. This file is intended to be overwritten by the
     * implementation.
     */
    bool shouldWriteRestartFile() const
    {
        return (this->timeManager().timeStepIndex() % freqRestart_ == 0
                || this->timeManager().episodeWillBeFinished()
                || this->timeManager().willBeFinished());
    }

    /*!
     * \brief Returns true if the current solution should be written to
     *        disk (i.e. as a VTK file)
     *
     * The default behavior is to write out the solution for
     * every time step. This function is intended to be overwritten by the
     * implementation.
     */
    bool shouldWriteOutput() const
    {
        return (this->timeManager().timeStepIndex() % freqOutput_ == 0
                || this->timeManager().episodeWillBeFinished()
                || this->timeManager().willBeFinished());
    }

private:
    bool verbose_;
    unsigned freqRestart_;
    unsigned freqOutput_;

    Scalar interfacePosY_;
    Scalar noDarcyX_;
    Scalar episodeLength_;
    Scalar initializationTime_;

    static constexpr Scalar eps_ = 1e-8;
};

} // end namespace Dumux

#endif // DUMUX_2CSTOKES2P2CPROBLEM_HH
