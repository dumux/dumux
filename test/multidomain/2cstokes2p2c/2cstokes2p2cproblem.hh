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
 * It uses the 2p2cCoupling model and the Stokes2ccoupling model and provides
 * the problem specifications for common parameters of the two submodels.
 * The initial and boundary conditions of the submodels are specified in the two subproblems,
 * 2p2csubproblem.hh and stokes2csubproblem.hh, which are accessible via the coupled problem.
 */

#ifndef DUMUX_2CSTOKES2P2CPROBLEM_HH
#define DUMUX_2CSTOKES2P2CPROBLEM_HH

#include <dune/common/float_cmp.hh>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/multidomaingrid.hh>

#include <dumux/material/fluidsystems/h2oair.hh>
#include <dumux/multidomain/2cstokes2p2c/localoperator.hh>
#include <dumux/multidomain/2cstokes2p2c/problem.hh>
#include <dumux/multidomain/2cstokes2p2c/propertydefaults.hh>

#ifdef HAVE_PARDISO
#include <dumux/linear/pardisobackend.hh>
#endif // HAVE_PARDISO

#include "2cstokes2p2cspatialparams.hh"
#include "stokes2csubproblem.hh"
#include "2p2csubproblem.hh"

namespace Dumux
{
template <class TypeTag>
class TwoCStokesTwoPTwoCTestProblem;

namespace Properties
{
NEW_TYPE_TAG(TwoCStokesTwoPTwoCTestProblem, INHERITS_FROM(TwoCStokesTwoPTwoC));

// Set the grid type
SET_TYPE_PROP(TwoCStokesTwoPTwoCTestProblem, Grid, Dune::YaspGrid<2, Dune::TensorProductCoordinates<typename GET_PROP_TYPE(TypeTag, Scalar), 2> >);

// Set the global problem
SET_TYPE_PROP(TwoCStokesTwoPTwoCTestProblem, Problem, Dumux::TwoCStokesTwoPTwoCTestProblem<TypeTag>);

// Set the local coupling operator
SET_TYPE_PROP(TwoCStokesTwoPTwoCTestProblem, MultiDomainCouplingLocalOperator,
              Dumux::TwoCStokesTwoPTwoCLocalOperator<TypeTag>);

// Set the two sub-problems of the global problem
SET_TYPE_PROP(TwoCStokesTwoPTwoCTestProblem, SubDomain1TypeTag, TTAG(Stokes2cSubProblem));
SET_TYPE_PROP(TwoCStokesTwoPTwoCTestProblem, SubDomain2TypeTag, TTAG(TwoPTwoCSubProblem));

// Set the global problem in the context of the two sub-problems
SET_TYPE_PROP(Stokes2cSubProblem, MultiDomainTypeTag, TTAG(TwoCStokesTwoPTwoCTestProblem));
SET_TYPE_PROP(TwoPTwoCSubProblem, MultiDomainTypeTag, TTAG(TwoCStokesTwoPTwoCTestProblem));

// Set the other sub-problem for each of the two sub-problems
SET_TYPE_PROP(Stokes2cSubProblem, OtherSubDomainTypeTag, TTAG(TwoPTwoCSubProblem));
SET_TYPE_PROP(TwoPTwoCSubProblem, OtherSubDomainTypeTag, TTAG(Stokes2cSubProblem));

// Set the spatial parameters used for the problems
SET_TYPE_PROP(TwoPTwoCSubProblem, SpatialParams, Dumux::TwoCStokesTwoPTwoCSpatialParams<TypeTag>);

// Set the fluid system to use simple relations (last argument)
SET_TYPE_PROP(TwoCStokesTwoPTwoCTestProblem, FluidSystem,
              FluidSystems::H2OAir<typename GET_PROP_TYPE(TypeTag, Scalar),
                                   Dumux::H2O<typename GET_PROP_TYPE(TypeTag, Scalar)>, false>);

// if you do not have PARDISO, the SuperLU solver is used:
#ifdef HAVE_PARDISO
SET_TYPE_PROP(TwoCStokesTwoPTwoCTestProblem, LinearSolver, PardisoBackend<TypeTag>);
#else
SET_TYPE_PROP(TwoCStokesTwoPTwoCTestProblem, LinearSolver, SuperLUBackend<TypeTag>);
#endif // HAVE_PARDISO
}

/*!
 * \ingroup ImplicitTestProblems
 * \ingroup TwoPTwoCStokesTwoCModel
 * \brief The problem class for the coupling of an isothermal two-component Stokes
 *        and an isothermal two-phase two-component Darcy model.
 *
 * The problem class for the coupling of an isothermal two-component Stokes (stokes2c)
 * and an isothermal two-phase two-component Darcy model (2p2c).
 * It uses the 2p2cCoupling model and the Stokes2ccoupling model and provides
 * the problem specifications for common parameters of the two submodels.
 * The initial and boundary conditions of the submodels are specified in the two subproblems,
 * 2p2csubproblem.hh and stokes2csubproblem.hh, which are accessible via the coupled problem.
 */
template <class TypeTag = TTAG(TwoCStokesTwoPTwoCTestProblem) >
class TwoCStokesTwoPTwoCTestProblem : public TwoCStokesTwoPTwoCProblem<TypeTag>
{
    typedef TwoCStokesTwoPTwoCTestProblem<TypeTag> ThisType;
    typedef TwoCStokesTwoPTwoCProblem<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

    typedef typename GET_PROP_TYPE(TypeTag, MultiDomainGrid) MDGrid;
    typedef typename MDGrid::LeafGridView MDGridView;
    enum { dim = MDGridView::dimension };
    typedef Dune::FieldVector<Scalar, dim> GlobalPosition;

    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

public:
    /*!
     * \brief The problem for the coupling of Stokes and Darcy flow
     *
     * \param timeManager The time manager
     * \param gridView The grid view
     */
    template<class GridView>
    TwoCStokesTwoPTwoCTestProblem(TimeManager &timeManager,
                                  GridView gridView)
    : ParentType(timeManager, gridView)
    {
        interfacePosY_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, InterfacePosY);
        noDarcyX_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, NoDarcyX);
        episodeLength_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, TimeManager, EpisodeLength);
        initializationTime_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, TimeManager, InitTime);

        // define output options
        freqRestart_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, Output, FreqRestart);
        freqOutput_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, Output, FreqOutput);

        stokes2c_ = this->sdID1();
        twoPtwoC_ = this->sdID2();

        initializeGrid();

        // initialize the tables of the fluid system
        FluidSystem::init(/*tempMin=*/273.15, /*tempMax=*/373.15, /*numTemp=*/200,
                          /*pMin=*/1e4, /*pMax=*/2e5, /*numP=*/200);

        if (initializationTime_ > 0.0)
            this->timeManager().startNextEpisode(initializationTime_);
        else
            this->timeManager().startNextEpisode(episodeLength_);
    }

    /*!
     * \brief Initialization of the grids
     *
     * This function splits the multidomain grid in the two
     * individual subdomain grids and takes care of parallelization.
     */
    void initializeGrid()
    {
        MDGrid& mdGrid = this->mdGrid();
        mdGrid.startSubDomainMarking();

        // subdivide grid in two subdomains
        for (const auto& element : elements(mdGrid.leafGridView()))
        {
            // this is required for parallelization, checks if element is within a partition
            if (element.partitionType() != Dune::InteriorEntity)
                continue;

            GlobalPosition globalPos = element.geometry().center();

            if (globalPos[1] > interfacePosY_)
                mdGrid.addToSubDomain(stokes2c_,element);
            else
                if(globalPos[0] > noDarcyX_)
                    mdGrid.addToSubDomain(twoPtwoC_,element);
        }
        mdGrid.preUpdateSubDomains();
        mdGrid.updateSubDomains();
        mdGrid.postUpdateSubDomains();

        gridinfo(this->sdGrid1());
        gridinfo(this->sdGrid2());
    }

    /*!
     * \brief Called by the time manager after the time integration to
     *        do some post processing on the solution.
     */
    void postTimeStep()
    {
        // call the postTimeStep function of the subproblems
        this->sdProblem1().postTimeStep();
        this->sdProblem2().postTimeStep();
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
                || this->timeManager().episodeWillBeOver());
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
                || this->timeManager().episodeWillBeOver());
    }

private:
    typename MDGrid::SubDomainIndex stokes2c_;
    typename MDGrid::SubDomainIndex twoPtwoC_;

    unsigned freqRestart_;
    unsigned freqOutput_;

    Scalar interfacePosY_;
    Scalar noDarcyX_;
    Scalar episodeLength_;
    Scalar initializationTime_;
};

} // namespace Dumux

#endif // DUMUX_2CSTOKES2P2CPROBLEM_HH
