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
 * \file
 * \brief The problem class for the coupling of an isothermal two-component ZeroEq
 *        and an isothermal two-phase two-component Darcy model.
 */
#ifndef DUMUX_TWOCZEROEQTWOPTWOCPROBLEM_HH
#define DUMUX_TWOCZEROEQTWOPTWOCPROBLEM_HH

#include <dune/grid/multidomaingrid.hh>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/io/file/dgfparser.hh>

#include <dumux/material/fluidsystems/h2oairfluidsystem.hh>
#include <dumux/multidomain/common/multidomainproblem.hh>
#include <dumux/multidomain/2cstokes2p2c/2cstokes2p2clocaloperator.hh>
#include <dumux/multidomain/2cstokes2p2c/2cstokes2p2cpropertydefaults.hh>

#include "2czeroeq2p2cspatialparameters.hh"
#include "zeroeq2csubproblem.hh"
#include "2p2csubproblem.hh"

namespace Dumux
{
template <class TypeTag>
class TwoCZeroEqTwoPTwoCProblem;

namespace Properties
{
NEW_TYPE_TAG(TwoCZeroEqTwoPTwoCProblem, INHERITS_FROM(TwoCStokesTwoPTwoC));

// Set the grid type
#if HAVE_UG
SET_TYPE_PROP(TwoCZeroEqTwoPTwoCProblem, Grid, Dune::UGGrid<2>);
#elif HAVE_ALUGRID
SET_TYPE_PROP(TwoCZeroEqTwoPTwoCProblem, Grid, Dune::ALUGrid<2, 2, Dune::cube, Dune::nonconforming>);
#else
#error Requires UG or ALUGrid.
#endif

// Set the global problem
SET_TYPE_PROP(TwoCZeroEqTwoPTwoCProblem, Problem, TwoCZeroEqTwoPTwoCProblem<TypeTag>);

// Set the two sub-problems of the global problem
SET_TYPE_PROP(TwoCZeroEqTwoPTwoCProblem, SubDomain1TypeTag, TTAG(ZeroEq2cSubProblem));
SET_TYPE_PROP(TwoCZeroEqTwoPTwoCProblem, SubDomain2TypeTag, TTAG(TwoPTwoCSubProblem));

// Set the local coupling operator
SET_TYPE_PROP(TwoCZeroEqTwoPTwoCProblem, MultiDomainCouplingLocalOperator,
              Dumux::TwoCStokesTwoPTwoCLocalOperator<TypeTag>);

// Set the global problem in the context of the two sub-problems
SET_TYPE_PROP(ZeroEq2cSubProblem, MultiDomainTypeTag, TTAG(TwoCZeroEqTwoPTwoCProblem));
SET_TYPE_PROP(TwoPTwoCSubProblem, MultiDomainTypeTag, TTAG(TwoCZeroEqTwoPTwoCProblem));

// Set the other sub-problem for each of the two sub-problems
SET_TYPE_PROP(ZeroEq2cSubProblem, OtherSubDomainTypeTag, TTAG(TwoPTwoCSubProblem));
SET_TYPE_PROP(TwoPTwoCSubProblem, OtherSubDomainTypeTag, TTAG(ZeroEq2cSubProblem));

// Set the same spatial parameters for both sub-problems
SET_TYPE_PROP(ZeroEq2cSubProblem, SpatialParams, Dumux::TwoCZeroEqTwoPTwoCSpatialParams<TypeTag>);
SET_TYPE_PROP(TwoPTwoCSubProblem, SpatialParams, Dumux::TwoCZeroEqTwoPTwoCSpatialParams<TypeTag>);

// Set the fluid system
SET_PROP(TwoCZeroEqTwoPTwoCProblem, FluidSystem)
{
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef Dumux::FluidSystems::H2OAir<Scalar,
                                        Dumux::H2O<Scalar>,
                                        /*useComplexrelations=*/false> type;
};

// If SuperLU is not available, the UMFPack solver is used:
#ifdef HAVE_SUPERLU
SET_TYPE_PROP(TwoCZeroEqTwoPTwoCProblem, LinearSolver, SuperLUBackend<TypeTag>);
#else
SET_TYPE_PROP(TwoCZeroEqTwoPTwoCProblem, LinearSolver, UMFPackBackend<TypeTag>);
#endif
}

/*!
 * \brief The problem class for the coupling of an isothermal two-component ZeroEq (zeroeq2c)
 *        and an isothermal two-phase two-component Darcy model (2p2c).
 *
 *        The problem class for the coupling of a non-isothermal two-component ZeroEq (zeroeq2c)
 *        and an isothermal two-phase two-component Darcy model (2p2c).
 *        It uses the 2p2cCoupling model and the ZeroEq2ccoupling model and provides
 *        the problem specifications for common parameters of the two submodels.
 *        The initial and boundary conditions of the submodels are specified in the two subproblems,
 *        2p2csubproblem.hh and zeroeq2csubproblem.hh, which are accessible via the coupled problem.
 */
template <class TypeTag = TTAG(TwoCZeroEqTwoPTwoCProblem) >
class TwoCZeroEqTwoPTwoCProblem : public MultiDomainProblem<TypeTag>
{
    typedef MultiDomainProblem<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

    typedef typename GET_PROP_TYPE(TypeTag, MultiDomainGrid) MDGrid;
    typedef typename MDGrid::LeafGridView MDGridView;
    enum { dim = MDGridView::dimension };
    typedef Dune::FieldVector<Scalar, dim> GlobalPosition;

    typedef typename MDGrid::template Codim<0>::LeafIterator ElementIterator;

    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;

public:
    /*!
     * \brief The problem for the coupling of ZeroEq and Darcy flow
     *
     * \param mdGrid The multidomain grid
     * \param timeManager The time manager
     */
    TwoCZeroEqTwoPTwoCProblem(MDGrid &mdGrid,
                              TimeManager &timeManager)
        : ParentType(mdGrid, timeManager)
    {
        dtInit_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, TimeManager, DtInitial);
        episodeLength_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, TimeManager, EpisodeLength);

        // define location of the interface
        interfacePos_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, InterfacePos);
        noDarcyX_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, NoDarcyX);

        // define output options
        freqRestart_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, Output, FreqRestart);
        freqOutput_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, Output, FreqOutput);

        zeroeq2c_ = this->sdID1();
        twoPtwoC_ = this->sdID2();

        initializeGrid();

        // initialize the tables of the fluid system
        FluidSystem::init(/*tempMin=*/273.15, /*tempMax=*/373.15, /*numTemp=*/200,
                          /*pMin=*/1e4, /*pMax=*/2e5, /*numP=*/200);

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
        ElementIterator eendit = mdGrid.template leafend<0>();
        for (ElementIterator elementIt = mdGrid.template leafbegin<0>();
             elementIt != eendit; ++elementIt)
        {
            // this is required for parallelization
            // checks if element is within a partition
            if (elementIt->partitionType() != Dune::InteriorEntity)
                continue;

            GlobalPosition globalPos = elementIt->geometry().center();

            if (globalPos[1] > interfacePos_)
                mdGrid.addToSubDomain(zeroeq2c_,*elementIt);
            else
                if(globalPos[0] > noDarcyX_)
                    mdGrid.addToSubDomain(twoPtwoC_,*elementIt);
        }
        mdGrid.preUpdateSubDomains();
        mdGrid.updateSubDomains();
        mdGrid.postUpdateSubDomains();

        gridinfo(this->sdGrid1());
        gridinfo(this->sdGrid2());
    }

    /*!
     * \brief Called when the end of an simulation episode is reached.
     *
     * Typically a new episode should be started in this method.
     */
    void episodeEnd()
    { this->timeManager().startNextEpisode(episodeLength_); }

    //! \copydoc Dumux::CoupledProblem::shouldWriteRestartFile()
    bool shouldWriteRestartFile() const
    {
        return ( ((this->timeManager().timeStepIndex() > 0)
                  && (this->timeManager().timeStepIndex() % freqRestart_ == 0))
                // also write a restart file at the end of each episode
                || this->timeManager().episodeWillBeOver());
    }

    //! \copydoc Dumux::CoupledProblem::shouldWriteOutput()
    bool shouldWriteOutput() const
    {
        return ( ((this->timeManager().timeStepIndex() > 0)
                  && (this->timeManager().timeStepIndex() % freqOutput_ == 0))
                // also write a restart file at the end of each episode
                || this->timeManager().episodeWillBeOver());
    }

private:
    typename MDGrid::SubDomainType zeroeq2c_;
    typename MDGrid::SubDomainType twoPtwoC_;

    unsigned freqRestart_;
    unsigned freqOutput_;

    Scalar interfacePos_;
    Scalar noDarcyX_;
    Scalar episodeLength_;
    Scalar initializationTime_;
    Scalar dtInit_;
};

} // namespace Dumux

#endif // DUMUX_TWOCZEROEQTWOPTWOCPROBLEM_HH
