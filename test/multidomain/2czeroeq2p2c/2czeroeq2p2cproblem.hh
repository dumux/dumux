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
 * \brief The problem which couples an isothermal two-component ZeroEq
 *        and an isothermal two-phase two-component Darcy model.
 */
#ifndef DUMUX_TWOCZEROEQTWOPTWOCPROBLEM_HH
#define DUMUX_TWOCZEROEQTWOPTWOCPROBLEM_HH

#include <dune/grid/multidomaingrid.hh>
#include <dune/grid/common/gridinfo.hh>

#include <dumux/material/fluidsystems/h2oair.hh>
#include <dumux/multidomain/2cstokes2p2c/localoperator.hh>
#include <dumux/multidomain/2cstokes2p2c/problem.hh>
#include <dumux/multidomain/2cstokes2p2c/propertydefaults.hh>

#include "2czeroeq2p2cspatialparameters.hh"
#include "zeroeq2csubproblem.hh"
#include "2p2csubproblem.hh"

namespace Dumux
{
template <class TypeTag>
class TwoCZeroEqTwoPTwoCTestProblem;

namespace Properties
{
NEW_TYPE_TAG(TwoCZeroEqTwoPTwoCTestProblem, INHERITS_FROM(TwoCStokesTwoPTwoC));

// Set the grid type
SET_TYPE_PROP(TwoCZeroEqTwoPTwoCTestProblem, Grid, Dune::YaspGrid<2, Dune::TensorProductCoordinates<typename GET_PROP_TYPE(TypeTag, Scalar), 2> >);

// Set the global problem
SET_TYPE_PROP(TwoCZeroEqTwoPTwoCTestProblem, Problem, TwoCZeroEqTwoPTwoCTestProblem<TypeTag>);

// Set the two sub-problems of the global problem
SET_TYPE_PROP(TwoCZeroEqTwoPTwoCTestProblem, SubDomain1TypeTag, TTAG(ZeroEq2cSubProblem));
SET_TYPE_PROP(TwoCZeroEqTwoPTwoCTestProblem, SubDomain2TypeTag, TTAG(TwoPTwoCSubProblem));

// Set the local coupling operator
SET_TYPE_PROP(TwoCZeroEqTwoPTwoCTestProblem, MultiDomainCouplingLocalOperator,
              Dumux::TwoCStokesTwoPTwoCLocalOperator<TypeTag>);

// Set the global problem in the context of the two sub-problems
SET_TYPE_PROP(ZeroEq2cSubProblem, MultiDomainTypeTag, TTAG(TwoCZeroEqTwoPTwoCTestProblem));
SET_TYPE_PROP(TwoPTwoCSubProblem, MultiDomainTypeTag, TTAG(TwoCZeroEqTwoPTwoCTestProblem));

// Set the other sub-problem for each of the two sub-problems
SET_TYPE_PROP(ZeroEq2cSubProblem, OtherSubDomainTypeTag, TTAG(TwoPTwoCSubProblem));
SET_TYPE_PROP(TwoPTwoCSubProblem, OtherSubDomainTypeTag, TTAG(ZeroEq2cSubProblem));

// Set the same spatial parameters for both sub-problems
SET_TYPE_PROP(TwoPTwoCSubProblem, SpatialParams, Dumux::TwoCZeroEqTwoPTwoCSpatialParams<TypeTag>);

// Set the fluid system to use simple relations (last argument)
SET_TYPE_PROP(TwoCZeroEqTwoPTwoCTestProblem, FluidSystem,
              FluidSystems::H2OAir<typename GET_PROP_TYPE(TypeTag, Scalar),
                                   Dumux::H2O<typename GET_PROP_TYPE(TypeTag, Scalar)>, false>);

// If SuperLU is not available, the UMFPack solver is used:
#ifdef HAVE_SUPERLU
SET_TYPE_PROP(TwoCZeroEqTwoPTwoCTestProblem, LinearSolver, SuperLUBackend<TypeTag>);
#else
SET_TYPE_PROP(TwoCZeroEqTwoPTwoCTestProblem, LinearSolver, UMFPackBackend<TypeTag>);
#endif
}


/*!
 * \ingroup TwoPTwoCZeroEqTwoCModel
 * \ingroup ImplicitTestProblems
 * \brief The problem which couples an isothermal two-component ZeroEq (zeroeq2c)
 *        and an isothermal two-phase two-component Darcy model (2p2c).
 *
 * It uses the multidomain problem and specifies parameters for the two submodels.
 * The initial and boundary conditions of the submodels are specified in the two subproblems,
 * 2p2csubproblem.hh and zeroeq2csubproblem.hh, which are accessible via the coupled problem.
 */
template <class TypeTag = TTAG(TwoCZeroEqTwoPTwoCTestProblem) >
class TwoCZeroEqTwoPTwoCTestProblem : public TwoCStokesTwoPTwoCProblem<TypeTag>
{
    typedef TwoCStokesTwoPTwoCProblem<TypeTag> ParentType;

    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, TimeManager) TimeManager;

    typedef typename GET_PROP_TYPE(TypeTag, GridCreator) GridCreator;
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
    TwoCZeroEqTwoPTwoCTestProblem(TimeManager &timeManager,
                                  GridView gridView)
    : ParentType(timeManager, gridView)
    {
        dtInit_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, TimeManager, DtInitial);
        episodeLength_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, TimeManager, EpisodeLength);

        // define location of the interface
        interfacePosY_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, InterfacePosY);
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
        for (const auto& element : elements(mdGrid.leafGridView()))
        {
            // this is required for parallelization
            // checks if element is within a partition
            if (element.partitionType() != Dune::InteriorEntity)
                continue;

            GlobalPosition globalPos = element.geometry().center();

            if (globalPos[1] > interfacePosY_)
                mdGrid.addToSubDomain(zeroeq2c_,element);
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
     * \brief Called when the end of an simulation episode is reached.
     *
     * Typically a new episode should be started in this method.
     */
    void episodeEnd()
    { this->timeManager().startNextEpisode(episodeLength_); }

    //! \copydoc Dumux::ImplicitProblem::shouldWriteRestartFile()
    bool shouldWriteRestartFile() const
    {
        return ( ((this->timeManager().timeStepIndex() > 0)
                  && (this->timeManager().timeStepIndex() % freqRestart_ == 0))
                // also write a restart file at the end of each episode
                || this->timeManager().episodeWillBeOver());
    }

    //! \copydoc Dumux::ImplicitProblem::shouldWriteOutput()
    bool shouldWriteOutput() const
    {
        return (this->timeManager().timeStepIndex() % freqOutput_ == 0
                || this->timeManager().episodeWillBeOver());
    }

private:
    typename MDGrid::SubDomainType zeroeq2c_;
    typename MDGrid::SubDomainType twoPtwoC_;

    unsigned freqRestart_;
    unsigned freqOutput_;

    Scalar interfacePosY_;
    Scalar noDarcyX_;
    Scalar episodeLength_;
    Scalar initializationTime_;
    Scalar dtInit_;
};

} // namespace Dumux

#endif // DUMUX_TWOCZEROEQTWOPTWOCPROBLEM_HH
