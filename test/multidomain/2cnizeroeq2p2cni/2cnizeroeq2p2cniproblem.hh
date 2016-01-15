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
 * \brief The problem which couples a non-isothermal two-component ZeroEq
 *        and a non-isothermal two-phase two-component Darcy model.
 */
#ifndef DUMUX_TWOCNIZEROEQTWOPTWOCNIPROBLEM_HH
#define DUMUX_TWOCNIZEROEQTWOPTWOCNIPROBLEM_HH

#include <dune/grid/multidomaingrid.hh>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/io/file/dgfparser.hh>

#include <dumux/material/fluidsystems/h2oairfluidsystem.hh>
#include <dumux/multidomain/common/multidomainproblem.hh>
#include <dumux/multidomain/2cnistokes2p2cni/2cnistokes2p2cnilocaloperator.hh>
#include <dumux/multidomain/2cnistokes2p2cni/2cnistokes2p2cnipropertydefaults.hh>

#include "2cnizeroeq2p2cnispatialparameters.hh"
#include "zeroeq2cnisubproblem.hh"
#include "2p2cnisubproblem.hh"

namespace Dumux
{
template <class TypeTag>
class TwoCNIZeroEqTwoPTwoCNIProblem;

namespace Properties
{
NEW_TYPE_TAG(TwoCNIZeroEqTwoPTwoCNIProblem, INHERITS_FROM(TwoCNIStokesTwoPTwoCNI));

// Set the grid type
SET_TYPE_PROP(TwoCNIZeroEqTwoPTwoCNIProblem, Grid, Dune::YaspGrid<2, Dune::TensorProductCoordinates<typename GET_PROP_TYPE(TypeTag, Scalar), 2> >);

// Set the global problem
SET_TYPE_PROP(TwoCNIZeroEqTwoPTwoCNIProblem, Problem, TwoCNIZeroEqTwoPTwoCNIProblem<TypeTag>);

// Set the local coupling operator
SET_TYPE_PROP(TwoCNIZeroEqTwoPTwoCNIProblem, MultiDomainCouplingLocalOperator,
              Dumux::TwoCNIStokesTwoPTwoCNILocalOperator<TypeTag>);

// Set the two sub-problems of the global problem
SET_TYPE_PROP(TwoCNIZeroEqTwoPTwoCNIProblem, SubDomain1TypeTag, TTAG(ZeroEq2cniSubProblem));
SET_TYPE_PROP(TwoCNIZeroEqTwoPTwoCNIProblem, SubDomain2TypeTag, TTAG(TwoPTwoCNISubProblem));

// Set the global problem in the context of the two sub-problems
SET_TYPE_PROP(ZeroEq2cniSubProblem, MultiDomainTypeTag, TTAG(TwoCNIZeroEqTwoPTwoCNIProblem));
SET_TYPE_PROP(TwoPTwoCNISubProblem, MultiDomainTypeTag, TTAG(TwoCNIZeroEqTwoPTwoCNIProblem));

// Set the other sub-problem for each of the two sub-problems
SET_TYPE_PROP(ZeroEq2cniSubProblem, OtherSubDomainTypeTag, TTAG(TwoPTwoCNISubProblem));
SET_TYPE_PROP(TwoPTwoCNISubProblem, OtherSubDomainTypeTag, TTAG(ZeroEq2cniSubProblem));

// Set the same spatial parameters for both sub-problems
SET_TYPE_PROP(TwoPTwoCNISubProblem, SpatialParams, Dumux::TwoCNIZeroEqTwoPTwoCNISpatialParams<TypeTag>);

// Set the fluid system to use complex relations (last argument)
SET_TYPE_PROP(TwoCNIZeroEqTwoPTwoCNIProblem, FluidSystem,
              FluidSystems::H2OAir<typename GET_PROP_TYPE(TypeTag, Scalar),
                                   Dumux::H2O<typename GET_PROP_TYPE(TypeTag, Scalar)>, true>);

// If SuperLU is not available, the UMFPack solver is used:
#ifdef HAVE_SUPERLU
SET_TYPE_PROP(TwoCNIZeroEqTwoPTwoCNIProblem, LinearSolver, SuperLUBackend<TypeTag>);
#else
SET_TYPE_PROP(TwoCNIZeroEqTwoPTwoCNIProblem, LinearSolver, UMFPackBackend<TypeTag>);
#endif
}

/*!
 * \ingroup ImplicitTestProblems
 * \ingroup TwoPTwoCNIZeroEqTwoCNIModel
 *
 * \brief The problem which couples a non-isothermal two-component ZeroEq (zeroeq2cni)
 *        and a non-isothermal two-phase two-component Darcy model (2p2cni).
 *
 * It uses the multidomain problem and specifies parameters for the two submodels.
 * The initial and boundary conditions of the submodels are specified in the two subproblems,
 * 2p2csubproblem.hh and zeroeq2csubproblem.hh, which are accessible via the coupled problem.
 */
template <class TypeTag = TTAG(TwoCNIZeroEqTwoPTwoCNIProblem) >
class TwoCNIZeroEqTwoPTwoCNIProblem : public MultiDomainProblem<TypeTag>
{
    typedef MultiDomainProblem<TypeTag> ParentType;

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
    TwoCNIZeroEqTwoPTwoCNIProblem(TimeManager &timeManager,
                                  GridView gridView)
    : ParentType(timeManager, gridView)
    {
        dtInit_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, TimeManager, DtInitial);
        episodeLength_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, TimeManager, EpisodeLength);

        // define location of the interface
        interfacePosY_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, InterfacePosY);
        noDarcyX1_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, NoDarcyX1);
        noDarcyX2_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, NoDarcyX2);

        // define output options
        freqRestart_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, Output, FreqRestart);
        freqOutput_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, Output, FreqOutput);

        zeroeq2cni_ = this->sdID1();
        twoPtwoCNI_ = this->sdID2();

        initializeGrid();

        // initialize the tables of the fluid system
        FluidSystem::init();
//         FluidSystem::init(/*tempMin=*/273.15, /*tempMax=*/373.15, /*numTemp=*/200,
//                           /*pMin=*/1e4, /*pMax=*/2e5, /*numP=*/200);

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
        for (const auto& element : Dune::elements(mdGrid.leafGridView()))
        {
            // this is required for parallelization checks if element is within a partition
            if (element.partitionType() != Dune::InteriorEntity)
                continue;

            GlobalPosition globalPos = element.geometry().center();

            if (globalPos[1] > interfacePosY_)
                mdGrid.addToSubDomain(zeroeq2cni_,element);
            else
                if(globalPos[0] > noDarcyX1_ && globalPos[0] < noDarcyX2_)
                    mdGrid.addToSubDomain(twoPtwoCNI_,element);
        }
        mdGrid.preUpdateSubDomains();
        mdGrid.updateSubDomains();
        mdGrid.postUpdateSubDomains();

        gridinfo(this->sdGrid1());
        gridinfo(this->sdGrid2());
    }

    //! \copydoc Dumux::ImplicitProblem::episodeEnd()
    void episodeEnd()
    { this->timeManager().startNextEpisode(episodeLength_); }

    //! \copydoc Dumux::ImplicitProblem::shouldWriteRestartFile()
    bool shouldWriteRestartFile() const
    {
        return (((this->timeManager().timeStepIndex() > 0)
                  && (this->timeManager().timeStepIndex() % freqRestart_ == 0))
                || this->timeManager().episodeWillBeOver());
    }

    //! \copydoc Dumux::ImplicitProblem::shouldWriteOutput()
    bool shouldWriteOutput() const
    {
        return (this->timeManager().timeStepIndex() % freqOutput_ == 0
                || this->timeManager().episodeWillBeOver());
    }

private:
    typename MDGrid::SubDomainType zeroeq2cni_;
    typename MDGrid::SubDomainType twoPtwoCNI_;

    unsigned freqRestart_;
    unsigned freqOutput_;

    Scalar interfacePosY_;
    Scalar noDarcyX1_;
    Scalar noDarcyX2_;
    Scalar episodeLength_;
    Scalar dtInit_;
};

} //end namespace

#endif // DUMUX_TWOCNIZEROEQTWOPTWOCNIPROBLEM_HH
