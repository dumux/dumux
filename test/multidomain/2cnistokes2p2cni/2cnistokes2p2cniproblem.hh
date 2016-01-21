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
 * \brief The problem class for the coupling of a non-isothermal two-component Stokes
 *        and a non-isothermal two-phase two-component Darcy model.
 *
 * The problem class for the coupling of a non-isothermal two-component Stokes (stokes2cn)
 * and a non-isothermal two-phase two-component Darcy model (2p2cni).
 * It uses the 2p2cniCoupling model and the Stokes2cnicoupling model and provides
 * the problem specifications for common parameters of the two submodels.
 * The initial and boundary conditions of the submodels are specified in the two subproblems,
 * 2p2cnisubproblem.hh and stokes2cnisubproblem.hh, which are accessible via the coupled problem.
 */

#ifndef DUMUX_2CNISTOKES2P2CNIPROBLEM_HH
#define DUMUX_2CNISTOKES2P2CNIPROBLEM_HH

#include <dune/common/float_cmp.hh>
#include <dune/grid/common/gridinfo.hh>
#include <dune/grid/multidomaingrid.hh>
#include <dune/grid/io/file/dgfparser.hh>

#include <dumux/material/fluidsystems/h2oairfluidsystem.hh>
#include <dumux/multidomain/common/problem.hh>
#include <dumux/multidomain/2cnistokes2p2cni/localoperator.hh>
#include <dumux/multidomain/2cnistokes2p2cni/propertydefaults.hh>

#include <dumux/linear/seqsolverbackend.hh>
#ifdef HAVE_PARDISO
#include <dumux/linear/pardisobackend.hh>
#endif // HAVE_PARDISO

#include "stokes2cnisubproblem.hh"
#include "2p2cnisubproblem.hh"

namespace Dumux
{
template <class TypeTag>
class TwoCNIStokesTwoPTwoCNIProblem;

namespace Properties
{
NEW_TYPE_TAG(TwoCNIStokesTwoPTwoCNIProblem, INHERITS_FROM(TwoCNIStokesTwoPTwoCNI));

// Set the grid type
SET_TYPE_PROP(TwoCNIStokesTwoPTwoCNIProblem, Grid, Dune::YaspGrid<2, Dune::TensorProductCoordinates<typename GET_PROP_TYPE(TypeTag, Scalar), 2> >);

// Set the global problem
SET_TYPE_PROP(TwoCNIStokesTwoPTwoCNIProblem, Problem, TwoCNIStokesTwoPTwoCNIProblem<TypeTag>);

// Set the local coupling operator
SET_TYPE_PROP(TwoCNIStokesTwoPTwoCNIProblem, MultiDomainCouplingLocalOperator,
              Dumux::TwoCNIStokesTwoPTwoCNILocalOperator<TypeTag>);

// Set the two sub-problems of the global problem
SET_TYPE_PROP(TwoCNIStokesTwoPTwoCNIProblem, SubDomain1TypeTag, TTAG(Stokes2cniSubProblem));
SET_TYPE_PROP(TwoCNIStokesTwoPTwoCNIProblem, SubDomain2TypeTag, TTAG(TwoPTwoCNISubProblem));

// Set the global problem in the context of the two sub-problems
SET_TYPE_PROP(Stokes2cniSubProblem, MultiDomainTypeTag, TTAG(TwoCNIStokesTwoPTwoCNIProblem));
SET_TYPE_PROP(TwoPTwoCNISubProblem, MultiDomainTypeTag, TTAG(TwoCNIStokesTwoPTwoCNIProblem));

// Set the other sub-problem for each of the two sub-problems
SET_TYPE_PROP(Stokes2cniSubProblem, OtherSubDomainTypeTag, TTAG(TwoPTwoCNISubProblem));
SET_TYPE_PROP(TwoPTwoCNISubProblem, OtherSubDomainTypeTag, TTAG(Stokes2cniSubProblem));

// Set the spatial parameters used for the problems
SET_TYPE_PROP(TwoPTwoCNISubProblem, SpatialParams, Dumux::TwoCNIStokesTwoPTwoCNISpatialParams<TypeTag>);

// Set the fluid system to use complex relations (last argument)
SET_TYPE_PROP(TwoCNIStokesTwoPTwoCNIProblem, FluidSystem,
              FluidSystems::H2OAir<typename GET_PROP_TYPE(TypeTag, Scalar),
                                   Dumux::H2O<typename GET_PROP_TYPE(TypeTag, Scalar)>, true>);

#ifdef HAVE_PARDISO
SET_TYPE_PROP(TwoCNIStokesTwoPTwoCNIProblem, LinearSolver, PardisoBackend<TypeTag>);
#else
SET_TYPE_PROP(TwoCNIStokesTwoPTwoCNIProblem, LinearSolver, SuperLUBackend<TypeTag>);
#endif // HAVE_PARDISO
}

/*!
 * \ingroup ImplicitTestProblems
 * \ingroup TwoPTwoCNIStokesTwoCNIModel
 * \brief The problem class for the coupling of a non-isothermal two-component Stokes
 *        and a non-isothermal two-phase two-component Darcy model.
 *
 * The problem class for the coupling of a non-isothermal two-component Stokes (stokes2cn)
 * and a non-isothermal two-phase two-component Darcy model (2p2cni).
 * It uses the 2p2cniCoupling model and the Stokes2cnicoupling model and provides
 * the problem specifications for common parameters of the two submodels.
 * The initial and boundary conditions of the submodels are specified in the two subproblems,
 * 2p2cnisubproblem.hh and stokes2cnisubproblem.hh, which are accessible via the coupled problem.
 */
template <class TypeTag = TTAG(TwoCNIStokesTwoPTwoCNIProblem) >
class TwoCNIStokesTwoPTwoCNIProblem : public MultiDomainProblem<TypeTag>
{
    typedef TwoCNIStokesTwoPTwoCNIProblem<TypeTag> ThisType;
    typedef MultiDomainProblem<TypeTag> ParentType;

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
    TwoCNIStokesTwoPTwoCNIProblem(TimeManager &timeManager,
                                  GridView gridView)
    : ParentType(timeManager, gridView)
    {
        interfacePosY_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, InterfacePosY);
        noDarcyX_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, NoDarcyX);
        episodeLength_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, TimeManager, EpisodeLength);
        initializationTime_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, TimeManager, InitTime);
        dtInit_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, TimeManager, DtInitial);

        // define output options
        freqRestart_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, Output, FreqRestart);
        freqOutput_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, Output, FreqOutput);

        stokes2cni_ = this->sdID1();
        twoPtwoCNI_ = this->sdID2();

        initializeGrid();

        // initialize the tables of the fluid system
        FluidSystem::init(/*tempMin=*/273.15, /*tempMax=*/373.15, /*numTemp=*/200,
                          /*pMin=*/1e3, /*pMax=*/2e5, /*numP=*/200);

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
        for (const auto& element : Dune::elements(mdGrid.leafGridView()))
        {
            // this is required for parallelization, checks if element is within a partition
            if (element.partitionType() != Dune::InteriorEntity)
                continue;

            GlobalPosition globalPos = element.geometry().center();

            if (globalPos[1] > interfacePosY_)
                mdGrid.addToSubDomain(stokes2cni_,element);
            else
                if(globalPos[0] > noDarcyX_)
                    mdGrid.addToSubDomain(twoPtwoCNI_,element);
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
    typename MDGrid::SubDomainIndex stokes2cni_;
    typename MDGrid::SubDomainIndex twoPtwoCNI_;

    unsigned freqRestart_;
    unsigned freqOutput_;

    Scalar interfacePosY_;
    Scalar noDarcyX_;
    Scalar episodeLength_;
    Scalar initializationTime_;
    Scalar dtInit_;
};

} // namespace Dumux

#endif // DUMUX_2CNISTOKES2P2CNIPROBLEM_HH
