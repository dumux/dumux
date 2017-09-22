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
 * \brief A windtunnel test problem (taken from Fetzer2017c)
 */
#ifndef DUMUX_WINDTUNNEL_TEST_PROBLEM_HH
#define DUMUX_WINDTUNNEL_TEST_PROBLEM_HH

#include "windtunnel_2p2cnisubproblem.hh"
#include "windtunnel_navierstokes2cnisubproblem.hh"

#include <dumux/multidomain/problem.hh>
#include <dumux/multidomain/couplingmanager.hh>
#include <dumux/multidomain/map.hh>

#include <dumux/multidomain/staggered-ccfv/properties.hh>
#include <dumux/multidomain/staggered-ccfv/model.hh>
#include <dumux/multidomain/staggered-ccfv/assembler.hh>
#include <dumux/multidomain/staggered-ccfv/newtoncontroller.hh>
#include <dumux/multidomain/staggered-ccfv/subproblemlocaljacobian.hh>

namespace Dumux
{
template <class TypeTag>
class WindTunnelProblem;

namespace Properties
{
NEW_TYPE_TAG(WindTunnelProblem, INHERITS_FROM(MultiDomain, CouplingStokesStaggeredModel));

// Set the problem property
SET_TYPE_PROP(WindTunnelProblem, Problem, Dumux::WindTunnelProblem<TypeTag>);

// Set the coupling manager
SET_TYPE_PROP(WindTunnelProblem, CouplingManager, Dumux::CouplingManagerStokesDarcy<TypeTag>);

//////////////////////////////////////////////////////////////////////////
// Set the two sub-problems of the global problem
SET_TYPE_PROP(WindTunnelProblem, DarcyProblemTypeTag, TTAG(DarcySubProblem));
SET_TYPE_PROP(WindTunnelProblem, StokesProblemTypeTag, TTAG(StokesSubProblem));
////////////////////////////////////////////////////////////////////////////

// publish this problem in the sub problems
SET_TYPE_PROP(DarcySubProblem, GlobalProblemTypeTag, TTAG(WindTunnelProblem));
SET_TYPE_PROP(StokesSubProblem, GlobalProblemTypeTag, TTAG(WindTunnelProblem));

// The subproblems inherit the parameter tree from this problem
SET_PROP(DarcySubProblem, ParameterTree) : GET_PROP(TTAG(WindTunnelProblem), ParameterTree) {};
SET_PROP(StokesSubProblem, ParameterTree) : GET_PROP(TTAG(WindTunnelProblem), ParameterTree) {};

// Set the fluid system to use complex relations (last argument)
SET_TYPE_PROP(WindTunnelProblem, FluidSystem,FluidSystems::H2OAir<typename GET_PROP_TYPE(TypeTag, Scalar)>);
SET_TYPE_PROP(DarcySubProblem, FluidSystem, typename GET_PROP_TYPE(TTAG(WindTunnelProblem), FluidSystem));
SET_TYPE_PROP(StokesSubProblem, FluidSystem, typename GET_PROP_TYPE(TTAG(WindTunnelProblem), FluidSystem));

// Use a smaller FluidSystem table
NEW_PROP_TAG(ProblemUseSmallFluidSystemTable);
SET_BOOL_PROP(WindTunnelProblem, ProblemUseSmallFluidSystemTable, false);

// Define whether mole(true) or mass (false) fractions are used
SET_BOOL_PROP(WindTunnelProblem, UseMoles, false);
SET_BOOL_PROP(DarcySubProblem, UseMoles, GET_PROP_VALUE(TTAG(WindTunnelProblem), UseMoles));
SET_BOOL_PROP(StokesSubProblem, UseMoles, GET_PROP_VALUE(TTAG(WindTunnelProblem), UseMoles));

// Set the output frequency
NEW_PROP_TAG(OutputFreqOutput);
SET_INT_PROP(WindTunnelProblem, OutputFreqOutput, 5);

// Set the default episode length
NEW_PROP_TAG(TimeManagerEpisodeLength);
SET_INT_PROP(WindTunnelProblem, TimeManagerEpisodeLength, 43200);

// TODO: Fetzer2017c: GridInterfaceProfile = flat, GridPorousMediumBoxType = rectangular

NEW_PROP_TAG(DarcyToStokesMapValue);
SET_TYPE_PROP(WindTunnelProblem, DarcyToStokesMapValue, Dumux::DarcyToStokesMap<TypeTag>);

}//end namespace properties

template <class TypeTag>
class WindTunnelProblem : public MultiDomainProblem<TypeTag>
{
    using ParentType = MultiDomainProblem<TypeTag>;
    using TimeManager = typename GET_PROP_TYPE(TypeTag, TimeManager);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);

    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    enum { dim = GridView::dimension };
    using GlobalPosition = Dune::FieldVector<Scalar, dim>;

    // obtain the type tags of the sub problems
    using StokesProblemTypeTag = typename GET_PROP_TYPE(TypeTag, StokesProblemTypeTag);
    using DarcyProblemTypeTag = typename GET_PROP_TYPE(TypeTag, DarcyProblemTypeTag);

    // obtain types from the sub problem type tags
    using StokesProblem = typename GET_PROP_TYPE(StokesProblemTypeTag, Problem);
    using DarcyProblem = typename GET_PROP_TYPE(DarcyProblemTypeTag, Problem);

    using StokesGridView = typename GET_PROP_TYPE(StokesProblemTypeTag, GridView);
    using DarcyGridView = typename GET_PROP_TYPE(DarcyProblemTypeTag, GridView);

    // TODO
//    constexpr static unsigned int stokesSubDomainIdx = MultiDomainIndices::stokesSubDomainIdx;
//    constexpr static unsigned int darcySubDomainIdx = MultiDomainIndices::darcySubDomainIdx;

public:
    /*!
     * \brief Windtunnel problem
     *
     * \param timeManager The TimeManager which is used by the simulation
     * \stokesGridView The GridView of the Stokes domain
     * \darcyGridView The GridView of the Darcy domain
     */
    WindTunnelProblem(TimeManager &timeManager, const StokesGridView &stokesGridView, const DarcyGridView &darcygridView)
    : ParentType(timeManager, stokesGridView, darcygridView)
    {
    // initialize the tables of the fluid system
    if (GET_PARAM_FROM_GROUP(TypeTag, bool, Problem, UseSmallFluidSystemTable))
    {
        FluidSystem::init(/*tempMin=*/273.15, /*tempMax=*/323.15, /*numTemp=*/50,
                          /*pMin=*/5e4, /*pMax=*/1.5e5, /*numP=*/100);
    }
    else
    {
        FluidSystem::init(/*tempMin=*/273.15, /*tempMax=*/343.15, /*numTemp=*/140,
                          /*pMin=*/5e4, /*pMax=*/1e7, /*numP=*/995);
    }

    // spatial parameter stuff, if requested
    this->sdProblemDarcy().spatialParams().plotMaterialLaw();

    freqOutput_ = GET_PARAM_FROM_GROUP(TypeTag, int, Output, FreqOutput);
    episodeLength_ = GET_PARAM_FROM_GROUP(TypeTag, Scalar, TimeManager, EpisodeLength);
    this->timeManager().startNextEpisode(episodeLength_);
    }

    // TODO necessary?
//    //! \copydoc Dumux::ImplicitProblem::name()
//    const std::string name() const
//    { return GET_RUNTIME_PARAM(TypeTag, std::string, Output.Name) + "_staggered"; }

    //! \copydoc Dumux::CoupledProblem::episodeEnd()
    void episodeEnd()
    { this->timeManager().startNextEpisode(episodeLength_); }

    //! \copydoc Dumux::CoupledProblem::shouldWriteOutput()
    bool shouldWriteOutput() const
    {
        return this->timeManager().timeStepIndex() % freqOutput_ == 0
               || this->timeManager().episodeWillBeFinished()
               || this->timeManager().willBeFinished();
    }

    // TODO initializeGrid not needed here, --> where?
//    /*!
//     * \brief Initialization multi-domain and the sub-domain grids
//     */
//    void initializeGrid()
//    {
//        std::string porousMediumBoxType = GET_PARAM_FROM_GROUP(TypeTag, std::string, Grid, PorousMediumBoxType);
//        this->mdGrid().startSubDomainMarking();
//
//        for (auto eIt = this->mdGrid().template leafbegin<0>();
//              eIt != this->mdGrid().template leafend<0>(); ++eIt)
//        {
//            auto globalPos = eIt->geometry().center();
//            // only for interior entities, required for parallelization
//            if (eIt->partitionType() == Dune::InteriorEntity)
//            {
//                if (std::strcmp(porousMediumBoxType.c_str(), "rectangular") == 0 || dim == 2)
//                {
//                    if (globalPos[0] > darcyXLeftFront()[0]
//                        && globalPos[0] < darcyXRightBack()[0]
//#if DUMUX_MULTIDOMAIN_DIM > 2
//                        && globalPos[2] > darcyXLeftFront()[2]
//                        && globalPos[2] < darcyXRightBack()[2]
//#endif
//                        && (globalPos[1] < interfaceVerticalPos(globalPos)))
//                    {
//                        this->mdGrid().addToSubDomain(darcySubDomainIdx, *eIt);
//                    }
//                    else if ((globalPos[1] > interfaceVerticalPos() && globalPos[0] < darcyXLeftFront()[0])
//                            || (globalPos[1] > interfaceVerticalPos() && globalPos[0] > darcyXRightBack()[0])
//                            || (globalPos[0] > darcyXLeftFront()[0] && globalPos[0] < darcyXRightBack()[0]
//#if DUMUX_MULTIDOMAIN_DIM < 3
//                            )
//#else
//                              && globalPos[2] > darcyXLeftFront()[2] && globalPos[2] < darcyXRightBack()[2])
//                            || (globalPos[1] > interfaceVerticalPos() && globalPos[2] < darcyXLeftFront()[2])
//                            || (globalPos[1] > interfaceVerticalPos() && globalPos[2] > darcyXRightBack()[2])
//#endif
//                            )
//                    {
//                        this->mdGrid().addToSubDomain(stokesSubDomainIdx, *eIt);
//                    }
//                }
//                else if (std::strcmp(porousMediumBoxType.c_str(), "cylindrical") == 0)
//                {
//                    Scalar radius = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, Radius);
//                    GlobalPosition center = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, GlobalPosition, Grid, CircleCenter);
//
//                    using std::pow;
//                    if ((pow(globalPos[0] - center[0], 2) + pow(globalPos[2] - center[2],2)) < pow(radius, 2)
//                        && globalPos[1] < interfaceVerticalPos())
//                    {
//                        this->mdGrid().addToSubDomain(darcySubDomainIdx, *eIt);
//                    }
//                    else if (globalPos[1] > interfaceVerticalPos())
//                    {
//                        this->mdGrid().addToSubDomain(stokesSubDomainIdx, *eIt);
//                    }
//                }
//            }
//        }
//        this->mdGrid().preUpdateSubDomains();
//        this->mdGrid().updateSubDomains();
//        this->mdGrid().postUpdateSubDomains();
//
//        this->darcyElementIndices_.resize(this->sdGridViewDarcy().size(0));
//        Dune::MultipleCodimMultipleGeomTypeMapper<MultiDomainGridView, Dune::MCMGElementLayout>
//            multidomainDofMapper(this->mdGridView());
//        Dune::MultipleCodimMultipleGeomTypeMapper<SubDomainGridView, Dune::MCMGElementLayout>
//            subdomainDofMapper(this->sdGridViewDarcy());
//        for (auto eIt = this->sdGridDarcy().template leafbegin<0>();
//              eIt != this->sdGridDarcy().template leafend<0>(); ++eIt)
//        {
//            this->darcyElementIndices_[subdomainDofMapper.index(*eIt)] =
//                multidomainDofMapper.index(this->mdGrid().multiDomainEntity(*eIt));
//        }
//    }

    // TODO not needed for flat interface profile
//    //! \brief Returns the vertical position of the interface
//    Scalar interfaceVerticalPos(GlobalPosition globalPos)
//    {
//        Scalar amplitude = std::numeric_limits<Scalar>::quiet_NaN();
//        Scalar baseline = std::numeric_limits<Scalar>::quiet_NaN();
//        Scalar offset = std::numeric_limits<Scalar>::quiet_NaN();
//        Scalar scaling = std::numeric_limits<Scalar>::quiet_NaN();
//        std::string interfaceProfile = GET_PARAM_FROM_GROUP(TypeTag, std::string, Grid, InterfaceProfile);
//        if (std::strcmp(interfaceProfile.c_str(), "sinus") == 0
//            || std::strcmp(interfaceProfile.c_str(), "rectangle") == 0
//            || std::strcmp(interfaceProfile.c_str(), "triangle") == 0)
//        {
//            amplitude = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, ObstacleAmplitude);
//            baseline = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, ObstacleBaseline);
//            offset = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, ObstacleOffset);
//            scaling = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, ObstacleScaling);
//        }
//
//        if (std::strcmp(interfaceProfile.c_str(), "sinus") == 0)
//            return amplitude * std::sin((globalPos[0]-offset) / scaling * 2.0 * M_PI) + baseline;
//        else if (std::strcmp(interfaceProfile.c_str(), "rectangle") == 0)
//            return std::sin((globalPos[0]-offset) / scaling * 2.0 * M_PI) > 0 ? baseline+amplitude : baseline-amplitude;
//        else if (std::strcmp(interfaceProfile.c_str(), "triangle") == 0)
//        {
//            GlobalPosition relativePos(globalPos);
//            relativePos[0] = fmod((globalPos[0]-offset) / scaling + 0.25, 0.5) * 2.0;
//            if (relativePos[0] < 0)
//                relativePos[0] = 1.0 + relativePos[0];
//            std::cout << globalPos[0]
//                      << " " << (globalPos[0]-offset)
//                      << " " << (globalPos[0]-offset) / scaling
//                      << " " << relativePos[0]
//                      << " " << (std::cos((globalPos[0]-offset) / scaling * 2.0 * M_PI) > 0)
//                      << std::endl;
//            if (std::cos((globalPos[0]-offset) / scaling * 2.0 * M_PI) > 1e-6)
//                return baseline-amplitude + 2.0 * amplitude * relativePos[0];
//            else
//                return baseline+amplitude - 2.0 * amplitude * relativePos[0];
//        }
//
//        return interfaceVerticalPos();
//    }

//    //! \brief Returns the lower bound of the wind tunnel
//    Scalar interfaceVerticalPos()
//    {
//        return GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, InterfacePosY);
//    }

//    //! \brief Returns the left front position of the darcy domain
//    GlobalPosition darcyXLeftFront()
//    {
//        GlobalPosition darcyXLeftFront(0.0);
//        darcyXLeftFront[0] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, NoDarcyX1);
//#if DUMUX_MULTIDOMAIN_DIM > 2
//        darcyXLeftFront[2] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, NoDarcyZ1);
//#endif
//        return darcyXLeftFront;
//    }

//    //! \brief Returns the right back position of the darcy domain
//    GlobalPosition darcyXRightBack()
//    {
//        GlobalPosition darcyXRightBack(0.0);
//        darcyXRightBack[0] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, NoDarcyX2);
//#if DUMUX_MULTIDOMAIN_DIM > 2
//        darcyXRightBack[2] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, NoDarcyZ2);
//#endif
//        return darcyXRightBack;
//    }

private:
    unsigned freqOutput_;
    Scalar episodeLength_;
};

} //end namespace

#endif // DUMUX_WINDTUNNEL_TEST_PROBLEM_HH
