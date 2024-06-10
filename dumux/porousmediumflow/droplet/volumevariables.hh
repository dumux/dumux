// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPModel
 * \brief Contains the quantities which are constant within a
 *        finite volume in the two-phase model.
 */

#ifndef DUMUX_2P_VOLUME_VARIABLES_DROPLET_HH
#define DUMUX_2P_VOLUME_VARIABLES_DROPLET_HH

// #include <dumux/porousmediumflow/volumevariables.hh>
// #include <dumux/porousmediumflow/nonisothermal/volumevariables.hh>
// #include <dumux/material/solidstates/updatesolidvolumefractions.hh>
// #include <dumux/porousmediumflow/2p/formulation.hh>
#include <dumux/porousmediumflow/2p/volumevariables.hh>

namespace Dumux {

/*!
 * \ingroup TwoPModel
 * \brief Contains the quantities which are are constant within a
 *        finite volume in the two-phase model.
 */
template <class Traits>
class TwoPVolumeVariablesDroplet
: public TwoPVolumeVariables <Traits>
{
    using ParentType = TwoPVolumeVariables<Traits>;
    using EnergyVolVars = EnergyVolumeVariables<Traits, TwoPVolumeVariables<Traits> >;
    using PermeabilityType = typename Traits::PermeabilityType;
    using ModelTraits = typename Traits::ModelTraits;
    using Idx = typename ModelTraits::Indices;
    using Scalar = typename Traits::PrimaryVariables::value_type;
    using FS = typename Traits::FluidSystem;
    static constexpr int numFluidComps = ParentType::numFluidComponents();
    enum
    {
        pressureIdx = Idx::pressureIdx,
        saturationIdx = Idx::saturationIdx,

        phase0Idx = FS::phase0Idx,
        phase1Idx = FS::phase1Idx
    };

    static constexpr auto formulation = ModelTraits::priVarFormulation();

public:
    //! Export type of fluid system
    using FluidSystem = typename Traits::FluidSystem;
    //! Export type of fluid state
    using FluidState = typename Traits::FluidState;
    //! Export the indices
    using Indices = typename ModelTraits::Indices;
    //! Export type of solid state
    using SolidState = typename Traits::SolidState;
    //! Export type of solid system
    using SolidSystem = typename Traits::SolidSystem;

    /*!
     * \brief Updates all quantities for a given control volume.
     *
     * \param elemSol A vector containing all primary variables connected to the element
     * \param problem The object specifying the problem which ought to
     *                be simulated
     * \param element An element which contains part of the control volume
     * \param scv The sub control volume
    */
    template<class ElemSol, class Problem, class Element, class Scv>
    void update(const ElemSol &elemSol,
                const Problem &problem,
                const Element &element,
                const Scv& scv)
    {
        ParentType::update(elemSol, problem, element, scv);

        const auto& globalPos = scv.dofPosition();
        if (problem.dropSolver()->isCoupledWithDroplet(globalPos))
        {
            const auto& droplet = problem.dropSolver()->droplet();
            dropletVolume_ = droplet.volume();
            dropletRadius_ = droplet.radius();
            dropletContactRadius_ = droplet.contactRadius();
            dropletContactAngle_ = droplet.contactAngle();
            dropletPc_ = problem.dropSolver()->Pc(droplet);
        }

    }

    Scalar dropletVolume() const
    { return dropletVolume_; }

    Scalar dropletRadius() const
    { return dropletRadius_; }

    Scalar dropletContactRadius() const
    { return dropletContactRadius_; }

    Scalar dropletContactAngle() const
    { return 180 * dropletContactAngle_ / M_PI; }

    Scalar dropletPc() const
    { return dropletPc_; }

private:
    Scalar dropletVolume_ = 0.0;
    Scalar dropletRadius_ = 0.0;
    Scalar dropletContactRadius_ = 0.0;
    Scalar dropletContactAngle_ = 0.0;
    Scalar dropletPc_ = 0.0;

};

} // end namespace Dumux

#endif
