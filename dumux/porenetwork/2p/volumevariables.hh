// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PNMTwoPModel
 * \copydoc Dumux::PoreNetwork::TwoPVolumeVariables
 */

#ifndef DUMUX_PNM_2P_VOLUME_VARIABLES_HH
#define DUMUX_PNM_2P_VOLUME_VARIABLES_HH

#include <dumux/common/parameters.hh>
#include <dumux/porousmediumflow/2p/volumevariables.hh>

namespace Dumux::PoreNetwork {

/*!
 * \ingroup PNMTwoPModel
 * \brief Contains the quantities which are are constant within a
 *        finite volume (the pore body) in the two-phase model.
 */
template <class Traits>
class TwoPVolumeVariables
: public Dumux::TwoPVolumeVariables<Traits>
{
    using ParentType = Dumux::TwoPVolumeVariables<Traits>;
    using ModelTraits = typename Traits::ModelTraits;
    using Scalar = typename Traits::PrimaryVariables::value_type;
    using FS = typename Traits::FluidSystem;
    static constexpr int numFluidComps = ParentType::numFluidComponents();
    static constexpr auto formulation = ModelTraits::priVarFormulation();
public:

    //! Export type of fluid system
    using FluidSystem = typename Traits::FluidSystem;
    //! Export type of fluid state
    using FluidState = typename Traits::FluidState;
    //! Export type of solid state
    using SolidState = typename Traits::SolidState;
    //! Export type of solid system
    using SolidSystem = typename Traits::SolidSystem;
    //! Export the indices
    using Indices = typename ModelTraits::Indices;

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
    void update(const ElemSol& elemSol,
                const Problem& problem,
                const Element& element,
                const Scv& scv)
    {
        ParentType::update(elemSol, problem, element, scv);
        poreInscribedRadius_ = problem.spatialParams().poreInscribedRadius(element, scv, elemSol);
        poreVolume_ = problem.gridGeometry().poreVolume(scv.dofIndex()) * this->porosity();
        surfaceTension_ = problem.spatialParams().surfaceTension(element, scv, elemSol);
    }

    /*!
     * \brief Returns the pore's inscribed radius.
     */
    Scalar poreInscribedRadius() const
    { return poreInscribedRadius_; }

    /*!
     * \brief Returns the pore volume. // TODO should this be a fraction only?
     */
    Scalar poreVolume() const
    { return poreVolume_; }

    /*!
     * \brief Returns the surface tension.
     */
    Scalar surfaceTension() const
    { return surfaceTension_; }

protected:
    Scalar poreInscribedRadius_;
    Scalar poreVolume_;
    Scalar surfaceTension_;
};

} // end namespace Dumux::PoreNetwork

#endif
