// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FreeflowModels
 * \copydoc Dumux::OneEqMomentumProblem
 */
#ifndef DUMUX_RANS_ONEEQ_MOMENTUM_PROBLEM_HH
#define DUMUX_RANS_ONEEQ_MOMENTUM_PROBLEM_HH

#include <dumux/common/properties.hh>
#include <dumux/freeflow/rans/common/ransmomentumproblem.hh>

namespace Dumux {

/*!
 * \ingroup FreeflowModels
 * \brief Mixin that adds the momentum-domain side of the one-equation (Spalart-Allmaras)
 *        RANS turbulence closure to a face-centered-staggered Navier-Stokes momentum
 *        problem: the wall-distance, velocity-gradient and vorticity data inherited from
 *        Dumux::RANSMomentumProblem, plus the turbulent (eddy) viscosity in
 *        effectiveViscosity(). The eddy viscosity is computed from the working viscosity ν̃,
 *        which is solved as an extra transport equation on the mass sub-model, see
 *        dumux/freeflow/rans/oneeq/volumevariables.hh.
 *
 * The eddy viscosity is read from the current Newton iterate through the coupling manager's
 * turbulentViscosity(). Only the production-term inputs (velocity gradients and vorticity,
 * provided by RANSMomentumProblem) are lagged.
 *
 * Usage: derive the problem from this class instead of directly from
 * Dumux::NavierStokesMomentumProblem<TypeTag> or Dumux::RANSMomentumProblem<TypeTag>.
 */
template<class TypeTag>
class OneEqMomentumProblem : public RANSMomentumProblem<TypeTag>
{
    using ParentType = RANSMomentumProblem<TypeTag>;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

public:
    using ParentType::ParentType;

    //! Adds the turbulent viscosity (derived from ν̃ on the mass domain and obtained
    //! through the coupling manager) to the molecular one.
    Scalar effectiveViscosity(const Element& element,
                              const FVElementGeometry& fvGeometry,
                              const SubControlVolumeFace& scvf) const
    { return ParentType::effectiveViscosity(element, fvGeometry, scvf) + this->couplingManager().turbulentViscosity(element, fvGeometry, scvf); }

    //! \copydoc effectiveViscosity(const Element&,const FVElementGeometry&,const SubControlVolumeFace&) const
    Scalar effectiveViscosity(const Element& element,
                              const FVElementGeometry& fvGeometry,
                              const SubControlVolume& scv) const
    { return ParentType::effectiveViscosity(element, fvGeometry, scv) + this->couplingManager().turbulentViscosity(element, fvGeometry, scv); }
};

} // end namespace Dumux

#endif
