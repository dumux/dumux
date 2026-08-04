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
 *        problem: the wall-distance/velocity-gradient/vorticity bookkeeping inherited from
 *        Dumux::RANSMomentumProblem, plus forwarding the turbulent (eddy) viscosity - which
 *        is computed from the working viscosity ν̃ solved as an extra transport equation on
 *        the *mass* sub-model, see dumux/freeflow/rans/oneeq/volumevariables.hh - into
 *        effectiveViscosity(), analogous to how Dumux::ZeroEqProblem adds its own (purely
 *        algebraic, not transport-equation-based) eddy viscosity.
 *
 * Unlike the zero-equation model, the eddy viscosity here is read live (i.e. from the
 * *current* Newton iterate, not lagged) through the coupling manager's turbulentViscosity().
 * Only the production-term inputs (velocity gradients/vorticity, handled by the
 * RANSMomentumProblem base) are lagged; this feedback link is not.
 *
 * Usage: derive your test problem from this class instead of directly from
 * Dumux::NavierStokesMomentumProblem<TypeTag>/Dumux::RANSMomentumProblem<TypeTag>, exactly
 * like Dumux::ZeroEqProblem.
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

    //! Adds the (current-iterate) turbulent viscosity, read from the mass domain's ν̃-derived
    //! eddy viscosity through the coupling manager, to the molecular one.
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
