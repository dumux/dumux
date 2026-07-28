// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FreeflowModels
 * \copydoc Dumux::KOmegaMomentumProblem
 */
#ifndef DUMUX_RANS_KOMEGA_MOMENTUM_PROBLEM_HH
#define DUMUX_RANS_KOMEGA_MOMENTUM_PROBLEM_HH

#include <dumux/common/properties.hh>
#include <dumux/freeflow/rans/common/ransmomentumproblem.hh>

namespace Dumux {

/*!
 * \ingroup FreeflowModels
 * \brief Mixin that adds the momentum-domain side of the two-equation k-omega RANS turbulence
 *        closure to a face-centered-staggered Navier-Stokes momentum problem: the
 *        wall-distance/velocity-gradient/vorticity/stress-tensor bookkeeping inherited from
 *        Dumux::RANSMomentumProblem, plus forwarding the turbulent (eddy) viscosity and the
 *        turbulent kinetic energy k - both computed on the *mass* sub-model, see
 *        dumux/freeflow/rans/komega/volumevariables.hh - into effectiveViscosity() and the
 *        normal-stress correction respectively. Analogous to Dumux::OneEqMomentumProblem/
 *        Dumux::ZeroEqProblem.
 */
template<class TypeTag>
class KOmegaMomentumProblem : public RANSMomentumProblem<TypeTag>
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

    //! Adds the (current-iterate) turbulent viscosity, read from the mass domain's k/omega
    //! through the coupling manager, to the molecular one.
    Scalar effectiveViscosity(const Element& element,
                              const FVElementGeometry& fvGeometry,
                              const SubControlVolumeFace& scvf) const
    { return ParentType::effectiveViscosity(element, fvGeometry, scvf) + this->couplingManager().turbulentViscosity(element, fvGeometry, scvf); }

    //! \copydoc effectiveViscosity(const Element&,const FVElementGeometry&,const SubControlVolumeFace&) const
    Scalar effectiveViscosity(const Element& element,
                              const FVElementGeometry& fvGeometry,
                              const SubControlVolume& scv) const
    { return ParentType::effectiveViscosity(element, fvGeometry, scv) + this->couplingManager().turbulentViscosity(element, fvGeometry, scv); }

    //! The turbulent kinetic energy k, read live from the mass domain through the coupling
    //! manager - detected (duck-typed) by
    //! NavierStokesMomentumFCStaggeredFluxVariables::turbulentKineticEnergyContribution() to
    //! add the isotropic 2/dim*rho*k normal-stress correction to the momentum flux.
    Scalar turbulentKineticEnergy(const Element& element,
                                  const FVElementGeometry& fvGeometry,
                                  const SubControlVolumeFace& scvf) const
    { return this->couplingManager().turbulentKineticEnergy(element, fvGeometry, scvf); }
};

} // end namespace Dumux

#endif
