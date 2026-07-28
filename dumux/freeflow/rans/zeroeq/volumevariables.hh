// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FreeflowModels
 * \copydoc Dumux::ZeroEqMassVolumeVariables
 */
#ifndef DUMUX_RANS_ZEROEQ_MASS_VOLUME_VARIABLES_HH
#define DUMUX_RANS_ZEROEQ_MASS_VOLUME_VARIABLES_HH

#include <dumux/freeflow/navierstokes/mass/1p/volumevariables.hh>
#include <dumux/freeflow/rans/common/thermalconductivitymodel.hh>

namespace Dumux {

/*!
 * \ingroup FreeflowModels
 * \brief Volume variables for the single-phase, nonisothermal Navier-Stokes mass balance, with
 *        a `dynamicEddyViscosity()` forwarded (read-only) from the zero-equation model's
 *        momentum-side eddy viscosity (Dumux::ZeroEqProblem), so the energy equation's eddy
 *        thermal conductivity can be computed - see dumux/freeflow/rans/zeroeq/massproblem.hh.
 *        Used only by the nonisothermal TypeTag (Dumux::NavierStokesMassOneZeroEqNI); the
 *        isothermal zero-equation test needs no mass-side turbulence data at all.
 */
template<class Traits>
class ZeroEqMassVolumeVariables
: public NavierStokesMassOnePVolumeVariables<Traits>
{
    using ParentType = NavierStokesMassOnePVolumeVariables<Traits>;
    using Scalar = typename Traits::PrimaryVariables::value_type;

public:
    template<class ElementSolution, class Problem, class Element, class SubControlVolume>
    void update(const ElementSolution& elemSol,
                const Problem& problem,
                const Element& element,
                const SubControlVolume& scv)
    {
        ParentType::update(elemSol, problem, element, scv);

        dynamicEddyViscosity_ = problem.dynamicEddyViscosity(scv.elementIndex());

        // Nonisothermal only: see dumux/freeflow/rans/common/thermalconductivitymodel.hh's
        // class docs for why this must be done explicitly here rather than automatically.
        if constexpr (Traits::ModelTraits::enableEnergyBalance())
            this->lambdaEff_ = RANSThermalConductivityModel::turbulentEffectiveThermalConductivity(*this);
    }

    Scalar dynamicEddyViscosity() const
    { return dynamicEddyViscosity_; }

private:
    Scalar dynamicEddyViscosity_ = 0.0;
};

} // end namespace Dumux

#endif
