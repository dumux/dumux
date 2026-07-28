// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FreeflowModels
 * \copydoc Dumux::ZeroEqMassProblem
 */
#ifndef DUMUX_RANS_ZEROEQ_MASS_PROBLEM_HH
#define DUMUX_RANS_ZEROEQ_MASS_PROBLEM_HH

#include <memory>

#include <dumux/common/properties.hh>
#include <dumux/freeflow/navierstokes/mass/problem.hh>

namespace Dumux {

/*!
 * \ingroup FreeflowModels
 * \brief Mixin that adds the mass-domain side of the zero-equation RANS turbulence closure -
 *        used only by the nonisothermal ("...NI") TypeTag (Dumux::NavierStokesMassOneZeroEqNI),
 *        since the isothermal zero-equation test needs no mass-side wiring at all (it uses the
 *        plain Dumux::NavierStokesMassProblem directly, see test/freeflow/rans/zeroeq/channel).
 *
 * Unlike every other RANS model, zero-eq computes and adds its eddy viscosity entirely on the
 * *momentum* side (Dumux::ZeroEqProblem::effectiveViscosity()), with no mass-domain
 * representation and no coupling-manager `turbulentViscosity()` involvement. For the
 * nonisothermal energy equation to pick up an eddy thermal conductivity, the mass domain needs
 * read-only access to that momentum-side eddy viscosity - forwarded here the same way every
 * other RANS model's mass problem forwards *momentum*-owned data (wallDistance,
 * stressTensorScalarProduct, ...), just for a different quantity.
 */
template<class TypeTag, class MomentumProblem>
class ZeroEqMassProblem : public NavierStokesMassProblem<TypeTag>
{
    using ParentType = NavierStokesMassProblem<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

public:
    using ParentType::ParentType;

    //! Must be called once, after both sub-problems have been constructed.
    void setMomentumProblem(std::shared_ptr<const MomentumProblem> momentumProblem)
    { momentumProblem_ = momentumProblem; }

    Scalar dynamicEddyViscosity(std::size_t eIdx) const
    { return momentumProblem_->dynamicEddyViscosity(eIdx); }

private:
    std::shared_ptr<const MomentumProblem> momentumProblem_;
};

} // end namespace Dumux

#endif
