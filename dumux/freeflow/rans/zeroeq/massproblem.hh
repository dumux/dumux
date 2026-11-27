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
 * \brief Mixin that adds the mass-domain side of the zero-equation RANS turbulence closure.
 *        It is only used by the nonisothermal ("...NI") TypeTag
 *        (Dumux::NavierStokesMassOneZeroEqNI); the isothermal model uses
 *        Dumux::NavierStokesMassProblem directly.
 *
 * The zero-equation model computes the eddy viscosity on the momentum side
 * (Dumux::ZeroEqProblem::effectiveViscosity()). For the energy equation to account for an eddy
 * thermal conductivity, the mass domain needs read-only access to it, which is provided by
 * this class.
 */
template<class TypeTag, class MomentumProblem>
class ZeroEqMassProblem : public NavierStokesMassProblem<TypeTag>
{
    using ParentType = NavierStokesMassProblem<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

public:
    using ParentType::ParentType;

    //! Sets the momentum problem. Must be called once, after both sub-problems have been constructed.
    void setMomentumProblem(std::shared_ptr<const MomentumProblem> momentumProblem)
    { momentumProblem_ = momentumProblem; }

    Scalar dynamicEddyViscosity(std::size_t eIdx) const
    { return momentumProblem_->dynamicEddyViscosity(eIdx); }

private:
    std::shared_ptr<const MomentumProblem> momentumProblem_;
};

} // end namespace Dumux

#endif
