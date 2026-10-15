// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Hyperelastic
 * \brief Variables of the hyperelastic model
 */
#ifndef DUMUX_SOLIDMECHANICS_HYPERELASTIC_VARIABLES_HH
#define DUMUX_SOLIDMECHANICS_HYPERELASTIC_VARIABLES_HH

#include <dumux/common/concepts/ipdata_.hh>

namespace Dumux {

/*!
 * \ingroup Hyperelastic
 * \brief Variables of the hyperelastic model: the displacement at a local dof
 */
template <class Traits>
class HyperelasticVariables
{
    using Scalar = typename Traits::PrimaryVariables::value_type;
public:
    using PrimaryVariables = typename Traits::PrimaryVariables;
    using Indices = typename Traits::ModelTraits::Indices;

    template<class ElementSolution, class Problem, class ElementDiscretization, Concept::LocalDofIpData IpData>
    void update(const ElementSolution& elemSol,
                const Problem& problem,
                const ElementDiscretization& elemDisc,
                const IpData& ipData)
    {
        priVars_ = elemSol[ipData.localDofIndex()];
        extrusionFactor_ = problem.spatialParams().extrusionFactor(elemDisc, ipData, elemSol);
    }

    Scalar extrusionFactor() const
    { return extrusionFactor_; }

    const PrimaryVariables& displacement() const
    { return priVars_; }

    Scalar displacement(const int i) const
    { return priVars_[i]; }

    Scalar priVar(const int pvIdx) const
    { return priVars_[pvIdx]; }

    const PrimaryVariables& priVars() const
    { return priVars_; }

private:
    PrimaryVariables priVars_;
    Scalar extrusionFactor_;
};

} // end namespace Dumux

#endif
