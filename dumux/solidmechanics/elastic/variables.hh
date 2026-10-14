// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Elastic
 * \brief Variables of the elastic model for the assembly with local degrees of freedom
 */
#ifndef DUMUX_SOLIDMECHANICS_ELASTIC_VARIABLES_HH
#define DUMUX_SOLIDMECHANICS_ELASTIC_VARIABLES_HH

#include <dumux/common/concepts/ipdata_.hh>

namespace Dumux {

/*!
 * \ingroup Elastic
 * \brief Variables of the elastic model at a local degree of freedom: the displacement.
 *        Used with `Experimental::GridVariables`; the classic finite-volume assembly uses
 *        `ElasticVolumeVariables`.
 */
template <class Traits>
class ElasticVariables
{
    using Scalar = typename Traits::PrimaryVariables::value_type;
public:
    using PrimaryVariables = typename Traits::PrimaryVariables;
    using DisplacementVector = typename Traits::DisplacementVector;
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

    Scalar displacement(unsigned int dir) const
    { return priVars_[Indices::momentum(dir)]; }

    DisplacementVector displacement() const
    {
        DisplacementVector d;
        for (int dir = 0; dir < d.size(); ++dir)
            d[dir] = displacement(dir);
        return d;
    }

    Scalar priVar(const int pvIdx) const
    { return priVars_[pvIdx]; }

    const PrimaryVariables& priVars() const
    { return priVars_; }

    Scalar extrusionFactor() const
    { return extrusionFactor_; }

private:
    PrimaryVariables priVars_;
    Scalar extrusionFactor_;
};

} // end namespace Dumux

#endif
