// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PhaseFieldModels
 * \brief Variables for single-equation (scalar) phase-field models, shared
 *        across models whose only primary variable is the phase field itself
 *        (e.g. Allen-Cahn, conservative level set).
 */
#ifndef DUMUX_PHASEFIELD_COMMON_VARIABLES_HH
#define DUMUX_PHASEFIELD_COMMON_VARIABLES_HH

#include <dumux/common/concepts/ipdata_.hh>

namespace Dumux::PhaseField {

/*!
 * \ingroup PhaseFieldModels
 * \brief Variables for single-equation (scalar) phase-field models.
 * \note The update interface takes interpolation point data (rather than a
 *       sub-control volume), since the hybrid CVFE local dofs without an
 *       associated sub-control volume (e.g. the edge dofs of the PQ2/PQ3
 *       schemes) must be updated the same way as the control-volume dofs.
 */
template<class Traits>
class ScalarVariables
{
    using Scalar = typename Traits::PrimaryVariables::value_type;
    static_assert(Traits::PrimaryVariables::dimension == Traits::ModelTraits::numEq());
    static_assert(Traits::ModelTraits::numEq() == 1, "ScalarVariables is for single-equation phase-field models only");

public:
    //! export the type used for the primary variables
    using PrimaryVariables = typename Traits::PrimaryVariables;
    //! export the indices type
    using Indices = typename Traits::ModelTraits::Indices;

    /*!
     * \brief Update all quantities for a given local dof
     *
     * \param elemSol A vector containing all primary variables connected to the element
     * \param problem The object specifying the problem which ought to be simulated
     * \param fvGeometry The local geometry
     * \param ipData The interpolation point data
     */
    template<class ElementSolution, class Problem, class FVElementGeometry, Concept::LocalDofIpData IpData>
    void update(const ElementSolution& elemSol,
                const Problem& problem,
                const FVElementGeometry& fvGeometry,
                const IpData& ipData)
    {
        priVars_ = elemSol[ipData.localDofIndex()];
    }

    Scalar phaseField() const
    { return priVars_[Indices::phaseFieldIdx]; }

    Scalar priVar(const int pvIdx) const
    { return priVars_[pvIdx]; }

    const PrimaryVariables& priVars() const
    { return priVars_; }

    Scalar extrusionFactor() const
    { return 1.0; }

private:
    PrimaryVariables priVars_;
};

} // end namespace Dumux::PhaseField

#endif
