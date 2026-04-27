// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_FREEFLOW_NAVIERSTOKES_MASS_2PVAPOR_VOLUME_VARIABLES_HH
#define DUMUX_FREEFLOW_NAVIERSTOKES_MASS_2PVAPOR_VOLUME_VARIABLES_HH

#include <dumux/freeflow/navierstokes/mass/2p/volumevariables.hh>

namespace Dumux {

/*!
 * \brief Volume variables for the 2p Cahn-Hilliard model extended with vapor transport.
 */
template <class Traits>
class NavierStokesMassTwoPVaporVolumeVariables
: public NavierStokesMassTwoPVolumeVariables<Traits>
{
    using ParentType = NavierStokesMassTwoPVolumeVariables<Traits>;
    using Scalar = typename Traits::PrimaryVariables::value_type;

public:
    using Indices = typename Traits::ModelTraits::Indices;
    using PrimaryVariables = typename Traits::PrimaryVariables;

    template<class ElementSolution, class Problem, class Element, class SubControlVolume>
    void update(const ElementSolution& elemSol,
                const Problem& problem,
                const Element& element,
                const SubControlVolume& scv)
    {
        ParentType::update(elemSol, problem, element, scv);
        vaporConcentration_ = this->priVar(Indices::vaporIdx);
    }

    //! Returns the vapor mass concentration c_v [kg/m³]
    Scalar vaporConcentration() const
    { return vaporConcentration_; }

private:
    Scalar vaporConcentration_;
};

} // end namespace Dumux

#endif
