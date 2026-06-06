// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoroElasticLargeDef
 * \brief Volume variables for the large-deformation poroelastic mixed u–p_s–p_f model.
 */
#ifndef DUMUX_POROMECHANICS_POROELASTIC_LARGE_DEFORMATIONS_VARIABLES_HH
#define DUMUX_POROMECHANICS_POROELASTIC_LARGE_DEFORMATIONS_VARIABLES_HH

#include <dumux/common/concepts/ipdata_.hh>

namespace Dumux {

/*!
 * \ingroup PoroElasticLargeDef
 * \brief Volume variables for the momentum subdomain (displacement).
 */
template <class Traits>
class PoroElasticLargeDefMomentumVariables
{
    using Scalar = typename Traits::PrimaryVariables::value_type;
    static_assert(Traits::PrimaryVariables::dimension == Traits::ModelTraits::numEq());

public:
    using PrimaryVariables = typename Traits::PrimaryVariables;
    using Indices = typename Traits::ModelTraits::Indices;

    template<class ElementSolution, class Problem, class FVElementGeometry, Concept::LocalDofIpData IpData>
    void update(const ElementSolution& elemSol, const Problem&,
                const FVElementGeometry&, const IpData& ipData)
    { priVars_ = elemSol[ipData.localDofIndex()]; }

    PrimaryVariables displacement() const { return priVars_; }
    Scalar displacement(int dirIdx) const { return priVars_[dirIdx]; }
    Scalar priVar(int pvIdx) const { return priVars_[pvIdx]; }
    const PrimaryVariables& priVars() const { return priVars_; }

private:
    PrimaryVariables priVars_;
};

/*!
 * \ingroup PoroElasticLargeDef
 * \brief Volume variables for the solid pressure subdomain.
 */
template <class Traits>
class PoroElasticLargeDefSolidPressureVariables
{
    using Scalar = typename Traits::PrimaryVariables::value_type;
    static_assert(Traits::PrimaryVariables::dimension == Traits::ModelTraits::numEq());

public:
    using PrimaryVariables = typename Traits::PrimaryVariables;
    using Indices = typename Traits::ModelTraits::Indices;

    template<class ElementSolution, class Problem, class FVElementGeometry, Concept::LocalDofIpData IpData>
    void update(const ElementSolution& elemSol, const Problem&,
                const FVElementGeometry&, const IpData& ipData)
    { priVars_ = elemSol[ipData.localDofIndex()]; }

    Scalar solidPressure() const { return priVars_[0]; }
    Scalar priVar(int pvIdx) const { return priVars_[pvIdx]; }
    const PrimaryVariables& priVars() const { return priVars_; }

private:
    PrimaryVariables priVars_;
};

/*!
 * \ingroup PoroElasticLargeDef
 * \brief Volume variables for the fluid pressure subdomain.
 */
template <class Traits>
class PoroElasticLargeDefFluidPressureVariables
{
    using Scalar = typename Traits::PrimaryVariables::value_type;
    static_assert(Traits::PrimaryVariables::dimension == Traits::ModelTraits::numEq());

public:
    using PrimaryVariables = typename Traits::PrimaryVariables;
    using Indices = typename Traits::ModelTraits::Indices;

    template<class ElementSolution, class Problem, class FVElementGeometry, Concept::LocalDofIpData IpData>
    void update(const ElementSolution& elemSol, const Problem&,
                const FVElementGeometry&, const IpData& ipData)
    { priVars_ = elemSol[ipData.localDofIndex()]; }

    Scalar fluidPressure() const { return priVars_[0]; }
    Scalar priVar(int pvIdx) const { return priVars_[pvIdx]; }
    const PrimaryVariables& priVars() const { return priVars_; }

private:
    PrimaryVariables priVars_;
};

} // end namespace Dumux

#endif
