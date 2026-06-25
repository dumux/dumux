// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_BOUSSINESQ_VORTICITY_VOLUME_VARIABLES_HH
#define DUMUX_BOUSSINESQ_VORTICITY_VOLUME_VARIABLES_HH

#include <array>

namespace Dumux {

/*!
 * \brief Volume variables for the vector-potential Boussinesq model.
 *
 * Stores nPot vector-potential components A_k and nComp concentrations C_i,
 * plus the local porosity.  No fluid-system interaction — all physical
 * coefficients (K, μ, ρ₀, β) are read from the spatial parameters at flux time.
 *
 * \tparam Traits  must provide PrimaryVariables and ModelTraits (with Indices,
 *                 numPotentialEqs, numFluidComponents()).
 */
template<class Traits>
class BoussinesqVorticityVolumeVariables
{
    using Scalar  = typename Traits::PrimaryVariables::value_type;
    using MT      = typename Traits::ModelTraits;
    using Indices = typename MT::Indices;

    static constexpr int nPot  = MT::numPotentialEqs;
    static constexpr int nComp = MT::numFluidComponents();

public:
    using PrimaryVariables = typename Traits::PrimaryVariables;
    using ModelTraits      = MT;

    static constexpr int numFluidPhases() { return 1; }

    template<class ElemSol, class Problem, class Element, class Scv>
    void update(const ElemSol& elemSol,
                const Problem& problem,
                const Element& element,
                const Scv&     scv)
    {
        priVars_ = elemSol[scv.localDofIndex()];
        for (int k = 0; k < nPot; ++k)
            A_[k] = priVars_[Indices::vectorPotentialIdx(k)];
        for (int i = 0; i < nComp; ++i)
            c_[i] = priVars_[Indices::transportIdx(i)];
        porosity_ = problem.spatialParams().porosity(element, scv, elemSol);
    }

    //! k-th vector-potential component (k=0 is ψ in 2D)
    Scalar vectorPotential(int k = 0) const { return A_[k]; }
    //! Alias kept for 2D single-component backward compat
    Scalar streamfunction()            const { return A_[0]; }

    //! i-th transported concentration
    Scalar concentration(int i = 0)    const { return c_[i]; }

    Scalar porosity()        const { return porosity_; }
    Scalar extrusionFactor() const { return 1.0; }

    //! Full primary variable vector — required by CVFELocalAssembler Dirichlet path
    const PrimaryVariables& priVars() const { return priVars_; }

private:
    PrimaryVariables          priVars_{};
    std::array<Scalar, nPot>  A_{};
    std::array<Scalar, nComp> c_{};
    Scalar porosity_ = 1.0;
};

} // namespace Dumux

#endif
