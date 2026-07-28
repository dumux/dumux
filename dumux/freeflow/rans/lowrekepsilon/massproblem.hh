// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FreeflowModels
 * \copydoc Dumux::LowReKEpsilonMassProblem
 */
#ifndef DUMUX_RANS_LOWREKEPSILON_MASS_PROBLEM_HH
#define DUMUX_RANS_LOWREKEPSILON_MASS_PROBLEM_HH

#include <memory>

#include <dumux/common/properties.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/freeflow/navierstokes/mass/problem.hh>

namespace Dumux {

/*!
 * \ingroup FreeflowModels
 * \brief Mixin that adds the mass-domain side of the two-equation, low-Reynolds-number k-epsilon
 *        (Chien 1982) RANS turbulence closure: the production/destruction/damping-function
 *        source term (through the standard problem.source() extension point) and the
 *        wallDistance()/stressTensorScalarProduct()/kinematicViscosity()/yPlus() forwarding
 *        used by Dumux::LowReKEpsilonMassVolumeVariables.
 *
 * Unlike the wall-function k-epsilon model (dumux/freeflow/rans/kepsilon/massproblem.hh), this
 * model needs no near-wall-region/matching-point bookkeeping, no internal Dirichlet constraints,
 * and no per-time-step mass-side update at all: k and epsilon-tilde are solved by their ordinary
 * transport PDE all the way to the wall, with a plain (boundary-face) Dirichlet condition k=0,
 * epsilon-tilde=0 there (see the test problem), and the eddy viscosity is always evaluated live
 * from the current Newton iterate (releases/3.10's optional RANS.UseStoredEddyViscosity lagging
 * opt-out is not ported, the same simplification already made for the one-equation and k-omega
 * models). \c yPlus(eIdx) is computed the same way k-epsilon's *local* (non-matching-point) y+
 * already is - from the wall-adjacent element's velocity gradient, via
 * Dumux::RANSMomentumProblem's existing wallElementIndex()/velocityGradient()/
 * wallNormalAxis()/flowDirectionAxis() accessors - no new momentum-side bookkeeping is needed.
 *
 * \tparam TypeTag The mass sub-model's TypeTag.
 * \tparam MomentumProblem The concrete momentum problem type (a Dumux::LowReKEpsilonMomentumProblem),
 *         explicit template parameter, exactly as for Dumux::KOmegaMassProblem.
 *
 * Ported from the deleted releases/3.10:dumux/freeflow/rans/twoeq/lowrekepsilon/problem.hh
 * (RANSProblemImpl<TypeTag, TurbulenceModel::lowrekepsilon>) and, for the wall condition,
 * releases/3.10:test/freeflow/rans/problem.hh (which has no lowrekepsilon-specific branch at
 * all - plain face Dirichlet BCs are used, confirmed this session).
 */
template<class TypeTag, class MomentumProblem>
class LowReKEpsilonMassProblem : public NavierStokesMassProblem<TypeTag>
{
    using ParentType = NavierStokesMassProblem<TypeTag>;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename ModelTraits::Indices;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;

public:
    using ParentType::ParentType;

    //! Must be called once, after both sub-problems have been constructed.
    void setMomentumProblem(std::shared_ptr<const MomentumProblem> momentumProblem)
    { momentumProblem_ = momentumProblem; }

    Scalar wallDistance(std::size_t eIdx) const
    { return momentumProblem_->wallDistance(eIdx); }

    Scalar stressTensorScalarProduct(std::size_t eIdx) const
    { return momentumProblem_->stressTensorScalarProduct(eIdx); }

    //! \note Guarded against the transient 0.0/0.0 that RANSMomentumProblem's stored molecular
    //! density/viscosity briefly evaluate to between updateStaticWallProperties() (which only
    //! zero-initializes them) and the first call to updateDynamicWallProperties() (which
    //! populates them from a real solution) - this window is exercised once, when
    //! GridVariables::init() eagerly caches the t=0 volume variables in main.cc, before the
    //! time loop's first updateDynamicWallProperties() call. Returning 0.0 there (rather than
    //! NaN) keeps that transient state numerically harmless; it is never read by any solved
    //! time step.
    Scalar kinematicViscosity(std::size_t eIdx) const
    {
        const auto density = momentumProblem_->storedDensity(eIdx);
        return density > 0.0 ? momentumProblem_->kinematicViscosity(eIdx) : 0.0;
    }

    //! The wall unit y+ = y*uStar/nu, with uStar computed from the wall-adjacent element's
    //! velocity gradient - the same local (not matching-point-based) formula used by the
    //! wall-function k-epsilon model's yPlus(), reusing the momentum problem's already-existing
    //! flat-wall-bounded bookkeeping.
    Scalar yPlus(std::size_t eIdx) const
    {
        using std::abs;
        using std::sqrt;
        const auto wallElementIdx = momentumProblem_->wallElementIndex(eIdx);
        const auto nu = kinematicViscosity(eIdx);
        if (nu == 0.0)
            return 0.0;

        const auto wallNormalAxis = momentumProblem_->wallNormalAxis();
        const auto flowDirectionAxis = momentumProblem_->flowDirectionAxis();
        const auto uStar = sqrt(kinematicViscosity(wallElementIdx)
                                 *abs(momentumProblem_->velocityGradient(wallElementIdx, flowDirectionAxis, wallNormalAxis)));
        return wallDistance(eIdx)*uStar/nu;
    }

    /*!
     * \brief The low-Reynolds k-epsilon production/destruction/damping-function source term,
     *        ported from releases/3.10:dumux/freeflow/rans/twoeq/lowrekepsilon/staggered/
     *        localresidual.hh's computeSourceForCellCenter(). Added through the standard
     *        problem.source(...) extension point.
     */
    template<class ElementVolumeVariables>
    NumEqVector source(const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolume& scv) const
    {
        NumEqVector source(0.0);
        const auto& volVars = elemVolVars[scv];

        const Scalar productionTerm = 2.0*volVars.dynamicEddyViscosity()*volVars.stressTensorScalarProduct();
        source[Indices::turbulentKineticEnergyEqIdx] += productionTerm;
        source[Indices::dissipationEqIdx] += volVars.cOneEpsilon()*volVars.fOne()
                                            *volVars.dissipationTilde()/volVars.turbulentKineticEnergy()
                                            *productionTerm;

        source[Indices::turbulentKineticEnergyEqIdx] -= volVars.dissipationTilde()*volVars.density();
        source[Indices::dissipationEqIdx] -= volVars.cTwoEpsilon()*volVars.fTwo()*volVars.density()
                                            *volVars.dissipationTilde()*volVars.dissipationTilde()
                                            /volVars.turbulentKineticEnergy();

        source[Indices::turbulentKineticEnergyEqIdx] -= volVars.dValue()*volVars.density();
        source[Indices::dissipationEqIdx] += volVars.eValue()*volVars.density();

        return source;
    }

private:
    std::shared_ptr<const MomentumProblem> momentumProblem_;
};

} // end namespace Dumux

#endif
