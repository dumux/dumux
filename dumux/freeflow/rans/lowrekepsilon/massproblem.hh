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
 * k and epsilon-tilde are solved by their transport equations all the way to the wall, with
 * Dirichlet conditions k=0, epsilon-tilde=0 on the wall faces, so no near-wall-region or
 * matching-point bookkeeping and no internal Dirichlet constraints are needed (in contrast to
 * the wall-function model in dumux/freeflow/rans/kepsilon/massproblem.hh). The eddy viscosity is
 * evaluated from the current Newton iterate. \c yPlus(eIdx) is computed from the velocity
 * gradient of the wall-adjacent element, using the wallElementIndex()/velocityGradient()/
 * wallNormalAxis()/flowDirectionAxis() accessors of Dumux::RANSMomentumProblem.
 *
 * \tparam TypeTag The mass sub-model's TypeTag.
 * \tparam MomentumProblem The concrete momentum problem type (a Dumux::LowReKEpsilonMomentumProblem).
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

    //! \note Returns 0.0 as long as the stored density of RANSMomentumProblem is zero, i.e. before
    //! the first call to updateDynamicWallProperties(), to avoid a division by zero (NaN) when
    //! the volume variables are first cached.
    Scalar kinematicViscosity(std::size_t eIdx) const
    {
        const auto density = momentumProblem_->storedDensity(eIdx);
        return density > 0.0 ? momentumProblem_->kinematicViscosity(eIdx) : 0.0;
    }

    //! The wall unit y+ = y*uStar/nu, with uStar computed from the velocity gradient of the
    //! wall-adjacent element.
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
     * \brief The low-Reynolds k-epsilon production/destruction/damping-function source term.
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
