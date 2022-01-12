// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup TwoEqModel
 * \copydoc Dumux::TwoEqFluxVariablesImpl
 */
#ifndef DUMUX_CC_TWOEQ_FLUXVARIABLES_HH
#define DUMUX_CC_TWOEQ_FLUXVARIABLES_HH

#include <numeric>

#include <dumux/common/properties.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/flux/upwindscheme.hh>

namespace Dumux {

/*!
  * \ingroup TwoEqModel
  * \brief The flux variables class for the k-omega model using the staggered grid discretization.
 */

// forward declaration
template<class TypeTag, class BaseFluxVariables, class DiscretizationMethod>
class TwoEqFluxVariablesImpl;

template<class TypeTag, class BaseFluxVariables>
class TwoEqFluxVariablesImpl<TypeTag, BaseFluxVariables, DiscretizationMethods::CCTpfa>
: public BaseFluxVariables
{
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using UpwindScheme = Dumux::UpwindScheme<GridGeometry>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Extrusion = Extrusion_t<GridGeometry>;

    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using ElementFluxVariablesCache = typename GridVariables::GridFluxVariablesCache::LocalView;
    using ElementVolumeVariables = typename GridVariables::GridVolumeVariables::LocalView;;
    using NumEqVector = Dumux::NumEqVector<GetPropType<TypeTag, Properties::PrimaryVariables>>;

public:

    struct AdvectionType
    {
        static Scalar flux(const Problem& problem,
                           const Element& element,
                           const FVElementGeometry& fvGeometry,
                           const ElementVolumeVariables& elemVolVars,
                           const SubControlVolumeFace& scvf,
                           const int phaseIdx,
                           const ElementFluxVariablesCache& elemFluxVarsCache)
        {
            const auto velocity = problem.faceVelocity(element, fvGeometry, scvf);
            const Scalar volumeFlux = velocity*scvf.unitOuterNormal()*Extrusion::area(scvf)*extrusionFactor_(elemVolVars, scvf);
            return volumeFlux;
        }
    };

    /*!
     * \brief Returns the advective mass flux in kg/s
     *        or the advective mole flux in mole/s.
     */
    NumEqVector advectiveFlux(NumEqVector flux, int phaseIdx = 0) const
    {
        const auto upwind = [this](const auto& upwindTerm)
        {
            const auto& scvf = this->scvFace();
            const auto velocity = this->problem().faceVelocity(this->element(), this->fvGeometry(), scvf);
            const Scalar volumeFlux = velocity*scvf.unitOuterNormal()*Extrusion::area(scvf)*extrusionFactor_(this->elemVolVars(), scvf);
            return UpwindScheme::apply(*this, upwindTerm, volumeFlux, 0/*phaseIdx*/);
        };

        // Prepare upwinding terms:
        auto upwindTermMass = [](const auto& volVars)
        { return volVars.density(); };
        auto upwindTermK = [](const auto& volVars)
        { return volVars.turbulentKineticEnergy()*volVars.density(); };
        auto upwindTermOmega = [](const auto& volVars)
        { return volVars.dissipation()*volVars.density(); };

        flux[ModelTraits::Indices::conti0EqIdx] += upwind(upwindTermMass);
        flux[ModelTraits::Indices::turbulentKineticEnergyEqIdx] += upwind(upwindTermK);
        flux[ModelTraits::Indices::dissipationEqIdx] += upwind(upwindTermOmega);

        return flux;
    }

    NumEqVector diffusiveFlux(NumEqVector flux, int phaseIdx = 0) const
    {
        const auto& scvf = this->scvFace();
        const auto& fvGeometry = this->fvGeometry();
        const auto& elemVolVars = this->elemVolVars();

        // calculate diffusive flux
        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& outsideScv = fvGeometry.scv(scvf.outsideScvIdx());
        const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
        const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];

        // effective diffusion coefficients
        Scalar insideCoeff_k = (insideVolVars.viscosity() + (insideVolVars.sigmaK() * insideVolVars.dynamicEddyViscosity()))
                              * insideVolVars.extrusionFactor();
        Scalar outsideCoeff_k = (outsideVolVars.viscosity() + (outsideVolVars.sigmaK() * outsideVolVars.dynamicEddyViscosity())) // TODO: Call these coeffs correctly
                              * outsideVolVars.extrusionFactor();
        Scalar insideCoeff_w = (insideVolVars.viscosity() + (insideVolVars.sigmaOmega() * insideVolVars.dynamicEddyViscosity()))
                              * insideVolVars.extrusionFactor();
        Scalar outsideCoeff_w = (outsideVolVars.viscosity() + (outsideVolVars.sigmaOmega() * outsideVolVars.dynamicEddyViscosity()))
                              * outsideVolVars.extrusionFactor();
        // average and distance
        Scalar coeff_k = arithmeticMean(insideCoeff_k, outsideCoeff_k,
                                (outsideScv.dofPosition() - scvf.ipGlobal()).two_norm(),
                                (insideScv.dofPosition() - scvf.ipGlobal()).two_norm());
        Scalar coeff_w = arithmeticMean(insideCoeff_w, outsideCoeff_w,
                                (outsideScv.dofPosition() - scvf.ipGlobal()).two_norm(),
                                (insideScv.dofPosition() - scvf.ipGlobal()).two_norm());
        Scalar distance = (outsideScv.dofPosition() - insideScv.dofPosition()).two_norm();

        flux[ModelTraits::Indices::turbulentKineticEnergyEqIdx] += coeff_k / distance
                * (insideVolVars.turbulentKineticEnergy() - outsideVolVars.turbulentKineticEnergy()) * Extrusion::area(scvf);
        flux[ModelTraits::Indices::dissipationEqIdx] += coeff_w / distance
                * (insideVolVars.dissipation() - outsideVolVars.dissipation()) * Extrusion::area(scvf);

        return flux;
    }

    /*!
     * \brief Computes the flux for the cell center residual.
     */
    NumEqVector flux(int phaseIdx = 0) const
    {
        NumEqVector flux(0.0);

//        if constexpr (enableEnergyBalance)
//            add energyFluxes;

//        if constexpr (enableCompositionalBalance)
//            add compositionalFluxes;

        // Advective Flux
        flux += advectiveFlux(flux);
        // Diffusive Flux
        flux += diffusiveFlux(flux);

        return flux;
    }

// TODO: How to incorporate the weird dilitation term to the momentum flux?

private:

    static Scalar extrusionFactor_(const ElementVolumeVariables& elemVolVars, const SubControlVolumeFace& scvf)
    {
        const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
        const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];
        return harmonicMean(insideVolVars.extrusionFactor(), outsideVolVars.extrusionFactor());
    }
};

} // end namespace

#endif
