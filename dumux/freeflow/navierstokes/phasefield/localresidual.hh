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
 * \ingroup FreeflowNIModel
 * \copydoc Dumux::FreeFlowEnergyLocalResidual
 */
#ifndef DUMUX_NAVIERSTOKES_PHASEFIELD_LOCAL_RESIDUAL_HH
#define DUMUX_NAVIERSTOKES_PHASEFIELD_LOCAL_RESIDUAL_HH

#include <dumux/discretization/method.hh>
#include <dumux/flux/referencesystemformulation.hh>

namespace Dumux {


/*!
 * \ingroup FreeflowNIModel
 * \brief DOCME
 */
template<class ModelTraits>
class NavierStokesPhasefieldLocalResidual
{
    using Indices = typename ModelTraits::Indices;
    static constexpr bool enablePhasefield = ModelTraits::enablePhasefield();
public:

    //! The energy storage in the fluid phase
    template <class NumEqVector, class VolumeVariables>
    static void fluidPhaseStorage(NumEqVector& storage,
                                  const VolumeVariables& volVars)
    {
        if constexpr (enablePhasefield)
            storage[Indices::phasefieldEqIdx] = volVars.density() * volVars.internalEnergy();
            const static Scalar xi = getParam<Scalar>("Phasefield.xi");
            ////const static Scalar sigma = problem.getSigma();
            const static Scalar alpha = getParam<Scalar>("Phasefield.alpha");
            //const static Scalar alpha = problem.getAlpha();
            const static Scalar delta = getParam<Scalar>("Phasefield.delta");
            //const static Scalar rhoD = getParam<Scalar>("Problem.rhoD");
            //const static Scalar rhoP = getParam<Scalar>("Problem.rhoP");
            //const static std::array<Scalar, 3> b_D = { -1.0, 1.0, 0.0 };
            //const static std::array<Scalar, 3> b_P = { 0.0, 1.0, 1.0 };
            //NumEqVector storage{};
            // Phasefields
            storage[Indices::phiIdx]  = xi*xi*alpha * volVars.priVar(Indices::phiIdx);
            //storage[Indices::p2Idx]  = xi*xi*alpha * volVars.priVar(Indices::p2Idx);
            //storage[Indices::p3Idx]  = xi*xi*alpha * volVars.priVar(Indices::p3Idx);
            //storage[Indices::p11Idx] = xi*xi*alpha * volVars.priVar(Indices::p11Idx);
            //storage[Indices::p12Idx] = xi*xi*alpha * volVars.priVar(Indices::p12Idx);
            //// Storage of species in fluid and mineral phases
            //storage[Indices::uAIdx] = (volVars.priVar(Indices::p11Idx) + delta) * volVars.priVar(Indices::uAIdx)
            //        + rhoP * b_P[0] * (1 - volVars.priVars()[Indices::p11Idx] - volVars.priVars()[Indices::p12Idx])
            //        + rhoD * b_D[0] * volVars.priVars()[Indices::p12Idx]
            //        ;
            //storage[Indices::uBIdx] = (volVars.priVar(Indices::p11Idx) + delta) * volVars.priVar(Indices::uBIdx)
            //        + rhoP * b_P[1] * (1 - volVars.priVars()[Indices::p11Idx] - volVars.priVars()[Indices::p12Idx])
            //        + rhoD * b_D[1] * volVars.priVars()[Indices::p12Idx]
            //        ;
            //storage[Indices::uCIdx] = (volVars.priVar(Indices::p11Idx) + delta) * volVars.priVar(Indices::uCIdx)
            //        + rhoP * b_P[2] * (1 - volVars.priVars()[Indices::p11Idx] - volVars.priVars()[Indices::p12Idx])
            //        + rhoD * b_D[2] * volVars.priVars()[Indices::p12Idx]
            //        ;
            //storage[Indices::u3AIdx] = (volVars.priVar(Indices::p1Idx) + delta) *
            //    volVars.priVar(Indices::u3AIdx)
            //        + rhoP * b_P[0] * volVars.priVars()[Indices::p3Idx]
            //        + rhoD * b_D[0] * volVars.priVars()[Indices::p2Idx]
            //        ;
            //storage[Indices::u3BIdx] = (volVars.priVar(Indices::p1Idx) + delta) *
            //    volVars.priVar(Indices::u3BIdx)
            //        + rhoP * b_P[1] * volVars.priVars()[Indices::p3Idx]
            //        + rhoD * b_D[1] * volVars.priVars()[Indices::p2Idx]
            //        ;
            //storage[Indices::u3CIdx] = (volVars.priVar(Indices::p1Idx) + delta) *
            //    volVars.priVar(Indices::u3CIdx)
            //        + rhoP * b_P[2] * volVars.priVars()[Indices::p3Idx]
            //        + rhoD * b_D[2] * volVars.priVars()[Indices::p2Idx]
            //        ;
    }

    //! brief The advective phase energy flux
    template <class NumEqVector, class FluxVariables>
    static void phasefieldFlux(NumEqVector& flux,
                                   FluxVariables& fluxVars)
    {
        if constexpr (enablePhasefield)
        {
            const auto& scvf = this->scvFace();
            const auto& fvGeometry = this->fvGeometry();
            const auto& elemVolVars = this->elemVolVars();
            const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
            const auto& insideVolVars = elemVolVars[insideScv];
            const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];

            const static Scalar sigma = getParam<Scalar>("Phasefield.sigma");
            //const static Scalar sigma = problem.getSigma();
            const static Scalar xi = getParam<Scalar>("Phasefield.xi");
            const static Scalar delta = getParam<Scalar>("Phasefield.delta");
            //const static Scalar D_u = getParam<Scalar>("Problem.DiffCoeff");
            const static std::array<Scalar, numComponents> diffCoeff = {
                xi*xi*sigma//, xi*xi*sigma, xi*xi*sigma,// p1 - p3
            //    xi*xi*sigma, xi*xi*sigma,// p11, p12
            //    (insideVolVars.priVar(Indices::p11Idx)+delta) * D_u,//u_A
            //    (insideVolVars.priVar(Indices::p11Idx)+delta) * D_u,
            //    (insideVolVars.priVar(Indices::p11Idx)+delta) * D_u,
            //    (insideVolVars.priVar(Indices::p1Idx)+delta) * D_u,//u_3A
            //    (insideVolVars.priVar(Indices::p1Idx)+delta) * D_u,
            //    (insideVolVars.priVar(Indices::p1Idx)+delta) * D_u
            };

            const int numPhasefields = 1;
            for (int k = Indices::phiIdx; k < Indices::phiIdx + numPhasefields; ++k)
            {
                const auto valInside = insideVolVars.priVar(k);
                const auto valOutside = outsideVolVars.priVar(k);

                Scalar tij;

            //    const Scalar ti = computeTpfaTransmissibility(scvf, insideScv, diffCoeff[k],
            //                                                  insideVolVars.extrusionFactor());
                auto distanceVector = scvf.ipGlobal();
                distanceVector -= insideScv.center();
                distanceVector /= distanceVector.two_norm2();
                const Scalar ti = insideVolVars.extrusionFactor() * (distanceVector *
                        scvf.unitOuterNormal());

            //    if (scvf.boundary())
            //        tij = scvf.area()*ti;
            //    else
                {
                    const auto outsideScvIdx = scvf.outsideScvIdx();
                    const auto& outsideScv = fvGeometry.scv(outsideScvIdx);
            //        const Scalar tj = -computeTpfaTransmissibility(scvf, outsideScv, diffCoeff[k],
            //                                                       outsideVolVars.extrusionFactor());
                    distanceVector = scvf.ipGlobal();
                    distanceVector -= outsideScv.center();
                    distanceVector /= distanceVector.two_norm2();
                    const Scalar ti = -outsideVolVars.extrusionFactor() * (distanceVector *
                            scvf.unitOuterNormal());

                    if (ti*tj <= 0.0)
                        tij = 0.0;
                    else
                        tij = scvf.area()*(ti * tj)/(ti + tj);
                }

                flux[k] = tij*(valInside - valOutside);
            //}
            //    auto upwindTerm = [](const auto& volVars)
            //    { return volVars.density()*volVars.enthalpy(); };

            //    flux[Indices::energyEqIdx] += fluxVars.advectiveFlux(upwindTerm);
            //}
    }

    NumEqVector source( const Problem& problem,
                        const Element& element,
                        const FVGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const SubControlVolum& scv) const
    {
        NumEqVector source;
        // extract priVars
        const auto& priVars = elemVolVars[scv].priVars();
        // for ODE system activities f are updated globally when u is. Update both.
        // f and g for 2-phasefields formulation
       // Scalar f_D1 = 0.0;
       // Scalar f_P1 = 0.0;
       // //Scalar g_01 = 0.0;
       // if (priVars[uAIdx] > ua_threshold_)
       //     f_D1 = priVars[uBIdx] / priVars[uAIdx] - 1;
       // else
       //     f_D1 = priVars[uBIdx] / ua_threshold_ - 1;
       // f_P1 = priVars[uBIdx] * priVars[uCIdx] - 1;
       // Scalar f_D = 0.0;
       // Scalar f_P = 0.0;
       // if (priVars[u3AIdx] > ua_threshold_)
       //     f_D = priVars[u3BIdx] / priVars[u3AIdx] - 1;
       // else
       //     f_D = priVars[u3BIdx] / ua_threshold_ - 1;
       // f_P = priVars[u3BIdx] * priVars[u3CIdx] - 1;
       // f_D = f_P = f_D1 = f_P1 = 0.0;

        const Scalar sigma = getParams<Scalar> ("Phasefield.sigma");
        // define sources as in eqs (7-9) and (10-11)
        source[Indices::phiIdx] = 16.0 * sigma * priVars[Indices::phiIdx] *
            (1-priVars[Indices::phiIdx]) * ((1-2*priVars[Indices::phiIdx]);
    }
};


} // end namespace Dumux

#endif
