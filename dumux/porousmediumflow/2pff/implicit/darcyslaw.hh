// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
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
 * \brief This file contains the data which is required to calculate
 *        volume and mass fluxes of fluid phases over a face of a finite volume by means
 *        of the Darcy approximation. Specializations are provided for the different discretization methods.
 */
#ifndef DUMUX_DISCRETIZATION_CC_TPFA_DARCYS_LAW_FRACFLOW_HH
#define DUMUX_DISCRETIZATION_CC_TPFA_DARCYS_LAW_FRACFLOW_HH

#include <memory>

#include <dune/common/float_cmp.hh>

#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>
#include <math.h>

#include <dumux/discretization/cellcentered/tpfa/darcyslaw.hh>

namespace Dumux {

/*!
 * \ingroup DarcysLaw
 * \brief Specialization of Darcy's Law for the CCTpfa method.
 */
template <class ScalarType, class FVGridGeometry>
class FractionalFlowCCTpfaDarcysLaw
: public CCTpfaDarcysLaw<ScalarType, FVGridGeometry, false>
{
    using ThisType = FractionalFlowCCTpfaDarcysLaw<ScalarType, FVGridGeometry>;
    using GridView = typename FVGridGeometry::GridView;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolumeFace = typename FVGridGeometry::SubControlVolumeFace;
    using Element = typename GridView::template Codim<0>::Entity;

    enum {
        viscousFluxIdx = 0,
        gravityFluxIdx = 1,
        capillaryFluxIdx = 2,
    };

    enum {
        wPhaseIdx = 0,
        nPhaseIdx = 1
    };

    enum {
        transportEqIdx = 0
    };

public:
    //! state the scalar type of the law
    using Scalar = ScalarType;

    //! state the discretization method this implementation belongs to
    static const DiscretizationMethod discMethod = DiscretizationMethod::cctpfa;

    //! state the type for the corresponding cache
    using Cache = TpfaDarcysLawCache<ThisType, FVGridGeometry>;

    // the return type of the flux method
    using ReturnType = std::array<Scalar, 3>;

    template<class Problem, class ElementVolumeVariables, class ElementFluxVarsCache>
    static ReturnType flux(const Problem& problem,
                           const Element& element,
                           const FVElementGeometry& fvGeometry,
                           const ElementVolumeVariables& elemVolVars,
                           const SubControlVolumeFace& scvf,
                           int eqIdx,
                           const ElementFluxVarsCache& elemFluxVarsCache)
    {
        ReturnType fluxes;

        const auto& fluxVarsCache = elemFluxVarsCache[scvf];

        // Get the inside and outside volume variables
        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& insideVolVars = elemVolVars[insideScv];
        const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];

        if (eqIdx == transportEqIdx)
        {
            //////////////////////////////////////////////////////////////
            // we return three separate fluxes: the viscous, gravity and capillary flux
            //////////////////////////////////////////////////////////////

            // Get the total velocity from the input file
            using VelocityVector = Dune::FieldVector<Scalar, GridView::dimensionworld>;
//             static const VelocityVector v_t = getParamFromGroup<VelocityVector>(problem.paramGroup(), "Problem.TotalVelocity");

            if (eqIdx == transportEqIdx)
            {
                // viscous Flu
            //TODO: calculate the velocity (no more constant)


             const auto pnInside = insideVolVars.pressure(nPhaseIdx);
             const auto pnOutside = outsideVolVars.pressure(nPhaseIdx);


             const auto rho = [&](int phaseIdx)
                {
                    // boundaries
                    if (scvf.boundary())
                        return insideVolVars.density(phaseIdx);

                    // inner faces with two neighboring elements
                    else
                        return (insideVolVars.density(phaseIdx) + outsideVolVars.density(phaseIdx))*0.5;
                };

                // ask for the gravitational acceleration in the inside neighbor
                const auto xInside = insideScv.center();
                const auto gInside = problem.gravityAtPos(xInside);
                const auto xOutside = scvf.boundary() ? scvf.ipGlobal()
                                                      : fvGeometry.scv(scvf.outsideScvIdx()).center();
                const auto gOutside = problem.gravityAtPos(xOutside);

                // Sum rou*g*gradz
            Scalar Gravitywphasepotential= (xInside*gInside - xOutside*gOutside)*rho(wPhaseIdx);
            Scalar Gravitynphasepotential= (xInside*gInside - xOutside*gOutside)*rho(nPhaseIdx);
            Scalar Gravityphasepotential = Gravitywphasepotential+Gravitynphasepotential;
                // Average potential gradient at the interface

            Scalar AverageDp = (pnOutside - pnInside) + 0.5*Gravityphasepotential;
                //Calculate the coefficient Gammaij
            //1.Calculate the max|Gm,ij-0.5(g1,ij+g2,ij)|

            Scalar maxpart = std::max((Gravitywphasepotential-0.5*Gravityphasepotential), (Gravitynphasepotential-0.5*Gravityphasepotential));

            //in our case wenn S=0.5 we get e= 0.5, Sigma(suplambda)=2
            constexpr double pi = std::atan(1)*4 ;
            //avoid the issue maxpart very small
            Scalar Gammaij = std::min(10e6, pi*0.5/(2*maxpart));


                //Calculate the coefficient Beltaij
            Scalar Beltaij = 0.5 + 1/pi*std::atan(Gammaij*AverageDp);

            //Calculate the weighted mobility of each phase
            const auto mobW_WA = Beltaij*insideVolVars.mobility(wPhaseIdx) + (1-Beltaij)*outsideVolVars.mobility(wPhaseIdx);
            const auto mobN_WA = Beltaij*insideVolVars.mobility(nPhaseIdx) + (1-Beltaij)*outsideVolVars.mobility(nPhaseIdx);
            const auto mobT_WA = mobW_WA + mobN_WA;

            // the gravitational part of discretization of velocity
            Scalar WeightedGravityDp = mobW_WA*Gravitywphasepotential + mobN_WA*Gravitynphasepotential;

            // the capillary part of discretization of velocity
            // ! we first ignore the discontinuity of Pc and diffrent mobility relations
            Scalar WeighterCapillaryMobility = 0;
            Scalar WeightedCapillaryMobility = 2* (insideVolVars.mobility(nPhaseIdx)*outsideVolVars.mobility(nPhaseIdx))/(insideVolVars.mobility(nPhaseIdx) + outsideVolVars.mobility(nPhaseIdx) );
            // the capillary part of the discretization of velocity


            const auto pcInside = insideVolVars.capillaryPressure();
            const auto pcOutside = outsideVolVars.capillaryPressure();
            Scalar WeightedCapillaryDp = WeightedCapillaryMobility*(pcOutside - pcInside);

            static VelocityVector v_t = fluxVarsCache.advectionTij()*(mobT_WA*(pnInside - pnOutside)+ WeightedGravityDp + WeightedCapillaryDp);


            //////////////////////////////////////////////////////////////
            // The viscous flux before upwinding is the total velocity
            //////////////////////////////////////////////////////////////
            fluxes[viscousFluxIdx] = (v_t*scvf.unitOuterNormal())*scvf.area(); // TODO extrusion factor

            //////////////////////////////////////////////////////////////
            // Compute the capillary flux before upwinding
            //////////////////////////////////////////////////////////////
            static const bool useHybridUpwinding = getParamFromGroup<bool>(problem.paramGroup(), "Problem.UseHybridUpwinding");
            // for hybrid upwinding we want tij*delta_Sij
            if (useHybridUpwinding)
            {
                const auto sInside = insideVolVars.saturation(wPhaseIdx);
                const auto sOutside = outsideVolVars.saturation(wPhaseIdx);
                fluxes[capillaryFluxIdx] = fluxVarsCache.advectionTij()*(sInside-sOutside);
            }
            // for potetial phase upwinding we want tij*delta_Pcij
            else
            {
                const auto pcInside = insideVolVars.capillaryPressure();
                const auto pcOutside = outsideVolVars.capillaryPressure();
                fluxes[capillaryFluxIdx] = fluxVarsCache.advectionTij()*(pcOutside - pcInside);
            }

            //////////////////////////////////////////////////////////////
            // Compute the gravitational flux before upwinding
            //////////////////////////////////////////////////////////////

            static const bool enableGravity = getParamFromGroup<bool>(problem.paramGroup(), "Problem.EnableGravity");
            if (enableGravity)
            {
                // do averaging for the density over all neighboring elements
                const auto rho = [&](int phaseIdx)
                {
                    // boundaries
                    if (scvf.boundary())
                        return insideVolVars.density(phaseIdx);

                    // inner faces with two neighboring elements
                    else
                        return (insideVolVars.density(phaseIdx) + outsideVolVars.density(phaseIdx))*0.5;
                };

                // ask for the gravitational acceleration in the inside neighbor
                const auto xInside = insideScv.center();
                const auto gInside = problem.gravityAtPos(xInside);
                const auto xOutside = scvf.boundary() ? scvf.ipGlobal()
                                                      : fvGeometry.scv(scvf.outsideScvIdx()).center();
                const auto gOutside = problem.gravityAtPos(xOutside);

                fluxes[gravityFluxIdx] = fluxVarsCache.advectionTij()*(xInside*gInside - xOutside*gOutside)*(rho(nPhaseIdx) - rho(wPhaseIdx));
            }
            else // no gravity
            {
                fluxes[gravityFluxIdx] = 0.0;
            }
        }
        else
        {
            DUNE_THROW(Dune::InvalidStateException, "Unknown equation index!");
        }

        return fluxes;
    }
};

}; // end namespace Dumux

#endif
