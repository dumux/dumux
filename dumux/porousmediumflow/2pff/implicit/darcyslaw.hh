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


            Scalar WeightedDp, WeightedViscousDp, WeightedGravityDp, WeightedCapillaryDp, AverageDp = 0;
            Scalar Beltaij = 0;



             const auto pInside = insideVolVars.pressure(phaseIdx);
             const auto pOutside = outsideVolVars.pressure(phaseIdx);

             const auto xInside = insideScv.center();
             const auto xOutside = scvf.boundary() ? scvf.ipGlobal()
                                                      : fvGeometry.scv(scvf.outsideScvIdx()).center();

             const auto g = problem.gravityAtPos(xInside);

            Scalar AverageDp = (pOutside - pInside) + 0.5*g*(insideVolVars.density(0)*(xInside - xOutside) + insideVolVars.density(1)*(xInside - xOutside));


            constexpr double pi() { return std::atan(1)*4; };
            Scalar Beltaij = 0.5 + 1/pi*std::atan(Gammaij*AverageDp);
            const auto mobW_WA = Beltaij*insideVolVars.mobility(wPhaseIdx) + (1-Beltaij)*outsideVolVars.mobility(wPhaseIdx);
            const auto mobN_WA = Beltaij*insideVolVars.mobility(nPhaseIdx) + (1-Beltaij)*outsideVolVars.mobility(nPhaseIdx);
            //TODO:calculate Gammaij





            Scalar VelocityVektor v_t = fluxVarsCache.advectionTij()*((mobW_WA+mobN_WA)*(pInside - pOutside)+(mobW_WA*insideVolVars.density(0)*(xInside - xOutside)+mobN_WA*insideVolVars.density(1)*(xInside - xOutside)));

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

} // end namespace Dumux

#endif
