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
 * \ingroup BoxFlux
 * \brief Specialization of Darcy's Law for the box method.
 *
 * This file contains the data which is required to calculate
 * volume and mass fluxes of fluid phases over a face of a finite volume by means
 * of the Darcy approximation.
 */
#ifndef DUMUX_DISCRETIZATION_BOX_DARCYS_LAW_HH
#define DUMUX_DISCRETIZATION_BOX_DARCYS_LAW_HH

#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/method.hh>

namespace Dumux {

// forward declaration
template<class TypeTag, DiscretizationMethod discMethod>
class DarcysLawImplementation;

// forward declaration
template<class Scalar, class FVGridGeometry>
class BoxDarcysLaw;

/*!
 * \ingroup BoxFlux
 * \brief Specialization of Darcy's Law for the box method.
 */
template<class TypeTag>
class DarcysLawImplementation<TypeTag, DiscretizationMethod::box>
: public BoxDarcysLaw<GetPropType<TypeTag, Properties::Scalar>, GetPropType<TypeTag, Properties::FVGridGeometry>>
{ };

/*!
 * \ingroup BoxFlux
 * \brief Darcy's law for box schemes
 * \tparam Scalar the scalar type for scalar physical quantities
 * \tparam FVGridGeometry the grid geometry
 */
template<class Scalar, class FVGridGeometry>
class BoxDarcysLaw
{
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GridView = typename FVGridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using CoordScalar = typename GridView::ctype;

    enum { dim = GridView::dimension};
    enum { dimWorld = GridView::dimensionworld};

    using DimWorldMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;

public:

    template<class Problem, class ElementVolumeVariables, class ElementFluxVarsCache>
    static Scalar flux(const Problem& problem,
                       const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolumeFace& scvf,
                       const int phaseIdx,
                       const ElementFluxVarsCache& elemFluxVarCache)
    {
        const auto& fluxVarCache = elemFluxVarCache[scvf];
        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& outsideScv = fvGeometry.scv(scvf.outsideScvIdx());
        const auto& insideVolVars = elemVolVars[insideScv];
        const auto& outsideVolVars = elemVolVars[outsideScv];

        auto insideK = insideVolVars.permeability();
        auto outsideK = outsideVolVars.permeability();

        // scale with correct extrusion factor
        insideK *= insideVolVars.extrusionFactor();
        outsideK *= outsideVolVars.extrusionFactor();

        const auto K = problem.spatialParams().harmonicMean(insideK, outsideK, scvf.unitOuterNormal());
        static const bool enableGravity = getParamFromGroup<bool>(problem.paramGroup(), "Problem.EnableGravity");

        const auto& shapeValues = fluxVarCache.shapeValues();

        // evaluate gradP - rho*g at integration point
        Dune::FieldVector<Scalar, dimWorld> gradP(0.0);
        Scalar rho(0.0);
        for (auto&& scv : scvs(fvGeometry))
        {
            const auto& volVars = elemVolVars[scv];

            if (enableGravity)
                rho += volVars.density(phaseIdx)*shapeValues[scv.indexInElement()][0];

            // the global shape function gradient
            gradP.axpy(volVars.pressure(phaseIdx), fluxVarCache.gradN(scv.indexInElement()));
        }

        if (enableGravity)
            gradP.axpy(-rho, problem.gravityAtPos(scvf.center()));

        // apply the permeability and return the flux
        return -1.0*vtmv(scvf.unitOuterNormal(), K, gradP)*scvf.area();
    }

    // compute transmissibilities ti for analytical jacobians
    template<class Problem, class ElementVolumeVariables, class FluxVarCache>
    static std::vector<Scalar> calculateTransmissibilities(const Problem& problem,
                                                           const Element& element,
                                                           const FVElementGeometry& fvGeometry,
                                                           const ElementVolumeVariables& elemVolVars,
                                                           const SubControlVolumeFace& scvf,
                                                           const FluxVarCache& fluxVarCache)
    {
        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& outsideScv = fvGeometry.scv(scvf.outsideScvIdx());
        const auto& insideVolVars = elemVolVars[insideScv];
        const auto& outsideVolVars = elemVolVars[outsideScv];

        auto insideK = insideVolVars.permeability();
        auto outsideK = outsideVolVars.permeability();

        // scale with correct extrusion factor
        insideK *= insideVolVars.extrusionFactor();
        outsideK *= outsideVolVars.extrusionFactor();

        const auto K = problem.spatialParams().harmonicMean(insideK, outsideK, scvf.unitOuterNormal());

        std::vector<Scalar> ti(fvGeometry.numScv());
        for (const auto& scv : scvs(fvGeometry))
            ti[scv.indexInElement()] =
                -1.0*scvf.area()*vtmv(scvf.unitOuterNormal(), K, fluxVarCache.gradN(scv.indexInElement()));

        return ti;
    }
};

} // end namespace Dumux

#endif
