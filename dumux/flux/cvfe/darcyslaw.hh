// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup CVFEFlux
 * \brief Specialization of Darcy's Law for control-volume finite element schemes
 *
 * This file contains the data which is required to calculate
 * volume and mass fluxes of fluid phases over a face of a finite volume by means
 * of the Darcy approximation.
 */
#ifndef DUMUX_DISCRETIZATION_CVFE_DARCYS_LAW_HH
#define DUMUX_DISCRETIZATION_CVFE_DARCYS_LAW_HH

#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/discretization/cvfe/quadraturerules.hh>
#include <dumux/flux/facetensoraverage.hh>

namespace Dumux {

/*!
 * \ingroup CVFEFlux
 * \brief Darcy's law for control-volume finite element schemes
 * \tparam Scalar the scalar type for scalar physical quantities
 * \tparam GridGeometry the grid geometry
 */
template<class Scalar, class GridGeometry>
class CVFEDarcysLaw
{
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using Extrusion = Extrusion_t<GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;

public:

    /*!
     * \brief Returns the advective flux of a fluid phase
     *        across the given sub-control volume face.
     * \note This assembles the term
     *       \f$-|\sigma| \mathbf{n}^T \mathbf{K} \left( \nabla p - \rho \mathbf{g} \right)\f$,
     *       where \f$|\sigma|\f$ is the area of the face and \f$\mathbf{n}\f$ is the outer
     *       normal vector. Thus, the flux is given in N*m, and can be converted
     *       into a volume flux (m^3/s) or mass flux (kg/s) by applying an upwind scheme
     *       for the mobility or the product of density and mobility, respectively.
     */
    template<class Problem, class ElementVolumeVariables, class ElementFluxVarsCache>
    static Scalar flux(const Problem& problem,
                       const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolumeFace& scvf,
                       const int phaseIdx,
                       const ElementFluxVarsCache& elemFluxVarCache)
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

        // an average over the two control volumes the face separates, and thus a property of
        // the face rather than of a point on it
        const auto K = faceTensorAverage(insideK, outsideK, scvf.unitOuterNormal());
        static const bool enableGravity = getParamFromGroup<bool>(problem.paramGroup(), "Problem.EnableGravity");

        // The quantity the loop below accumulates is a flux density at an interpolation
        // point. The form to arrive at is a flux law that offers exactly that, a
        // (context, ipData) -> density, and leaves the quadrature to the caller: the state
        // shared by all points, above all the interpolation over the local degrees of
        // freedom, is then computed once for all phases and components instead of once per
        // flux law call, as the free-flow momentum residual already does.
        Scalar flux = 0.0;
        for (const auto& quadPoint : CVFE::quadratureRule(fvGeometry, scvf))
        {
            const auto& faceIpData = quadPoint.ipData();
            const auto& fluxVarCache = cache(elemFluxVarCache, faceIpData);
            const auto& shapeValues = fluxVarCache.shapeValues();

            // evaluate gradP - rho*g at the integration point
            Dune::FieldVector<Scalar, dimWorld> gradP(0.0);
            Scalar rho(0.0);
            for (const auto& localDof : localDofs(fvGeometry))
            {
                const auto& volVars = elemVolVars[localDof];

                if (enableGravity)
                    rho += volVars.density(phaseIdx)*shapeValues[localDof.index()][0];

                // the global shape function gradient
                gradP.axpy(volVars.pressure(phaseIdx), fluxVarCache.gradN(localDof.index()));
            }

            if (enableGravity)
                gradP.axpy(-rho, problem.spatialParams().gravity(faceIpData.global()));

            // the quadrature weight carries the area of the face
            flux -= quadPoint.weight()*vtmv(faceIpData.unitOuterNormal(), K, gradP);
        }

        return flux;
    }

    // compute transmissibilities ti for analytical Jacobians
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

        const auto K = faceTensorAverage(insideK, outsideK, scvf.unitOuterNormal());

        std::vector<Scalar> ti(fvGeometry.numScv());
        for (const auto& scv : scvs(fvGeometry))
            ti[scv.indexInElement()] =
                -1.0*Extrusion::area(fvGeometry, scvf)*vtmv(scvf.unitOuterNormal(), K, fluxVarCache.gradN(scv.indexInElement()));

        return ti;
    }
};

// forward declaration
template<class TypeTag, class DiscretizationMethod>
class DarcysLawImplementation;

/*!
 * \ingroup CVFEFlux
 * \brief Specialization of Darcy's Law for control-volume finite element schemes
 */
template<class TypeTag, class DM>
class DarcysLawImplementation<TypeTag, DiscretizationMethods::CVFE<DM>>
: public CVFEDarcysLaw<GetPropType<TypeTag, Properties::Scalar>, GetPropType<TypeTag, Properties::GridGeometry>>
{};

} // end namespace Dumux

#endif
