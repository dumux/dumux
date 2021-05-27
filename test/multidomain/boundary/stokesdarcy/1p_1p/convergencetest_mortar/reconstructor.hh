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
 * \ingroup TODO doc me.
 * \brief TODO doc me.
 */
#ifndef DUMUX_MORTAR_RECONSTRUCTOR_HH
#define DUMUX_MORTAR_RECONSTRUCTOR_HH

#include <algorithm>
#include <vector>
#include <utility>
#include <unordered_map>

#include <dune/common/exceptions.hh>

#include <dumux/common/parameters.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/cellcentered/tpfa/computetransmissibility.hh>

#include "reconstructionhelper.hh"

namespace Dumux {

/*!
 * \ingroup TODO doc me.
 * \brief TODO doc me.
 */
template< class SubDomainTraits,
          DiscretizationMethod dm = SubDomainTraits::GridGeometry::discMethod >
class MortarReconstructor;

/*!
 * \ingroup TODO doc me.
 * \brief TODO doc me.
 */
template<class SubDomainTraits>
class MortarReconstructor< SubDomainTraits, DiscretizationMethod::staggered>
{
    using GridGeometry = typename SubDomainTraits::GridGeometry;
    using GridVariables = typename SubDomainTraits::GridVariables;
    using SolutionVector = typename SubDomainTraits::SolutionVector;

    using ElementScvfIndexMap = typename MortarReconstructionHelper::ElementScvfIndexMap<GridGeometry>;

public:
    /*!
     * \brief TODO doc me.
     */
    template<class MortarSolutionVector>
    static MortarSolutionVector recoverNormalFlux(const GridGeometry& gridGeometry,
                                                  const GridVariables& gridVariables,
                                                  const SolutionVector& sol,
                                                  const ElementScvfIndexMap& map)
    {
        // set up vector with interface flux per cell
        MortarSolutionVector cellFluxes;
        cellFluxes.resize(gridGeometry.gridView().size(0));
        cellFluxes = 0.0;

        for (const auto& element : elements(gridGeometry.gridView()))
        {
            auto fvGeometry = localView(gridGeometry);
            fvGeometry.bindElement(element);

            const auto eIdx = gridGeometry.elementMapper().index(element);
            if (map.find(eIdx) != map.end())
                for (auto scvfIdx : map.at(eIdx))
                {
                    const auto scvf = fvGeometry.scvf(scvfIdx);
                    static constexpr auto faceIdx = GridGeometry::faceIdx();
                    cellFluxes[eIdx] = sol[faceIdx][scvf.dofIndex()]*scvf.directionSign();
                }
        }

        return cellFluxes;
    }

    /*!
     * \brief TODO doc me.
     */
    template<class MortarSolutionVector>
    static MortarSolutionVector recoverPressure(const GridGeometry& gridGeometry,
                                                const GridVariables& gridVariables,
                                                const SolutionVector& sol,
                                                const ElementScvfIndexMap& map)
    {
        // set up vector with interface pressure per cell
        MortarSolutionVector cellPressure;
        cellPressure.resize(gridGeometry.gridView().size(0));
        cellPressure = 0.0;

        for (const auto& element : elements(gridGeometry.gridView()))
        {
            auto fvGeometry = localView(gridGeometry);
            auto elemVolVars = localView(gridVariables.curGridVolVars());
            auto elemFaceVars = localView(gridVariables.curGridFaceVars());

            fvGeometry.bind(element);
            elemVolVars.bind(element, fvGeometry, sol);
            elemFaceVars.bind(element, fvGeometry, sol);

            const auto eIdx = gridGeometry.elementMapper().index(element);
            if (map.find(eIdx) != map.end())
                for (auto scvfIdx : map.at(eIdx))
                {
                    const auto scvf = fvGeometry.scvf(scvfIdx);
                    typename SubDomainTraits::FluxVariables fluxVars;
                    auto stress = fluxVars.computeMomentumFlux(gridVariables.curGridVolVars().problem(),
                                                               element,
                                                               scvf,
                                                               fvGeometry,
                                                               elemVolVars,
                                                               elemFaceVars,
                                                               gridVariables.gridFluxVarsCache());
                    stress /= (scvf.area()*elemVolVars[scvf.insideScvIdx()].extrusionFactor());
                    cellPressure[eIdx] = stress;
                }
        }

        return cellPressure;
    }
};

/*!
 * \ingroup TODO doc me.
 * \brief TODO doc me.
 */
template<class SubDomainTraits>
class MortarReconstructor< SubDomainTraits, DiscretizationMethod::box >
{
    using GridGeometry = typename SubDomainTraits::GridGeometry;
    using GridVariables = typename SubDomainTraits::GridVariables;
    using SolutionVector = typename SubDomainTraits::SolutionVector;

    using ElementScvfIndexMap = typename MortarReconstructionHelper::ElementScvfIndexMap<GridGeometry>;

public:
    /*!
     * \brief TODO doc me.
     */
    template<class MortarSolutionVector>
    static MortarSolutionVector recoverNormalFlux(const GridGeometry& gridGeometry,
                                                  const GridVariables& gridVariables,
                                                  const SolutionVector& sol,
                                                  const ElementScvfIndexMap& map)
    {
        // set up vector with interface flux per cell
        MortarSolutionVector cellFluxes;
        cellFluxes.resize(gridGeometry.gridView().size(0));
        cellFluxes = 0.0;

        for (const auto& entry : map)
        {
            const auto sdElemIdx = entry.first;
            const auto sdScvfIndices = entry.second;
            const auto sdElement = gridGeometry.element(sdElemIdx);

            auto fvGeometry = localView(gridGeometry);
            auto elemVolVars = localView(gridVariables.curGridVolVars());
            auto elemFluxVarsCache = localView(gridVariables.gridFluxVarsCache());

            fvGeometry.bind(sdElement);
            elemVolVars.bind(sdElement, fvGeometry, sol);
            elemFluxVarsCache.bind(sdElement, fvGeometry, elemVolVars);

            typename GridGeometry::GridView::ctype elemFluxArea = 0.0;
            for (const auto scvfIdx : sdScvfIndices)
            {
                const auto& scvf = fvGeometry.scvf(scvfIdx);
                const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
                const auto& insideVolVars = elemVolVars[insideScv];
                const auto& problem = elemVolVars.gridVolVars().problem();

                // flux is f = -rho/mu K * (gradp - rhog)
                static const bool enableGravity = getParamFromGroup<bool>(problem.paramGroup(), "Problem.EnableGravity");

                const auto& fluxVarCache = elemFluxVarsCache[scvf];
                const auto& shapeValues = fluxVarCache.shapeValues();

                // evaluate gradP - rho*g at integration point
                using Scalar = typename MortarSolutionVector::field_type;
                using GridView = typename GridGeometry::GridView;
                Dune::FieldVector<Scalar, GridView::dimensionworld> gradP(0.0);
                Scalar rho(0.0);
                for (auto&& scv : scvs(fvGeometry))
                {
                    const auto& volVars = elemVolVars[scv];
                    if (enableGravity)
                        rho += volVars.density()*shapeValues[scv.indexInElement()][0];

                    // the global shape function gradient
                    gradP.axpy(volVars.pressure(), fluxVarCache.gradN(scv.indexInElement()));
                }

                if (enableGravity)
                    gradP.axpy(-rho, problem.spatialParams().gravity(scvf.ipGlobal()));

                // apply the permeability
                Scalar flux = -1.0*vtmv(scvf.unitOuterNormal(), insideVolVars.permeability(), gradP)
                                  *scvf.area() *insideVolVars.mobility();

                cellFluxes[sdElemIdx] += flux;
                elemFluxArea += scvf.area();
            }

            // turn into velocities
            cellFluxes[sdElemIdx] /= elemFluxArea;
        }

        return cellFluxes;
    }


    /*!
     * \brief TODO doc me.
     */
    template<class MortarSolutionVector>
    static MortarSolutionVector
    recoverPressure(const GridGeometry& gridGeometry,
                    const GridVariables& gridVariables,
                    const SolutionVector& sol,
                    const ElementScvfIndexMap& map)
    {
        MortarSolutionVector result = sol;
        return result;
    }
};

/*!
 * \ingroup TODO doc me.
 * \brief TODO doc me.
 */
template<class SubDomainTraits>
class MortarReconstructor< SubDomainTraits, DiscretizationMethod::cctpfa >
{
    using GridGeometry = typename SubDomainTraits::GridGeometry;
    using GridVariables = typename SubDomainTraits::GridVariables;
    using SolutionVector = typename SubDomainTraits::SolutionVector;

    using ElementScvfIndexMap = typename MortarReconstructionHelper::ElementScvfIndexMap<GridGeometry>;

public:

    /*!
     * \brief TODO doc me.
     */
    template<class MortarSolutionVector>
    static MortarSolutionVector recoverPressure(const GridGeometry& gridGeometry,
                                                const GridVariables& gridVariables,
                                                const SolutionVector& sol,
                                                const ElementScvfIndexMap& map)
    {
        // set up pressure vector with reconstructed face pressures
        MortarSolutionVector sdPressures;
        sdPressures.resize(gridGeometry.gridView().size(0));
        sdPressures = 0.0;

        for (const auto& entry : map)
        {
            const auto sdElemIdx = entry.first;
            const auto scvfIdx = entry.second[0];
            const auto sdElement = gridGeometry.element(sdElemIdx);

            auto fvGeometry = localView(gridGeometry);
            auto elemVolVars = localView(gridVariables.curGridVolVars());
            auto elemFluxVarsCache = localView(gridVariables.gridFluxVarsCache());

            fvGeometry.bind(sdElement);
            elemVolVars.bind(sdElement, fvGeometry, sol);
            elemFluxVarsCache.bind(sdElement, fvGeometry, elemVolVars);

            const auto& scvf = fvGeometry.scvf(scvfIdx);
            const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
            const auto& insideVolVars = elemVolVars[insideScv];
            const auto& problem = elemVolVars.gridVolVars().problem();

            // flux = -ti*(pBoundary - pCell)
            auto ti = computeTpfaTransmissibility(scvf,
                                                  insideScv,
                                                  insideVolVars.permeability(),
                                                  insideVolVars.extrusionFactor());
            ti *= insideVolVars.density()*insideVolVars.mobility();

            auto flux = problem.neumann(sdElement, fvGeometry, elemVolVars, elemFluxVarsCache, scvf)[0];
            flux -= ti*insideVolVars.pressure();
            sdPressures[sdElemIdx] = -1.0*flux/ti;
        }

        return sdPressures;
    }

    /*!
     * \brief TODO doc me.
     */
    template<class MortarSolutionVector>
    static MortarSolutionVector recoverNormalFlux(const GridGeometry& gridGeometry,
                                                  const GridVariables& gridVariables,
                                                  const SolutionVector& sol,
                                                  const ElementScvfIndexMap& map)
    {
        // set up vector with interface flux per cell
        MortarSolutionVector cellFluxes;
        cellFluxes.resize(gridGeometry.gridView().size(0));
        cellFluxes = 0.0;

        for (const auto& entry : map)
        {
            const auto sdElemIdx = entry.first;
            const auto scvfIdx = entry.second[0];
            const auto sdElement = gridGeometry.element(sdElemIdx);

            auto fvGeometry = localView(gridGeometry);
            auto elemVolVars = localView(gridVariables.curGridVolVars());
            auto elemFluxVarsCache = localView(gridVariables.gridFluxVarsCache());

            fvGeometry.bind(sdElement);
            elemVolVars.bind(sdElement, fvGeometry, sol);
            elemFluxVarsCache.bind(sdElement, fvGeometry, elemVolVars);

            const auto& scvf = fvGeometry.scvf(scvfIdx);
            const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
            const auto& insideVolVars = elemVolVars[insideScv];
            const auto& problem = elemVolVars.gridVolVars().problem();

            typename SubDomainTraits::FluxVariables fluxVars;
            fluxVars.init(problem, sdElement, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);

            auto up = [] (const auto& vv) { return vv.mobility(); };
            cellFluxes[sdElemIdx] = fluxVars.advectiveFlux(0, up)
                                    /insideVolVars.extrusionFactor()
                                    /scvf.area();
        }

        return cellFluxes;
    }
};

} // end namespace Dumux

#endif
