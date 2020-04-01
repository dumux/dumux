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
 * \ingroup OnePModel
 * \brief Element-wise calculation of the residual and its derivatives
 *        for a single-phase, incompressible, test problem.
 */

#ifndef DUMUX_1P_INCOMPRESSIBLE_LOCAL_RESIDUAL_HH
#define DUMUX_1P_INCOMPRESSIBLE_LOCAL_RESIDUAL_HH

#include <dumux/discretization/method.hh>
#include <dumux/porousmediumflow/immiscible/localresidual.hh>

namespace Dumux {

/*!
 * \ingroup OnePModel
 * \brief Element-wise calculation of the residual and its derivatives
 *        for a single-phase, incompressible, test problem.
 */
template<class TypeTag>
class OnePIncompressibleLocalResidual : public ImmiscibleLocalResidual<TypeTag>
{
    using ParentType = ImmiscibleLocalResidual<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using FluxVariables = GetPropType<TypeTag, Properties::FluxVariables>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using ElementFluxVariablesCache = typename GetPropType<TypeTag, Properties::GridFluxVariablesCache>::LocalView;
    using FVElementGeometry = typename GetPropType<TypeTag, Properties::GridGeometry>::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    // first index for the mass balance
    enum { conti0EqIdx = Indices::conti0EqIdx };
    enum { pressureIdx = Indices::pressureIdx };

    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;

public:
    using ParentType::ParentType;

    template<class PartialDerivativeMatrix>
    void addStorageDerivatives(PartialDerivativeMatrix& partialDerivatives,
                               const Problem& problem,
                               const Element& element,
                               const FVElementGeometry& fvGeometry,
                               const VolumeVariables& curVolVars,
                               const SubControlVolume& scv) const {}

    template<class PartialDerivativeMatrix>
    void addSourceDerivatives(PartialDerivativeMatrix& partialDerivatives,
                              const Problem& problem,
                              const Element& element,
                              const FVElementGeometry& fvGeometry,
                              const VolumeVariables& curVolVars,
                              const SubControlVolume& scv) const
    {
        problem.addSourceDerivatives(partialDerivatives, element, fvGeometry, curVolVars, scv);
    }

    //! Flux derivatives for the cell-centered tpfa scheme
    template<class PartialDerivativeMatrices, class T = TypeTag>
    std::enable_if_t<GetPropType<T, Properties::GridGeometry>::discMethod == DiscretizationMethod::cctpfa, void>
    addFluxDerivatives(PartialDerivativeMatrices& derivativeMatrices,
                       const Problem& problem,
                       const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& curElemVolVars,
                       const ElementFluxVariablesCache& elemFluxVarsCache,
                       const SubControlVolumeFace& scvf) const
    {
        static_assert(!FluidSystem::isCompressible(0),
                      "1p/incompressiblelocalresidual.hh: Only incompressible fluids are allowed!");
        static_assert(FluidSystem::viscosityIsConstant(0),
                      "1p/incompressiblelocalresidual.hh: Only fluids with constant viscosities are allowed!");

        // we know the "upwind factor" is constant, get inner one here and compute derivatives
        static const Scalar up = curElemVolVars[scvf.insideScvIdx()].density()
                                 / curElemVolVars[scvf.insideScvIdx()].viscosity();

        // for network grids we have to do something else
        if (dim < dimWorld && scvf.numOutsideScvs() > 1)
        {
            const auto insideTi = elemFluxVarsCache[scvf].advectionTij();
            const auto sumTi = [&]
            {
                Scalar sumTi(insideTi);
                for (unsigned int i = 0; i < scvf.numOutsideScvs(); ++i)
                {
                    const auto& flippedScvf = fvGeometry.flipScvf(scvf.index(), i);
                    const auto& outsideFluxVarsCache = elemFluxVarsCache[flippedScvf];
                    sumTi += outsideFluxVarsCache.advectionTij();
                }
                return sumTi;
            }();

            // add partial derivatives to the respective given matrices
            derivativeMatrices[scvf.insideScvIdx()][conti0EqIdx][pressureIdx]
                += up*insideTi*(1.0 - insideTi/sumTi);
            for (unsigned int i = 0; i < scvf.numOutsideScvs(); ++i)
                derivativeMatrices[scvf.outsideScvIdx(i)][conti0EqIdx][pressureIdx]
                    += -up*insideTi*elemFluxVarsCache[fvGeometry.flipScvf(scvf.index(), i)].advectionTij()/sumTi;

        }
        else
        {
            const auto deriv = elemFluxVarsCache[scvf].advectionTij()*up;
            // add partial derivatives to the respective given matrices
            derivativeMatrices[scvf.insideScvIdx()][conti0EqIdx][pressureIdx] += deriv;
            derivativeMatrices[scvf.outsideScvIdx()][conti0EqIdx][pressureIdx] -= deriv;
        }
    }

    //! Flux derivatives for the cell-centered mpfa scheme
    template<class PartialDerivativeMatrices, class T = TypeTag>
    std::enable_if_t<GetPropType<T, Properties::GridGeometry>::discMethod == DiscretizationMethod::ccmpfa, void>
    addFluxDerivatives(PartialDerivativeMatrices& derivativeMatrices,
                       const Problem& problem,
                       const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& curElemVolVars,
                       const ElementFluxVariablesCache& elemFluxVarsCache,
                       const SubControlVolumeFace& scvf) const
    {
        static_assert(!FluidSystem::isCompressible(0),
                      "1p/incompressiblelocalresidual.hh: Only incompressible fluids are allowed!");
        static_assert(FluidSystem::viscosityIsConstant(0),
                      "1p/incompressiblelocalresidual.hh: Only fluids with constant viscosities are allowed!");

        // we know the "upwind factor" is constant, get inner one here and compute derivatives
        static const Scalar up = curElemVolVars[scvf.insideScvIdx()].density()
                                 / curElemVolVars[scvf.insideScvIdx()].viscosity();

        const auto& fluxVarsCache = elemFluxVarsCache[scvf];
        const auto localFaceIdx = fluxVarsCache.ivLocalFaceIndex();
        const auto usesSecondary = fluxVarsCache.usesSecondaryIv();
        const auto switchSign = fluxVarsCache.advectionSwitchFluxSign();

        const auto& stencil = fluxVarsCache.advectionStencil();
        const auto& tij = usesSecondary ? fluxVarsCache.advectionSecondaryDataHandle().T()[localFaceIdx]
                                        : fluxVarsCache.advectionPrimaryDataHandle().T()[localFaceIdx];

        // We assume same the tij are order as the stencil up to stencil.size()
        // any contribution of Dirichlet BCs is assumed to be placed afterwards
        assert(stencil.size() <= tij.size());
        for (unsigned int i = 0; i < stencil.size();++i)
            derivativeMatrices[stencil[i]][conti0EqIdx][pressureIdx] += switchSign ? -tij[i]*up
                                                                                   :  tij[i]*up;
    }

    //! Flux derivatives for the box scheme
    template<class JacobianMatrix, class T = TypeTag>
    std::enable_if_t<GetPropType<T, Properties::GridGeometry>::discMethod == DiscretizationMethod::box, void>
    addFluxDerivatives(JacobianMatrix& A,
                       const Problem& problem,
                       const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& curElemVolVars,
                       const ElementFluxVariablesCache& elemFluxVarsCache,
                       const SubControlVolumeFace& scvf) const
    {
        static_assert(!FluidSystem::isCompressible(0),
                      "1p/incompressiblelocalresidual.hh: Only incompressible fluids are allowed!");
        static_assert(FluidSystem::viscosityIsConstant(0),
                      "1p/incompressiblelocalresidual.hh: Only fluids with constant viscosities are allowed!");

        using AdvectionType = GetPropType<T, Properties::AdvectionType>;
        const auto ti = AdvectionType::calculateTransmissibilities(problem,
                                                                   element,
                                                                   fvGeometry,
                                                                   curElemVolVars,
                                                                   scvf,
                                                                   elemFluxVarsCache[scvf]);

        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& outsideScv = fvGeometry.scv(scvf.outsideScvIdx());

        // we know the "upwind factor" is constant, get inner one here and compute derivatives
        static const Scalar up = curElemVolVars[scvf.insideScvIdx()].density()
                                 / curElemVolVars[scvf.insideScvIdx()].viscosity();
        for (const auto& scv : scvs(fvGeometry))
        {
            auto d = up*ti[scv.indexInElement()];
            A[insideScv.dofIndex()][scv.dofIndex()][conti0EqIdx][pressureIdx] += d;
            A[outsideScv.dofIndex()][scv.dofIndex()][conti0EqIdx][pressureIdx] -= d;
        }
    }

    //! Dirichlet flux derivatives for the cell-centered tpfa scheme
    template<class PartialDerivativeMatrices, class T = TypeTag>
    std::enable_if_t<GetPropType<T, Properties::GridGeometry>::discMethod == DiscretizationMethod::cctpfa, void>
    addCCDirichletFluxDerivatives(PartialDerivativeMatrices& derivativeMatrices,
                                  const Problem& problem,
                                  const Element& element,
                                  const FVElementGeometry& fvGeometry,
                                  const ElementVolumeVariables& curElemVolVars,
                                  const ElementFluxVariablesCache& elemFluxVarsCache,
                                  const SubControlVolumeFace& scvf) const
    {
        // we know the "upwind factor" is constant, get inner one here
        static const Scalar up = curElemVolVars[scvf.insideScvIdx()].density()
                                 / curElemVolVars[scvf.insideScvIdx()].viscosity();
        const auto deriv = elemFluxVarsCache[scvf].advectionTij()*up;

        // compute and add partial derivative to the respective given matrices
        derivativeMatrices[scvf.insideScvIdx()][conti0EqIdx][pressureIdx] += deriv;
    }

    //! Dirichlet flux derivatives for the cell-centered mpfa scheme
    template<class PartialDerivativeMatrices, class T = TypeTag>
    std::enable_if_t<GetPropType<T, Properties::GridGeometry>::discMethod == DiscretizationMethod::ccmpfa, void>
    addCCDirichletFluxDerivatives(PartialDerivativeMatrices& derivativeMatrices,
                                  const Problem& problem,
                                  const Element& element,
                                  const FVElementGeometry& fvGeometry,
                                  const ElementVolumeVariables& curElemVolVars,
                                  const ElementFluxVariablesCache& elemFluxVarsCache,
                                  const SubControlVolumeFace& scvf) const
    {
        addFluxDerivatives(derivativeMatrices, problem, element, fvGeometry, curElemVolVars, elemFluxVarsCache, scvf);
    }

    //! Robin-type flux derivatives
    template<class PartialDerivativeMatrices>
    void addRobinFluxDerivatives(PartialDerivativeMatrices& derivativeMatrices,
                                 const Problem& problem,
                                 const Element& element,
                                 const FVElementGeometry& fvGeometry,
                                 const ElementVolumeVariables& curElemVolVars,
                                 const ElementFluxVariablesCache& elemFluxVarsCache,
                                 const SubControlVolumeFace& scvf) const
    {
        //! Robin-type boundary conditions are problem-specific.
        //! We can't put a general implementation here - users defining Robin-type BCs
        //! while using analytical Jacobian assembly must overload this function!
    }
};

} // end namespace Dumux

#endif
