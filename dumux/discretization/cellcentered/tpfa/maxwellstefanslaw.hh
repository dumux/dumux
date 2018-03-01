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
 *        diffusive mass fluxes due to molecular diffusion with Fick's law.
 */
#ifndef DUMUX_DISCRETIZATION_CC_TPFA_MAXWELL_STEFAN_LAW_HH
#define DUMUX_DISCRETIZATION_CC_TPFA_MAXWELL_STEFAN_LAW_HH

#include <dune/common/float_cmp.hh>

#include <dumux/common/math.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/discretization/methods.hh>
#include <dumux/discretization/fluxvariablescaching.hh>

namespace Dumux {

// forward declaration
template <class TypeTag, DiscretizationMethod discMethod>
class MaxwellStefansLawImplementation;

/*!
 * \ingroup CCTpfaMaxwellStefansLaw
 * \brief Specialization of Maxwell Stefan's Law for the CCTpfa method.
 */
template <class TypeTag>
class MaxwellStefansLawImplementation<TypeTag, DiscretizationMethod::cctpfa >
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using EffDiffModel = typename GET_PROP_TYPE(TypeTag, EffectiveDiffusivityModel);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using IndexType = typename GridView::IndexSet::IndexType;
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);
    using FluxVariablesCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);

    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;
    static const int numPhases = GET_PROP_VALUE(TypeTag, NumPhases);
    static const int numComponents = GET_PROP_VALUE(TypeTag,NumComponents);

    using DimWorldMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
    using ComponentFluxVector = Dune::FieldVector<Scalar, numComponents>;
    using ReducedComponentVector = Dune::FieldVector<Scalar, numComponents-1>;
    using ReducedComponentMatrix = Dune::FieldMatrix<Scalar, numComponents-1, numComponents-1>;

public:
    // state the discretization method this implementation belongs to
    static const DiscretizationMethod discMethod = DiscretizationMethod::cctpfa;

    //! state the type for the corresponding cache and its filler
    //! We don't cache anything for this law
    using Cache = FluxVariablesCaching::EmptyDiffusionCache;
    using CacheFiller = FluxVariablesCaching::EmptyCacheFiller;

    static ComponentFluxVector flux(const Problem& problem,
                                    const Element& element,
                                    const FVElementGeometry& fvGeometry,
                                    const ElementVolumeVariables& elemVolVars,
                                    const SubControlVolumeFace& scvf,
                                    const int phaseIdx,
                                    const ElementFluxVariablesCache& elemFluxVarsCache)
    {
        //this is to calculate the maxwellStefan diffusion in a multicomponent system.
        //see: Multicomponent Mass Transfer. R. Taylor u. R. Krishna. J. Wiley & Sons, New York 1993
        ComponentFluxVector componentFlux(0.0);
        ReducedComponentVector moleFracInside(0.0);
        ReducedComponentVector moleFracOutside(0.0);
        ReducedComponentVector reducedFlux(0.0);
        ReducedComponentMatrix reducedDiffusionMatrixInside(0.0);
        ReducedComponentMatrix reducedDiffusionMatrixOutside(0.0);

        // get inside/outside volume variables
        const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
        const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];
        const auto rhoInside = insideVolVars.molarDensity(phaseIdx);
        const auto rhoOutside = outsideVolVars.molarDensity(phaseIdx);
        //calculate the mole fraction vectors
        for (int compIdx = 0; compIdx < numComponents-1; compIdx++)
        {
            //calculate x_inside
            const auto xInside = insideVolVars.moleFraction(phaseIdx, compIdx);
            //calculate outside molefraction with the respective transmissibility
            const auto xOutside = outsideVolVars.moleFraction(phaseIdx, compIdx);

            moleFracInside[compIdx] = xInside;
            moleFracOutside[compIdx] = xOutside;
        }

        //we cannot solve that if the matrix is 0 everywhere
        if(!(insideVolVars.saturation(phaseIdx) == 0 || outsideVolVars.saturation(phaseIdx) == 0))
        {
            const auto insideScvIdx = scvf.insideScvIdx();
            const auto& insideScv = fvGeometry.scv(insideScvIdx);
            const Scalar omegai = calculateOmega_(scvf,
                                                insideScv,
                                                insideVolVars.extrusionFactor());

            //now we have to do the tpfa: J_i = J_j which leads to: tij(xi -xj) = -rho Bi^-1 omegai(x*-xi) with x* = (omegai Bi^-1 + omegaj Bj^-1)^-1 (xi omegai Bi^-1 + xj omegaj Bj^-1) with i inside and j outside
            reducedDiffusionMatrixInside = setupMSMatrix_(problem, element, fvGeometry, insideVolVars, insideScv, phaseIdx);

            //if on boundary
            if (scvf.boundary() || scvf.numOutsideScvs() > 1)
            {
                moleFracOutside -= moleFracInside;
                reducedDiffusionMatrixInside.solve(reducedFlux, moleFracOutside);
                reducedFlux *= omegai;


            }
            else //we need outside cells as well if we are not on the boundary
            {
                Scalar omegaj;
                const auto outsideScvIdx = scvf.outsideScvIdx();
                const auto& outsideScv = fvGeometry.scv(outsideScvIdx);

                reducedDiffusionMatrixOutside = setupMSMatrix_(problem, element, fvGeometry, outsideVolVars, outsideScv, phaseIdx);

                if (dim == dimWorld)
                    // assume the normal vector from outside is anti parallel so we save flipping a vector
                    omegaj = -1.0*calculateOmega_(scvf,
                                            outsideScv,
                                            outsideVolVars.extrusionFactor());
                else
                    omegaj = calculateOmega_(fvGeometry.flipScvf(scvf.index()),
                                            outsideScv,
                                            outsideVolVars.extrusionFactor());

                reducedDiffusionMatrixInside.invert();
                reducedDiffusionMatrixOutside.invert();
                reducedDiffusionMatrixInside *= omegai;
                reducedDiffusionMatrixOutside *= omegaj;

                //in the helpervector we store the values for x*
                ReducedComponentVector helperVector(0.0);
                ReducedComponentVector gradientVectori(0.0);
                ReducedComponentVector gradientVectorj(0.0);

                reducedDiffusionMatrixInside.mv(moleFracInside, gradientVectori);
                reducedDiffusionMatrixOutside.mv(moleFracOutside, gradientVectorj);

                auto gradientVectorij = (gradientVectori + gradientVectorj);

                //add the two matrixes to each other
                reducedDiffusionMatrixOutside += reducedDiffusionMatrixInside;

                reducedDiffusionMatrixOutside.solve(helperVector, gradientVectorij);

                //Bi^-1 omegai rho(x*-xi)
                helperVector -=moleFracInside;
                reducedDiffusionMatrixInside.mv(helperVector, reducedFlux);
            }

            reducedFlux *= -0.5*(rhoInside+rhoOutside)*scvf.area();
            for (int compIdx = 0; compIdx < numComponents-1; compIdx++)
            {
                componentFlux[compIdx] = reducedFlux[compIdx];
                componentFlux[numComponents-1] -=reducedFlux[compIdx];
            }
        }
        return componentFlux ;
    }

private:
   static Scalar calculateOmega_(const SubControlVolumeFace& scvf,
                                 const SubControlVolume &scv,
                                 const Scalar extrusionFactor)
    {
        auto distanceVector = scvf.ipGlobal();
        distanceVector -= scv.center();
        distanceVector /= distanceVector.two_norm2();

        Scalar omega = (distanceVector * scvf.unitOuterNormal());
        omega *= extrusionFactor;

        return omega;
    }

    static ReducedComponentMatrix setupMSMatrix_(const Problem& problem,
                                                const Element& element,
                                                const FVElementGeometry& fvGeometry,
                                                const VolumeVariables& volVars,
                                                const SubControlVolume& scv,
                                                const int phaseIdx)
    {
        ReducedComponentMatrix reducedDiffusionMatrix(0.0);

        //this is to not devide by 0 if the saturation in 0 and the effectiveDiffusivity becomes zero due to that
        if(volVars.saturation(phaseIdx) == 0)
            return reducedDiffusionMatrix;

        for (int compIIdx = 0; compIIdx < numComponents-1; compIIdx++)
        {
            const auto xi = volVars.moleFraction(phaseIdx, compIIdx);

            //calculate diffusivity for i,numComponents
            Scalar tin = getDiffusionCoefficient(phaseIdx, compIIdx, numComponents-1, problem, element, volVars, scv);
            tin = EffDiffModel::effectiveDiffusivity(volVars.porosity(), volVars.saturation(phaseIdx), tin);
            //set the entrys of the diffusion matrix of the diagonal
            reducedDiffusionMatrix[compIIdx][compIIdx] += xi/tin;
            for (int compkIdx = 0; compkIdx < numComponents; compkIdx++)
            {
                if (compkIdx == compIIdx)
                            continue;

                const auto xk = volVars.moleFraction(phaseIdx, compkIdx);
                Scalar tik = getDiffusionCoefficient(phaseIdx, compIIdx, compkIdx, problem, element, volVars, scv);
                tik = EffDiffModel::effectiveDiffusivity(volVars.porosity(), volVars.saturation(phaseIdx), tik);
                reducedDiffusionMatrix[compIIdx][compIIdx] += xk/tik;
            }

            // now set the rest of the entries (off-diagonal)
            for (int compJIdx = 0; compJIdx < numComponents-1; compJIdx++)
            {
                //we don't want to calculate e.g. water in water diffusion
                if (compIIdx == compJIdx)
                    continue;
                //calculate diffusivity for compIIdx, compJIdx
                Scalar tij = getDiffusionCoefficient(phaseIdx, compIIdx, compJIdx, problem, element, volVars, scv);
                tij = EffDiffModel::effectiveDiffusivity(volVars.porosity(), volVars.saturation(phaseIdx), tij);
                reducedDiffusionMatrix[compIIdx][compJIdx] +=xi*(1/tin - 1/tij);
            }
        }
        return reducedDiffusionMatrix;
    }

    template <class T = TypeTag, typename std::enable_if_t<GET_PROP_TYPE(T, FluidSystem)::isTracerFluidSystem(), int> =0 >
    static Scalar getDiffusionCoefficient(const int phaseIdx,
                            const int compIIdx,
                            const int compJIdx,
                            const Problem& problem,
                            const Element& element,
                            const VolumeVariables& volVars,
                            const SubControlVolume& scv)
    {
        return FluidSystem::binaryDiffusionCoefficient(compIIdx,
                                                       compJIdx,
                                                       problem,
                                                       element,
                                                       scv);
    }

    template <class T = TypeTag, typename std::enable_if_t<!GET_PROP_TYPE(T, FluidSystem)::isTracerFluidSystem(), int> =0 >
    static Scalar getDiffusionCoefficient(const int phaseIdx,
                            const int compIIdx,
                            const int compJIdx,
                            const Problem& problem,
                            const Element& element,
                            const VolumeVariables& volVars,
                            const SubControlVolume& scv)
    {
        auto fluidState = volVars.fluidState();
        typename FluidSystem::ParameterCache paramCache;
        paramCache.updateAll(fluidState);
        return FluidSystem::binaryDiffusionCoefficient(fluidState,
                                                       paramCache,
                                                       phaseIdx,
                                                       compIIdx,
                                                       compJIdx);
    }



};
} // end namespace

#endif
