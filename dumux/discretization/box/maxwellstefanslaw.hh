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
#ifndef DUMUX_DISCRETIZATION_BOX_MAXWELL_STEFAN_LAW_HH
#define DUMUX_DISCRETIZATION_BOX_MAXWELL_STEFAN_LAW_HH

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
 * \brief Specialization of Maxwell Stefan's Law for the Box method.
 */
template <class TypeTag>
class MaxwellStefansLawImplementation<TypeTag, DiscretizationMethod::box >
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, GridVolumeVariables)::LocalView;
    using ElementFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, GridFluxVariablesCache)::LocalView;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using IndexType = typename GridView::IndexSet::IndexType;
    using Indices = typename GET_PROP_TYPE(TypeTag, ModelTraits)::Indices;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    enum { dim = GridView::dimension} ;
    enum { dimWorld = GridView::dimensionworld} ;
    enum
    {
        numPhases = GET_PROP_TYPE(TypeTag, ModelTraits)::numPhases(),
        numComponents = GET_PROP_TYPE(TypeTag, ModelTraits)::numComponents()
    };
    using DimWorldMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;
    using ComponentFluxVector = Dune::FieldVector<Scalar, numComponents>;
    using ReducedComponentVector = Dune::FieldVector<Scalar, numComponents-1>;
    using ReducedComponentMatrix = Dune::FieldMatrix<Scalar, numComponents-1, numComponents-1>;

public:
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
        ReducedComponentMatrix reducedDiffusionMatrix(0.0);
        ReducedComponentVector reducedFlux(0.0);
        ComponentFluxVector moleFrac(0.0);
        ReducedComponentVector normalX(0.0);

        // evaluate gradX at integration point and interpolate density
        const auto& fluxVarsCache = elemFluxVarsCache[scvf];
        const auto& shapeValues = fluxVarsCache.shapeValues();

        Scalar rho(0.0);
        for (auto&& scv : scvs(fvGeometry))
        {
            const auto& volVars = elemVolVars[scv];

            // density interpolation
            rho +=  volVars.molarDensity(phaseIdx)*shapeValues[scv.indexInElement()][0];

            //interpolate the mole fraction for the diffusion matrix
            for (int compIdx = 0; compIdx < numComponents; compIdx++)
            {
              moleFrac[compIdx] += volVars.moleFraction(phaseIdx, compIdx)*shapeValues[scv.indexInElement()][0];
            }
        }

        reducedDiffusionMatrix = setupMSMatrix_(problem, element, fvGeometry, elemVolVars, scvf, phaseIdx, moleFrac);

        for (int compIdx = 0; compIdx < numComponents-1; compIdx++)
        {
            Dune::FieldVector<Scalar, dimWorld> gradX(0.0);
            for (auto&& scv : scvs(fvGeometry))
            {
                const auto& volVars = elemVolVars[scv];

                // the mole/mass fraction gradient
                gradX.axpy(volVars.moleFraction(phaseIdx, compIdx), fluxVarsCache.gradN(scv.indexInElement()));
            }

           normalX[compIdx] = gradX *scvf.unitOuterNormal();
        }
         reducedDiffusionMatrix.solve(reducedFlux,normalX);
         reducedFlux *= -1.0*rho*scvf.area();

        for (int compIdx = 0; compIdx < numComponents-1; compIdx++)
        {
            componentFlux[compIdx] = reducedFlux[compIdx];
            componentFlux[numComponents-1] -= reducedFlux[compIdx];
        }
        return componentFlux ;
    }

private:

    static ReducedComponentMatrix setupMSMatrix_(const Problem& problem,
                                                 const Element& element,
                                                 const FVElementGeometry& fvGeometry,
                                                 const ElementVolumeVariables& elemVolVars,
                                                 const SubControlVolumeFace& scvf,
                                                 const int phaseIdx,
                                                 const ComponentFluxVector moleFrac)
    {
        ReducedComponentMatrix reducedDiffusionMatrix(0.0);

        const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
        const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];
        const auto insideScvIdx = scvf.insideScvIdx();
        const auto outsideScvIdx = scvf.outsideScvIdx();

        //this is to not devide by 0 if the saturation in 0 and the effectiveDiffusivity becomes zero due to that
        if(insideVolVars.saturation(phaseIdx) == 0 || outsideVolVars.saturation(phaseIdx) == 0)
            return reducedDiffusionMatrix;

        for (int compIIdx = 0; compIIdx < numComponents-1; compIIdx++)
        {
            // effective diffusion tensors
            using EffDiffModel = typename GET_PROP_TYPE(TypeTag, EffectiveDiffusivityModel);

            const auto xi = moleFrac[compIIdx];

            //calculate diffusivity for i,numComponents
            auto tinInside = getDiffusionCoefficient(phaseIdx, compIIdx, numComponents-1, problem, element, insideVolVars, fvGeometry.scv(insideScvIdx));
            tinInside = EffDiffModel::effectiveDiffusivity(insideVolVars.porosity(), insideVolVars.saturation(phaseIdx), tinInside);
            auto tinOutside = getDiffusionCoefficient(phaseIdx, compIIdx, numComponents-1, problem, element, outsideVolVars, fvGeometry.scv(outsideScvIdx));
            tinOutside = EffDiffModel::effectiveDiffusivity(insideVolVars.porosity(), insideVolVars.saturation(phaseIdx), tinOutside);

                            // scale by extrusion factor
            tinInside *= insideVolVars.extrusionFactor();
            tinOutside *= outsideVolVars.extrusionFactor();

            // the resulting averaged diffusion tensor
            const auto tin = problem.spatialParams().harmonicMean(tinInside, tinOutside, scvf.unitOuterNormal());

            //begin the entrys of the diffusion matrix of the diagonal
            reducedDiffusionMatrix[compIIdx][compIIdx] += xi/tin;

            // now set the rest of the entries (off-diagonal and additional entries for diagonal)
            for (int compJIdx = 0; compJIdx < numComponents; compJIdx++)
            {
                //we don't want to calculate e.g. water in water diffusion
                if (compIIdx == compJIdx)
                    continue;

                //calculate diffusivity for compIIdx, compJIdx
                const auto xj = moleFrac[compJIdx];
                auto tijInside = getDiffusionCoefficient(phaseIdx, compIIdx, compJIdx, problem, element, insideVolVars, fvGeometry.scv(insideScvIdx));
                tijInside = EffDiffModel::effectiveDiffusivity(insideVolVars.porosity(), insideVolVars.saturation(phaseIdx), tijInside);
                auto tijOutside = getDiffusionCoefficient(phaseIdx, compIIdx, compJIdx, problem, element, outsideVolVars, fvGeometry.scv(outsideScvIdx));
                tijOutside = EffDiffModel::effectiveDiffusivity(insideVolVars.porosity(), outsideVolVars.saturation(phaseIdx), tijOutside);

                // scale by extrusion factor
                tijInside *= insideVolVars.extrusionFactor();
                tijOutside *= outsideVolVars.extrusionFactor();

                // the resulting averaged diffusion tensor
                const auto tij = problem.spatialParams().harmonicMean(tijInside, tijOutside, scvf.unitOuterNormal());

                reducedDiffusionMatrix[compIIdx][compIIdx] += xj/tij;
                if (compJIdx < numComponents-1)
                    reducedDiffusionMatrix[compIIdx][compJIdx] +=xi*(1/tin - 1/tij);
            }
        }
        return reducedDiffusionMatrix;
    }

private:
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
