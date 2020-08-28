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
 * \brief This file contains the data which is required to calculate
 *        diffusive molar fluxes due to molecular diffusion with Maxwell Stefan
 */
#ifndef DUMUX_DISCRETIZATION_BOX_MAXWELL_STEFAN_LAW_HH
#define DUMUX_DISCRETIZATION_BOX_MAXWELL_STEFAN_LAW_HH

#include <dune/common/float_cmp.hh>
#include <dune/common/fmatrix.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/extrusion.hh>

#include <dumux/flux/maxwellstefandiffusioncoefficients.hh>
#include <dumux/flux/fluxvariablescaching.hh>
#include <dumux/flux/referencesystemformulation.hh>

namespace Dumux {

// forward declaration
template <class TypeTag, DiscretizationMethod discMethod, ReferenceSystemFormulation referenceSystem>
class MaxwellStefansLawImplementation;

/*!
 * \ingroup BoxFlux
 * \brief Specialization of Maxwell Stefan's Law for the Box method.
 */
template <class TypeTag, ReferenceSystemFormulation referenceSystem>
class MaxwellStefansLawImplementation<TypeTag, DiscretizationMethod::box, referenceSystem>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using Extrusion = Extrusion_t<GridGeometry>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using ElementFluxVariablesCache = typename GetPropType<TypeTag, Properties::GridFluxVariablesCache>::LocalView;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

    static constexpr auto numFluidPhases = GetPropType<TypeTag, Properties::ModelTraits>::numFluidPhases();
    static constexpr auto numComponents = GetPropType<TypeTag, Properties::ModelTraits>::numFluidComponents();

    using ComponentFluxVector = Dune::FieldVector<Scalar, numComponents>;
    using ReducedComponentVector = Dune::FieldVector<Scalar, numComponents-1>;
    using ReducedComponentMatrix = Dune::FieldMatrix<Scalar, numComponents-1, numComponents-1>;

    static_assert(referenceSystem == ReferenceSystemFormulation::massAveraged, "only the mass averaged reference system is supported for the Maxwell-Stefan formulation");

public:
    using DiffusionCoefficientsContainer = MaxwellStefanDiffusionCoefficients<Scalar, numFluidPhases, numComponents>;

    //return the reference system
    static constexpr ReferenceSystemFormulation referenceSystemFormulation()
    { return referenceSystem; }

    /*!
     * \brief Returns the diffusive fluxes of all components within
     *        a fluid phase across the given sub-control volume face.
     *        The computed fluxes are given in kg/s.
     */
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
        Scalar avgMolarMass(0.0);
        for (auto&& scv : scvs(fvGeometry))
        {
            const auto& volVars = elemVolVars[scv];

            // density interpolation
            rho +=  volVars.density(phaseIdx)*shapeValues[scv.indexInElement()][0];
            //average molar mass interpolation
            avgMolarMass += volVars.averageMolarMass(phaseIdx)*shapeValues[scv.indexInElement()][0];
            //interpolate the mole fraction for the diffusion matrix
            for (int compIdx = 0; compIdx < numComponents; compIdx++)
            {
              moleFrac[compIdx] += volVars.moleFraction(phaseIdx, compIdx)*shapeValues[scv.indexInElement()][0];
            }
        }

        reducedDiffusionMatrix = setupMSMatrix_(problem, element, fvGeometry, elemVolVars, scvf, phaseIdx, moleFrac, avgMolarMass);

        for (int compIdx = 0; compIdx < numComponents-1; compIdx++)
        {
            Dune::FieldVector<Scalar, GridView::dimensionworld> gradX(0.0);
            for (auto&& scv : scvs(fvGeometry))
            {
                const auto& volVars = elemVolVars[scv];

                // the mole/mass fraction gradient
                gradX.axpy(volVars.moleFraction(phaseIdx, compIdx), fluxVarsCache.gradN(scv.indexInElement()));
            }

           normalX[compIdx] = gradX *scvf.unitOuterNormal();
        }
         reducedDiffusionMatrix.solve(reducedFlux,normalX);
         reducedFlux *= -1.0*rho*Extrusion::area(scvf);

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
                                                 const ComponentFluxVector moleFrac,
                                                 const Scalar avgMolarMass)
    {
        ReducedComponentMatrix reducedDiffusionMatrix(0.0);

        const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
        const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];

        //this is to not devide by 0 if the saturation in 0 and the effectiveDiffusionCoefficient becomes zero due to that
        if(Dune::FloatCmp::eq<Scalar>(insideVolVars.saturation(phaseIdx), 0) || Dune::FloatCmp::eq<Scalar>(outsideVolVars.saturation(phaseIdx), 0))
            return reducedDiffusionMatrix;

        for (int compIIdx = 0; compIIdx < numComponents-1; compIIdx++)
        {
            //calculate diffusivity for i,numComponents
            const auto xi = moleFrac[compIIdx];
            const auto Mn = FluidSystem::molarMass(numComponents-1);

            auto tinInside = insideVolVars.effectiveDiffusionCoefficient(phaseIdx, compIIdx, numComponents-1);
            auto tinOutside = outsideVolVars.effectiveDiffusionCoefficient(phaseIdx, compIIdx, numComponents-1);

            // scale by extrusion factor
            tinInside *= insideVolVars.extrusionFactor();
            tinOutside *= outsideVolVars.extrusionFactor();

            // the resulting averaged diffusion tensor
            const auto tin = problem.spatialParams().harmonicMean(tinInside, tinOutside, scvf.unitOuterNormal());

            //begin the entrys of the diffusion matrix of the diagonal
            reducedDiffusionMatrix[compIIdx][compIIdx] += xi*avgMolarMass/(tin*Mn);

            // now set the rest of the entries (off-diagonal and additional entries for diagonal)
            for (int compJIdx = 0; compJIdx < numComponents; compJIdx++)
            {
                //we don't want to calculate e.g. water in water diffusion
                if (compIIdx == compJIdx)
                    continue;

                //calculate diffusivity for compIIdx, compJIdx
                const auto xj = moleFrac[compJIdx];
                const auto Mi = FluidSystem::molarMass(compIIdx);
                const auto Mj = FluidSystem::molarMass(compJIdx);

                // effective diffusion tensors
                auto tijInside = insideVolVars.effectiveDiffusionCoefficient(phaseIdx, compIIdx, compJIdx);
                auto tijOutside = outsideVolVars.effectiveDiffusionCoefficient(phaseIdx, compIIdx, compJIdx);

                // scale by extrusion factor
                tijInside *= insideVolVars.extrusionFactor();
                tijOutside *= outsideVolVars.extrusionFactor();

                // the resulting averaged diffusion tensor
                const auto tij = problem.spatialParams().harmonicMean(tijInside, tijOutside, scvf.unitOuterNormal());

                reducedDiffusionMatrix[compIIdx][compIIdx] += xj*avgMolarMass/(tij*Mi);
                if (compJIdx < numComponents-1)
                    reducedDiffusionMatrix[compIIdx][compJIdx] +=xi*(avgMolarMass/(tin*Mn) - avgMolarMass/(tij*Mj));
            }
        }
        return reducedDiffusionMatrix;
    }

};
} // end namespace Dumux

#endif
