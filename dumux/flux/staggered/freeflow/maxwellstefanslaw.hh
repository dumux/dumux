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
 * \ingroup StaggeredFlux
 * \brief Specialization of Maxwell Stefan's Law for the Staggered method.
 */
#ifndef DUMUX_DISCRETIZATION_STAGGERED_MAXWELL_STEFAN_LAW_HH
#define DUMUX_DISCRETIZATION_STAGGERED_MAXWELL_STEFAN_LAW_HH

#include <dune/common/fmatrix.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/extrusion.hh>

#include <dumux/flux/fluxvariablescaching.hh>
#include <dumux/flux/referencesystemformulation.hh>
#include <dumux/flux/maxwellstefandiffusioncoefficients.hh>

namespace Dumux {

// forward declaration
template <class TypeTag, DiscretizationMethod discMethod, ReferenceSystemFormulation referenceSystem>
class MaxwellStefansLawImplementation;

/*!
 * \ingroup StaggeredFlux
 * \brief Specialization of Maxwell Stefan's Law for the Staggered method.
 */
template <class TypeTag, ReferenceSystemFormulation referenceSystem>
class MaxwellStefansLawImplementation<TypeTag, DiscretizationMethod::staggered, referenceSystem>
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using Extrusion = Extrusion_t<GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using CellCenterPrimaryVariables = GetPropType<TypeTag, Properties::CellCenterPrimaryVariables>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;

    static const int numComponents = ModelTraits::numFluidComponents();
    static const int numPhases = ModelTraits::numFluidPhases();
    static constexpr bool useMoles = getPropValue<TypeTag, Properties::UseMoles>();

    using ReducedComponentVector = Dune::FieldVector<Scalar, numComponents-1>;
    using ReducedComponentMatrix = Dune::FieldMatrix<Scalar, numComponents-1, numComponents-1>;

    static_assert(ModelTraits::numFluidPhases() == 1, "Only one phase allowed supported!");

    static_assert(referenceSystem == ReferenceSystemFormulation::massAveraged, "only the mass averaged reference system is supported for the Maxwell-Stefan formulation");

public:
    // state the discretization method this implementation belongs to
    static const DiscretizationMethod discMethod = DiscretizationMethod::staggered;
    //return the reference system
    static constexpr ReferenceSystemFormulation referenceSystemFormulation()
    { return referenceSystem; }

    //! state the type for the corresponding cache and its filler
    //! We don't cache anything for this law
    using Cache = FluxVariablesCaching::EmptyDiffusionCache;
    using CacheFiller = FluxVariablesCaching::EmptyCacheFiller;

    using DiffusionCoefficientsContainer = MaxwellStefanDiffusionCoefficients<Scalar, numPhases, numComponents>;

    /*!
     * \brief Returns the diffusive fluxes of all components within
     *        a fluid phase across the given sub-control volume face.
     *        The computed fluxes are given in kg/s.
     */
    template<class ElementVolumeVariables>
    static CellCenterPrimaryVariables flux(const Problem& problem,
                                           const Element& element,
                                           const FVElementGeometry& fvGeometry,
                                           const ElementVolumeVariables& elemVolVars,
                                           const SubControlVolumeFace& scvf)
    {
        //this is to calculate the maxwellStefan diffusion in a multicomponent system.
        //see: Multicomponent Mass Transfer. R. Taylor u. R. Krishna. J. Wiley & Sons, New York 1993
        CellCenterPrimaryVariables componentFlux(0.0);
        ReducedComponentVector moleFracInside(0.0);
        ReducedComponentVector moleFracOutside(0.0);
        ReducedComponentVector reducedFlux(0.0);
        ReducedComponentMatrix reducedDiffusionMatrixInside(0.0);
        ReducedComponentMatrix reducedDiffusionMatrixOutside(0.0);

        // get inside/outside volume variables
        const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
        const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];
        const auto rhoInside = insideVolVars.density();
        const auto rhoOutside = outsideVolVars.density();

        //to implement outflow boundaries correctly we need to loop over all components but the main component as only for the transported ones we implement the outflow boundary. diffusion then is 0.
        for(int compIdx = 0; compIdx < numComponents; ++compIdx)
        {
            if(compIdx == FluidSystem::getMainComponent(0))
                continue;

            // get equation index
            const auto eqIdx = Indices::conti0EqIdx + compIdx;
            if(scvf.boundary())
            {
               const auto bcTypes = problem.boundaryTypes(element, scvf);
                 if((bcTypes.isOutflow(eqIdx)) || (bcTypes.isSymmetry()))
                    return componentFlux;
            }
        }

        //calculate the mole fraction vectors
        for (int compIdx = 0; compIdx < numComponents-1; compIdx++)
        {
            //calculate x_inside
            const auto xInside = insideVolVars.moleFraction(compIdx);
            //calculate outside molefraction with the respective transmissibility
            const auto xOutside = outsideVolVars.moleFraction(compIdx);

            moleFracInside[compIdx] = xInside;
            moleFracOutside[compIdx] = xOutside;
        }

        //now we have to do the tpfa: J_i = J_j which leads to: tij(xi -xj) = -rho Bi^-1 omegai(x*-xi) with x* = (omegai Bi^-1 + omegaj Bj^-1)^-1 (xi omegai Bi^-1 + xj omegaj Bj^-1) with i inside and j outside
        reducedDiffusionMatrixInside = setupMSMatrix_(problem, fvGeometry, insideVolVars, scvf);

        reducedDiffusionMatrixOutside = setupMSMatrix_(problem, fvGeometry, outsideVolVars, scvf);

        const auto insideScvIdx = scvf.insideScvIdx();
        const auto& insideScv = fvGeometry.scv(insideScvIdx);
        const auto outsideScvIdx = scvf.outsideScvIdx();
        const auto& outsideScv = fvGeometry.scv(outsideScvIdx);

        const Scalar insideDistance = (insideScv.dofPosition() - scvf.ipGlobal()).two_norm();

        const Scalar omegai = calculateOmega_(insideDistance, insideVolVars.extrusionFactor());

        //if on boundary
        if (scvf.boundary())
        {
            moleFracOutside -= moleFracInside;
            reducedDiffusionMatrixInside.solve(reducedFlux, moleFracOutside);
            reducedFlux *= omegai*rhoInside;
        }
        else
        {
            const Scalar outsideDistance = (outsideScv.dofPosition() - scvf.ipGlobal()).two_norm();
            const Scalar omegaj = calculateOmega_(outsideDistance, outsideVolVars.extrusionFactor());

            reducedDiffusionMatrixInside.invert();
            reducedDiffusionMatrixOutside.invert();
            reducedDiffusionMatrixInside *= omegai*rhoInside;
            reducedDiffusionMatrixOutside *= omegaj*rhoOutside;

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

        reducedFlux *= -Extrusion::area(scvf);

        for (int compIdx = 0; compIdx < numComponents-1; compIdx++)
        {
            componentFlux[compIdx] = reducedFlux[compIdx];
            componentFlux[numComponents-1] -= reducedFlux[compIdx];
        }

        return componentFlux ;
    }

private:
    static Scalar calculateOmega_(const Scalar distance,
                                  const Scalar extrusionFactor)
    {
        Scalar omega = 1/distance;
        omega *= extrusionFactor;

        return omega;
    }

    static ReducedComponentMatrix setupMSMatrix_(const Problem& problem,
                                                const FVElementGeometry& fvGeometry,
                                                const VolumeVariables& volVars,
                                                const SubControlVolumeFace& scvf)
    {
        ReducedComponentMatrix reducedDiffusionMatrix(0.0);

        for (int compIIdx = 0; compIIdx < numComponents-1; compIIdx++)
        {
            const auto xi = volVars.moleFraction(compIIdx);
            const auto avgMolarMass = volVars.averageMolarMass(0);
            const auto Mn = FluidSystem::molarMass(numComponents-1);
            const Scalar tin = getEffectiveDiffusionCoefficient_(volVars, compIIdx, numComponents-1);

            // set the entries of the diffusion matrix of the diagonal
            reducedDiffusionMatrix[compIIdx][compIIdx] +=  xi*avgMolarMass/(tin*Mn);

            for (int compJIdx = 0; compJIdx < numComponents; compJIdx++)
            {
                // we don't want to calculate e.g. water in water diffusion
                if (compJIdx == compIIdx)
                    continue;

                const auto xj = volVars.moleFraction(compJIdx);
                const auto Mi = FluidSystem::molarMass(compIIdx);
                const auto Mj = FluidSystem::molarMass(compJIdx);
                const Scalar tij = getEffectiveDiffusionCoefficient_(volVars, compIIdx, compJIdx);
                reducedDiffusionMatrix[compIIdx][compIIdx] +=  xj*avgMolarMass/(tij*Mi);
                if (compJIdx < numComponents-1)
                    reducedDiffusionMatrix[compIIdx][compJIdx] += xi*(avgMolarMass/(tin*Mn) - avgMolarMass/(tij*Mj));
            }
        }
        return reducedDiffusionMatrix;
    }

    static Scalar getEffectiveDiffusionCoefficient_(const VolumeVariables& volVars, const int phaseIdx, const int compIdx)
    {
        return volVars.effectiveDiffusionCoefficient(phaseIdx, FluidSystem::getMainComponent(phaseIdx), compIdx);
    }
};
} // end namespace

#endif
