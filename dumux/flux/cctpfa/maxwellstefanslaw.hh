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
 * \ingroup CCTpfaFlux
 * \brief This file contains the data which is required to calculate
 *        diffusive mass fluxes due to molecular diffusion with Maxwell Stefan
 */
#ifndef DUMUX_DISCRETIZATION_CC_TPFA_MAXWELL_STEFAN_LAW_HH
#define DUMUX_DISCRETIZATION_CC_TPFA_MAXWELL_STEFAN_LAW_HH

#include <dune/common/fmatrix.hh>
#include <dune/common/float_cmp.hh>

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
 * \ingroup CCTpfaFlux
 * \brief Specialization of Maxwell Stefan's Law for the CCTpfa method.
 */
template <class TypeTag, ReferenceSystemFormulation referenceSystem>
class MaxwellStefansLawImplementation<TypeTag, DiscretizationMethod::cctpfa, referenceSystem >
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using Extrusion = Extrusion_t<GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementFluxVariablesCache = typename GetPropType<TypeTag, Properties::GridFluxVariablesCache>::LocalView;
    using FluxVariablesCache = GetPropType<TypeTag, Properties::FluxVariablesCache>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;

    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;
    static const int numPhases = GetPropType<TypeTag, Properties::ModelTraits>::numFluidPhases();
    static const int numComponents = GetPropType<TypeTag, Properties::ModelTraits>::numFluidComponents();

    static_assert(referenceSystem == ReferenceSystemFormulation::massAveraged, "only the mass averaged reference system is supported for the Maxwell-Stefan formulation");


    using ComponentFluxVector = Dune::FieldVector<Scalar, numComponents>;
    using ReducedComponentVector = Dune::FieldVector<Scalar, numComponents-1>;
    using ReducedComponentMatrix = Dune::FieldMatrix<Scalar, numComponents-1, numComponents-1>;

public:
    // state the discretization method this implementation belongs to
    static const DiscretizationMethod discMethod = DiscretizationMethod::cctpfa;
    //return the reference system
    static constexpr ReferenceSystemFormulation referenceSystemFormulation()
    { return referenceSystem; }

    using DiffusionCoefficientsContainer = MaxwellStefanDiffusionCoefficients<Scalar, numPhases, numComponents>;

    //! state the type for the corresponding cache and its filler
    //! We don't cache anything for this law
    using Cache = FluxVariablesCaching::EmptyDiffusionCache;
    using CacheFiller = FluxVariablesCaching::EmptyCacheFiller;

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
        ReducedComponentVector moleFracInside(0.0);
        ReducedComponentVector moleFracOutside(0.0);
        ReducedComponentVector reducedFlux(0.0);
        ReducedComponentMatrix reducedDiffusionMatrixInside(0.0);
        ReducedComponentMatrix reducedDiffusionMatrixOutside(0.0);

        // get inside/outside volume variables
        const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
        const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];
        const auto rhoInside = insideVolVars.density(phaseIdx);
        const auto rhoOutside = outsideVolVars.density(phaseIdx);
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
        if(!(Dune::FloatCmp::eq<Scalar>(insideVolVars.saturation(phaseIdx), 0) || Dune::FloatCmp::eq<Scalar>(outsideVolVars.saturation(phaseIdx), 0)))
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
                reducedFlux *= omegai*rhoInside;
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

        //this is to not devide by 0 if the saturation in 0 and the effectiveDiffusionCoefficient becomes zero due to that
        if(Dune::FloatCmp::eq<Scalar>(volVars.saturation(phaseIdx), 0))
            return reducedDiffusionMatrix;

        const auto avgMolarMass = volVars.averageMolarMass(phaseIdx);

        for (int compIIdx = 0; compIIdx < numComponents-1; compIIdx++)
        {
            const auto xi = volVars.moleFraction(phaseIdx, compIIdx);

            const auto Mn = FluidSystem::molarMass(numComponents-1);
            Scalar tin = volVars.effectiveDiffusionCoefficient(phaseIdx, compIIdx, numComponents-1);

            // set the entries of the diffusion matrix of the diagonal
            reducedDiffusionMatrix[compIIdx][compIIdx] += xi*avgMolarMass/(tin*Mn);

            // now set the rest of the entries (off-diagonal and additional entries for diagonal)
            for (int compJIdx = 0; compJIdx < numComponents; compJIdx++)
            {
                // we don't want to calculate e.g. water in water diffusion
                if (compJIdx == compIIdx)
                    continue;

                const auto xj = volVars.moleFraction(phaseIdx, compJIdx);
                const auto Mi = FluidSystem::molarMass(compIIdx);
                const auto Mj = FluidSystem::molarMass(compJIdx);
                Scalar tij = volVars.effectiveDiffusionCoefficient(phaseIdx, compIIdx, compJIdx);
                reducedDiffusionMatrix[compIIdx][compIIdx] += xj*avgMolarMass/(tij*Mi);

                if (compJIdx < numComponents-1)
                    reducedDiffusionMatrix[compIIdx][compJIdx] += xi*(avgMolarMass/(tin*Mn) - avgMolarMass/(tij*Mj));
            }
        }
        return reducedDiffusionMatrix;
    }

};
} // end namespace

#endif
