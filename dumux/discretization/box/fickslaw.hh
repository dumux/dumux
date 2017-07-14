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
#ifndef DUMUX_DISCRETIZATION_BOX_FICKS_LAW_HH
#define DUMUX_DISCRETIZATION_BOX_FICKS_LAW_HH

#include <dune/common/float_cmp.hh>

#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>

#include <dumux/implicit/properties.hh>
#include <dumux/discretization/methods.hh>

namespace Dumux
{

namespace Properties
{
// forward declaration of properties
NEW_PROP_TAG(NumPhases);
NEW_PROP_TAG(FluidSystem);
NEW_PROP_TAG(EffectiveDiffusivityModel);
}

/*!
 * \ingroup BoxFicksLaw
 * \brief Specialization of Fick's Law for the box method.
 */
template <class TypeTag>
class FicksLawImplementation<TypeTag, DiscretizationMethods::Box>
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using EffDiffModel = typename GET_PROP_TYPE(TypeTag, EffectiveDiffusivityModel);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using ElementFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using IndexType = typename GridView::IndexSet::IndexType;
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using Element = typename GridView::template Codim<0>::Entity;

    enum { dim = GridView::dimension} ;
    enum { dimWorld = GridView::dimensionworld} ;
    enum
    {
        numPhases = GET_PROP_VALUE(TypeTag, NumPhases),
        numComponents = GET_PROP_VALUE(TypeTag,NumComponents)
    };
    using DimWorldMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
    using ComponentFluxVector = Dune::FieldVector<Scalar, numComponents>;

public:

    static ComponentFluxVector flux(const Problem& problem,
                                    const Element& element,
                                    const FVElementGeometry& fvGeometry,
                                    const ElementVolumeVariables& elemVolVars,
                                    const SubControlVolumeFace& scvf,
                                    const int phaseIdx,
                                    const ElementFluxVariablesCache& elemFluxVarsCache)
    {
        ComponentFluxVector componentFlux(0.0);
        for (int compIdx = 0; compIdx < numComponents; compIdx++)
        {
            // get inside and outside diffusion tensors and calculate the harmonic mean
            const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
            const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];

            // effective diffusion tensors
            auto insideD = EffDiffModel::effectiveDiffusivity(insideVolVars.porosity(),
                                                            insideVolVars.saturation(phaseIdx),
                                                            insideVolVars.diffusionCoefficient(phaseIdx, compIdx));
            auto outsideD = EffDiffModel::effectiveDiffusivity(outsideVolVars.porosity(),
                                                            outsideVolVars.saturation(phaseIdx),
                                                            outsideVolVars.diffusionCoefficient(phaseIdx, compIdx));

            // scale by extrusion factor
            insideD *= insideVolVars.extrusionFactor();
            outsideD *= outsideVolVars.extrusionFactor();

            // the resulting averaged diffusion tensor
            const auto D = problem.spatialParams().harmonicMean(insideD, outsideD, scvf.unitOuterNormal());

            // evaluate gradX at integration point and interpolate density
            const auto& fluxVarsCache = elemFluxVarsCache[scvf];
            const auto& jacInvT = fluxVarsCache.jacInvT();
            const auto& shapeJacobian = fluxVarsCache.shapeJacobian();
            const auto& shapeValues = fluxVarsCache.shapeValues();

            GlobalPosition gradX(0.0);
            Scalar rho(0.0);
            for (auto&& scv : scvs(fvGeometry))
            {
                const auto& volVars = elemVolVars[scv];

                // density interpolation
                rho +=  volVars.molarDensity(phaseIdx)*shapeValues[scv.indexInElement()][0];

                // the mole/mass fraction gradient
                GlobalPosition gradN;
                jacInvT.mv(shapeJacobian[scv.indexInElement()][0], gradN);
                gradX.axpy(volVars.moleFraction(phaseIdx, compIdx), gradN);
            }

            // apply the diffusion tensor and return the flux
            auto DGradX = applyDiffusionTensor_(D, gradX);
            componentFlux[compIdx] = -1.0*rho*(DGradX*scvf.unitOuterNormal())*scvf.area();
        }
        return componentFlux;
    }

private:
    static GlobalPosition applyDiffusionTensor_(const DimWorldMatrix& D, const GlobalPosition& gradI)
    {
        GlobalPosition result(0.0);
        D.mv(gradI, result);
        return result;
    }

    static GlobalPosition applyDiffusionTensor_(const Scalar d, const GlobalPosition& gradI)
    {
        GlobalPosition result(gradI);
        result *= d;
        return result;
    }
};

} // end namespace Dumux

#endif
