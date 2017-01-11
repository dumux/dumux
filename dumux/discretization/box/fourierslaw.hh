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
 *        energy fluxes due to molecular diffusion with Fourier's law.
 */
#ifndef DUMUX_DISCRETIZATION_BOX_FOURIERS_LAW_HH
#define DUMUX_DISCRETIZATION_BOX_FOURIERS_LAW_HH

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
NEW_PROP_TAG(FluidState);
NEW_PROP_TAG(FluidSystem);
NEW_PROP_TAG(ThermalConductivityModel);
}

/*!
 * \ingroup BoxFouriersLaw
 * \brief Specialization of Fourier's Law for the box method.
 */
template <class TypeTag>
class FouriersLawImplementation<TypeTag, DiscretizationMethods::Box>
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using FluidState = typename GET_PROP_TYPE(TypeTag, FluidState);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using ThermalConductivityModel = typename GET_PROP_TYPE(TypeTag, ThermalConductivityModel);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using ElementFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using IndexType = typename GridView::IndexSet::IndexType;
    using Stencil = typename std::vector<IndexType>;

    using Element = typename GridView::template Codim<0>::Entity;

    enum { dim = GridView::dimension} ;
    enum { dimWorld = GridView::dimensionworld} ;
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases)} ;

    using DimWorldMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

public:
    static Scalar flux(const Problem& problem,
                       const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolumeFace& scvf,
                       const ElementFluxVariablesCache& elemFluxVarsCache)
    {
        // get inside and outside diffusion tensors and calculate the harmonic mean
        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& outsideScv = fvGeometry.scv(scvf.outsideScvIdx());
        const auto& insideVolVars = elemVolVars[insideScv];
        const auto& outsideVolVars = elemVolVars[outsideScv];

        // effective diffusion tensors
        auto insideLambda = ThermalConductivityModel::effectiveThermalConductivity(insideVolVars, problem.spatialParams(), element, fvGeometry, insideScv);
        auto outsideLambda = ThermalConductivityModel::effectiveThermalConductivity(outsideVolVars, problem.spatialParams(), element, fvGeometry, outsideScv);

        // scale by extrusion factor
        insideLambda *= insideVolVars.extrusionFactor();
        outsideLambda *= outsideVolVars.extrusionFactor();

        // the resulting averaged diffusion tensor
        const auto lambda = problem.spatialParams().harmonicMean(insideLambda, outsideLambda, scvf.unitOuterNormal());

        // evaluate gradTemp at integration point
        const auto& fluxVarsCache = elemFluxVarsCache[scvf];
        const auto& jacInvT = fluxVarsCache.jacInvT();
        const auto& shapeJacobian = fluxVarsCache.shapeJacobian();

        GlobalPosition gradTemp(0.0);
        for (auto&& scv : scvs(fvGeometry))
        {
            const auto& volVars = elemVolVars[scv];

            // the mole/mass fraction gradient
            GlobalPosition gradI;
            jacInvT.mv(shapeJacobian[scv.index()][0], gradI);
            gradI *= volVars.temperature();
            gradTemp += gradI;
        }

        // apply the diffusion tensor and return the flux
        auto lambdaGradT = applyHeatConductivityTensor_(lambda, gradTemp);
        return -1.0*(lambdaGradT*scvf.unitOuterNormal())*scvf.area();
    }

    // This is for compatibility with the cc methods. The flux stencil info is obsolete for the box method.
    static Stencil stencil(const Problem& problem,
                           const Element& element,
                           const FVElementGeometry& fvGeometry,
                           const SubControlVolumeFace& scvf)
    { return Stencil(0); }

private:
    static GlobalPosition applyHeatConductivityTensor_(const DimWorldMatrix& Lambda, const GlobalPosition& gradI)
    {
        GlobalPosition result(0.0);
        Lambda.mv(gradI, result);
        return result;
    }

    static GlobalPosition applyHeatConductivityTensor_(const Scalar lambda, const GlobalPosition& gradI)
    {
        GlobalPosition result(gradI);
        result *= lambda;
        return result;
    }
};

} // end namespace Dumux

#endif
