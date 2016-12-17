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
#ifndef DUMUX_DISCRETIZATION_CC_TPFA_FICKS_LAW_HH
#define DUMUX_DISCRETIZATION_CC_TPFA_FICKS_LAW_HH

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
 * \ingroup CCTpfaFicksLaw
 * \brief Specialization of Fick's Law for the CCTpfa method.
 */
template <class TypeTag>
class FicksLawImplementation<TypeTag, DiscretizationMethods::CCTpfa >
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using EffDiffModel = typename GET_PROP_TYPE(TypeTag, EffectiveDiffusivityModel);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using IndexType = typename GridView::IndexSet::IndexType;
    using Stencil = typename std::vector<IndexType>;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);

    enum { dim = GridView::dimension} ;
    enum { dimWorld = GridView::dimensionworld} ;
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases)} ;

    using DimWorldMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

public:
    // state the discretization method this implementation belongs to
    static const DiscretizationMethods myDiscretizationMethod = DiscretizationMethods::CCTpfa;

    static Scalar flux(const Problem& problem,
                       const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolumeFace& scvf,
                       int phaseIdx, int compIdx,
                       const ElementFluxVariablesCache& elemFluxVarsCache,
                       bool useMoles = true)
    {
        // diffusion tensors are always solution dependent
        Scalar tij = calculateTransmissibility_(problem, element, fvGeometry, elemVolVars, scvf, phaseIdx, compIdx);

        // Get the inside volume variables
        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& insideVolVars = elemVolVars[insideScv];

        // and the outside volume variables
        const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];

        // compute the diffusive flux using mole fractions
        if (useMoles)
        {
            const auto xInside = insideVolVars.moleFraction(phaseIdx, compIdx);
            const auto xOutside = outsideVolVars.moleFraction(phaseIdx, compIdx);
            const auto rho = 0.5*(insideVolVars.molarDensity(phaseIdx) + outsideVolVars.molarDensity(phaseIdx));

            return rho*tij*(xInside - xOutside);
        }
        // compute the diffusive flux using mass fractions
        else
        {
            const auto xInside = insideVolVars.massFraction(phaseIdx, compIdx);
            const auto xOutside = outsideVolVars.massFraction(phaseIdx, compIdx);
            const auto rho = 0.5*(insideVolVars.density(phaseIdx) + outsideVolVars.density(phaseIdx));

            return rho*tij*(xInside - xOutside);
        }
    }

    static Stencil stencil(const Problem& problem,
                           const Element& element,
                           const FVElementGeometry& fvGeometry,
                           const SubControlVolumeFace& scvf)
    {
        if (!scvf.boundary())
            return Stencil({scvf.insideScvIdx(), scvf.outsideScvIdx()});
        else
            return Stencil({scvf.insideScvIdx()});
    }

private:


    static Scalar calculateTransmissibility_(const Problem& problem,
                                             const Element& element,
                                             const FVElementGeometry& fvGeometry,
                                             const ElementVolumeVariables& elemVolVars,
                                             const SubControlVolumeFace& scvf,
                                             const int phaseIdx, const int compIdx)
    {
        Scalar tij;

        const auto insideScvIdx = scvf.insideScvIdx();
        const auto& insideScv = fvGeometry.scv(insideScvIdx);
        const auto& insideVolVars = elemVolVars[insideScvIdx];

        auto insideD = insideVolVars.diffusionCoefficient(phaseIdx, compIdx);
        insideD = EffDiffModel::effectiveDiffusivity(insideVolVars.porosity(), insideVolVars.saturation(phaseIdx), insideD);
        Scalar ti = calculateOmega_(problem, element, scvf, insideD, insideScv);

        if (!scvf.boundary())
        {
            const auto outsideScvIdx = scvf.outsideScvIdx();
            const auto& outsideScv = fvGeometry.scv(outsideScvIdx);
            const auto& outsideVolVars = elemVolVars[outsideScvIdx];

            auto outsideD = outsideVolVars.diffusionCoefficient(phaseIdx, compIdx);
            outsideD = EffDiffModel::effectiveDiffusivity(outsideVolVars.porosity(), outsideVolVars.saturation(phaseIdx), outsideD);
            Scalar tj = -1.0*calculateOmega_(problem, element, scvf, outsideD, outsideScv);

            // check if we are dividing by zero!
            if (ti*tj <= 0.0)
                tij = 0;
            else
                tij = scvf.area()*(ti * tj)/(ti + tj);
        }
        else
        {
            tij = scvf.area()*ti;
        }

        return tij;
    }

    static Scalar calculateOmega_(const Problem& problem,
                                  const Element& element,
                                  const SubControlVolumeFace& scvf,
                                  const DimWorldMatrix &D,
                                  const SubControlVolume &scv)
    {
        GlobalPosition Dnormal;
        D.mv(scvf.unitOuterNormal(), Dnormal);

        auto distanceVector = scvf.center();
        distanceVector -= scv.center();
        distanceVector /= distanceVector.two_norm2();

        Scalar omega = Dnormal * distanceVector;
        omega *= problem.boxExtrusionFactor(element, scv);

        return omega;
    }

    static Scalar calculateOmega_(const Problem& problem,
                                  const Element& element,
                                  const SubControlVolumeFace& scvf,
                                  Scalar D,
                                  const SubControlVolume &scv)
    {
        auto distanceVector = scvf.center();
        distanceVector -= scv.center();
        distanceVector /= distanceVector.two_norm2();

        Scalar omega = D * (distanceVector * scvf.unitOuterNormal());
        omega *= problem.boxExtrusionFactor(element, scv);

        return omega;
    }
};
} // end namespace

#endif
