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
 *        heat conduction fluxes with Fourier's law.
 */
#ifndef DUMUX_DISCRETIZATION_CC_TPFA_FOURIERS_LAW_HH
#define DUMUX_DISCRETIZATION_CC_TPFA_FOURIERS_LAW_HH

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
NEW_PROP_TAG(ThermalConductivityModel);
}

/*!
 * \ingroup FouriersLaw
 * \brief Specialization of Fourier's Law for the CCTpfa method.
 */
template <class TypeTag>
class FouriersLawImplementation<TypeTag, DiscretizationMethods::CCTpfa>
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using FluidSystem = typename GET_PROP_TYPE(TypeTag, FluidSystem);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using IndexType = typename GridView::IndexSet::IndexType;
    using Stencil = typename std::vector<IndexType>;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementFluxVarsCache = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);

    enum { dim = GridView::dimension} ;
    enum { dimWorld = GridView::dimensionworld} ;
    enum { numPhases = GET_PROP_VALUE(TypeTag, NumPhases)} ;

    using DimWorldMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

    using ThermalConductivityModel = typename GET_PROP_TYPE(TypeTag, ThermalConductivityModel);

public:

    static Scalar flux(const Problem& problem,
                       const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolumeFace& scvFace,
                       const ElementFluxVarsCache& elemFluxVarsCache)
    {
        // heat conductivities are always solution dependent (?)
        Scalar tij = calculateTransmissibility_(problem, element, fvGeometry, elemVolVars, scvFace);

        // Get the inside volume variables
        const auto& insideScv = fvGeometry.scv(scvFace.insideScvIdx());
        const auto& insideVolVars = elemVolVars[insideScv];

        // and the outside volume variables
        const auto& outsideVolVars = elemVolVars[scvFace.outsideScvIdx()];

        // compute the diffusive flux
        const auto tInside = insideVolVars.temperature();
        const auto tOutside = outsideVolVars.temperature();

        return tij*(tInside - tOutside);
    }

    static Stencil stencil(const Problem& problem,
                           const Element& element,
                           const FVElementGeometry& fvGeometry,
                           const SubControlVolumeFace& scvFace)
    {
        if (!scvFace.boundary())
            return Stencil({scvFace.insideScvIdx(), scvFace.outsideScvIdx()});
        else
            return Stencil({scvFace.insideScvIdx()});
    }

private:

    static Scalar calculateTransmissibility_(const Problem& problem,
                                             const Element& element,
                                             const FVElementGeometry& fvGeometry,
                                             const ElementVolumeVariables& elemVolVars,
                                             const SubControlVolumeFace& scvFace)
    {
        Scalar tij;

        const auto insideScvIdx = scvFace.insideScvIdx();
        const auto& insideScv = fvGeometry.scv(insideScvIdx);
        const auto& insideVolVars = elemVolVars[insideScvIdx];

        auto insideLambda = ThermalConductivityModel::effectiveThermalConductivity(insideVolVars, problem.spatialParams(), element, fvGeometry, insideScv);
        Scalar ti = calculateOmega_(problem, element, scvFace, insideLambda, insideScv);

        if (!scvFace.boundary())
        {
            const auto outsideScvIdx = scvFace.outsideScvIdx();
            const auto& outsideScv = fvGeometry.scv(outsideScvIdx);
            const auto& outsideVolVars = elemVolVars[outsideScvIdx];

            auto outsideLambda = ThermalConductivityModel::effectiveThermalConductivity(outsideVolVars, problem.spatialParams(), element, fvGeometry, outsideScv);
            Scalar tj = -1.0*calculateOmega_(problem, element, scvFace, outsideLambda, outsideScv);

            // check for division by zero!
            if (ti*tj <= 0.0)
                tij = 0;
            else
                tij = scvFace.area()*(ti * tj)/(ti + tj);
        }
        else
        {
            tij = scvFace.area()*ti;
        }

        return tij;
    }

    static Scalar calculateOmega_(const Problem& problem,
                                  const Element& element,
                                  const SubControlVolumeFace& scvFace,
                                  const DimWorldMatrix &lambda,
                                  const SubControlVolume &scv)
    {
        GlobalPosition lambdaNormal;
        lambda.mv(scvFace.unitOuterNormal(), lambdaNormal);

        auto distanceVector = scvFace.center();
        distanceVector -= scv.center();
        distanceVector /= distanceVector.two_norm2();

        Scalar omega = lambdaNormal * distanceVector;
        omega *= problem.boxExtrusionFactor(element, scv);

        return omega;
    }

    static Scalar calculateOmega_(const Problem& problem,
                                  const Element& element,
                                  const SubControlVolumeFace& scvFace,
                                  Scalar lambda,
                                  const SubControlVolume &scv)
    {
        auto distanceVector = scvFace.center();
        distanceVector -= scv.center();
        distanceVector /= distanceVector.two_norm2();

        Scalar omega = lambda * (distanceVector * scvFace.unitOuterNormal());
        omega *= problem.boxExtrusionFactor(element, scv);

        return omega;
    }
};

} // end namespace Dumux

#endif
