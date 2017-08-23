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
#include <dumux/discretization/fluxvariablescaching.hh>

namespace Dumux
{

namespace Properties
{
// forward declaration of properties
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
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using IndexType = typename GridView::IndexSet::IndexType;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using Element = typename GridView::template Codim<0>::Entity;
    using ElementFluxVarsCache = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);
    using FluxVariablesCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);

    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;

    using DimWorldMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

    using ThermalConductivityModel = typename GET_PROP_TYPE(TypeTag, ThermalConductivityModel);

public:
    // state the discretization method this implementation belongs to
    static const DiscretizationMethods myDiscretizationMethod = DiscretizationMethods::CCTpfa;

    //! state the type for the corresponding cache and its filler
    //! We don't cache anything for this law
    using Cache = FluxVariablesCaching::EmptyHeatConductionCache;
    using CacheFiller = FluxVariablesCaching::EmptyCacheFiller<TypeTag>;

    static Scalar flux(const Problem& problem,
                       const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolumeFace& scvf,
                       const ElementFluxVarsCache& elemFluxVarsCache)
    {
        // heat conductivities are always solution dependent (?)
        Scalar tij = calculateTransmissibility_(problem, element, fvGeometry, elemVolVars, scvf);

        // get the inside/outside temperatures
        const auto tInside = elemVolVars[scvf.insideScvIdx()].temperature();
        const auto tOutside = scvf.numOutsideScvs() == 1 ? elemVolVars[scvf.outsideScvIdx()].temperature()
                              : branchingFacetTemperature_(problem, element, fvGeometry, elemVolVars, scvf, tInside, tij);

        return tij*(tInside - tOutside);
    }

private:

    //! compute the temperature at branching facets for network grids
    static Scalar branchingFacetTemperature_(const Problem& problem,
                                             const Element& element,
                                             const FVElementGeometry& fvGeometry,
                                             const ElementVolumeVariables& elemVolVars,
                                             const SubControlVolumeFace& scvf,
                                             Scalar insideTemperature,
                                             Scalar insideTi)
    {
        Scalar sumTi(insideTi);
        Scalar sumTempTi(insideTi*insideTemperature);

        for (unsigned int i = 0; i < scvf.numOutsideScvs(); ++i)
        {
            const auto outsideScvIdx = scvf.outsideScvIdx(i);
            const auto& outsideVolVars = elemVolVars[outsideScvIdx];
            const auto outsideElement = fvGeometry.fvGridGeometry().element(outsideScvIdx);
            const auto& flippedScvf = fvGeometry.flipScvf(scvf.index(), i);

            auto outsideTi = calculateTransmissibility_(problem, outsideElement, fvGeometry, elemVolVars, flippedScvf);
            sumTi += outsideTi;
            sumTempTi += outsideTi*outsideVolVars.temperature();
        }
        return sumTempTi/sumTi;
    }

    static Scalar calculateTransmissibility_(const Problem& problem,
                                             const Element& element,
                                             const FVElementGeometry& fvGeometry,
                                             const ElementVolumeVariables& elemVolVars,
                                             const SubControlVolumeFace& scvf)
    {
        Scalar tij;

        const auto insideScvIdx = scvf.insideScvIdx();
        const auto& insideScv = fvGeometry.scv(insideScvIdx);
        const auto& insideVolVars = elemVolVars[insideScvIdx];

        auto insideLambda = ThermalConductivityModel::effectiveThermalConductivity(insideVolVars, problem.spatialParams(), element, fvGeometry, insideScv);
        Scalar ti = calculateOmega_(scvf, insideLambda, insideScv, insideVolVars.extrusionFactor());

        // for the boundary (dirichlet) or at branching points we only need ti
        if (scvf.boundary() || scvf.numOutsideScvs() > 1)
        {
            tij = scvf.area()*ti;
        }
        // otherwise we compute a tpfa harmonic mean
        else
        {
            const auto outsideScvIdx = scvf.outsideScvIdx();
            const auto& outsideScv = fvGeometry.scv(outsideScvIdx);
            const auto& outsideVolVars = elemVolVars[outsideScvIdx];
            const auto outsideElement = fvGeometry.fvGridGeometry().element(outsideScvIdx);

            auto outsideLambda = ThermalConductivityModel::effectiveThermalConductivity(outsideVolVars,
                                                                                        problem.spatialParams(),
                                                                                        outsideElement,
                                                                                        fvGeometry,
                                                                                        outsideScv);
            Scalar tj;
            if (dim == dimWorld)
                // assume the normal vector from outside is anti parallel so we save flipping a vector
                tj = -1.0*calculateOmega_(scvf, outsideLambda, outsideScv, outsideVolVars.extrusionFactor());
            else
                tj = calculateOmega_(fvGeometry.flipScvf(scvf.index()), outsideLambda, outsideScv, outsideVolVars.extrusionFactor());

            // check for division by zero!
            if (ti*tj <= 0.0)
                tij = 0;
            else
                tij = scvf.area()*(ti * tj)/(ti + tj);
        }

        return tij;
    }

    static Scalar calculateOmega_(const SubControlVolumeFace& scvf,
                                  const DimWorldMatrix& lambda,
                                  const SubControlVolume& scv,
                                  Scalar extrusionFactor)
    {
        GlobalPosition lambdaNormal;
        lambda.mv(scvf.unitOuterNormal(), lambdaNormal);

        auto distanceVector = scvf.ipGlobal();
        distanceVector -= scv.center();
        distanceVector /= distanceVector.two_norm2();

        Scalar omega = lambdaNormal * distanceVector;
        return omega*extrusionFactor;
    }

    static Scalar calculateOmega_(const SubControlVolumeFace& scvf,
                                  Scalar lambda,
                                  const SubControlVolume &scv,
                                  Scalar extrusionFactor)
    {
        auto distanceVector = scvf.ipGlobal();
        distanceVector -= scv.center();
        distanceVector /= distanceVector.two_norm2();

        Scalar omega = lambda * (distanceVector * scvf.unitOuterNormal());
        return omega*extrusionFactor;
    }
};

} // end namespace Dumux

#endif
