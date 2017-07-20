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
#ifndef DUMUX_DISCRETIZATION_MIMETIC_FOURIERS_LAW_HH
#define DUMUX_DISCRETIZATION_MIMETIC_FOURIERS_LAW_HH

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
 * \brief Specialization of Fourier's Law for the Mimetic method.
 */
template <class TypeTag>
class FouriersLawImplementation<TypeTag, DiscretizationMethods::Mimetic>
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
    using ThermalConductivityModel = typename GET_PROP_TYPE(TypeTag, ThermalConductivityModel);
    using GlobalFaceVars = typename GET_PROP_TYPE(TypeTag, GlobalFaceVars);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);

    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;

    using DimWorldMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

    using DynamicMatrix = Dune::DynamicMatrix<Scalar>;

    enum {
        // indices of the primary variables
        temperatureIdx = Indices::temperatureIdx,
    };

    //! The cache used in conjunction with the mpfa Fourier's Law
    class MimeticFouriersLawCache
    {
    public:
        // update cached objects for heat conduction
        void updateHeatConduction(const Problem& problem,
                                    const Element& element,
                                    const FVElementGeometry& fvGeometry,
                                    const ElementVolumeVariables& elemVolVars,
                                    const SubControlVolumeFace &scvf)
        {
            W_ = calculateMatrix(problem, element, fvGeometry, elemVolVars, scvf);
        }

        const DynamicMatrix& WTemp() const
        { return W_; }

        const Scalar& tij() const
        { return ti_; }

    private:
        DynamicMatrix W_;
        Scalar ti_;
    };

    //! Class that fills the cache corresponding to mpfa Darcy's Law
    class MimeticFouriersLawCacheFiller
    {
    public:
        //! Function to fill an MpfaDarcysLawCache of a given scvf
        //! This interface has to be met by any cache filler class for heat conduction quantities
        template<class FluxVariablesCacheFiller>
        static void fill(FluxVariablesCache& scvfFluxVarsCache,
                         const Problem& problem,
                         const Element& element,
                         const FVElementGeometry& fvGeometry,
                         const ElementVolumeVariables& elemVolVars,
                         const SubControlVolumeFace& scvf,
                         const FluxVariablesCacheFiller& fluxVarsCacheFiller)
        {
            scvfFluxVarsCache.updateHeatConduction(problem, element, fvGeometry, elemVolVars, scvf);
        }
    };

public:
    // state the discretization method this implementation belongs to
    static const DiscretizationMethods myDiscretizationMethod = DiscretizationMethods::Mimetic;

    // state the type for the corresponding cache and its filler
    using Cache = MimeticFouriersLawCache;
    using CacheFiller = MimeticFouriersLawCacheFiller;

    static Scalar flux(const Problem& problem,
                        const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const GlobalFaceVars& globalFaceVars,
                        const SubControlVolumeFace& scvf,
                        const ElementFluxVarsCache& elemFluxVarsCache)
    {
        // heat conductivities are always solution dependent (?)
        //Scalar tij = calculateTransmissibility_(problem, element, fvGeometry, elemVolVars, scvf);

        const auto& fluxVarsCache = elemFluxVarsCache[scvf];

        const auto insideScvIdx = scvf.insideScvIdx();
        const auto& insideScv = fvGeometry.scv(insideScvIdx);
        const auto& insideVolVars = elemVolVars[insideScvIdx];

        auto insideLambda = ThermalConductivityModel::effectiveThermalConductivity(insideVolVars, problem.spatialParams(), element, fvGeometry, insideScv);
        const DynamicMatrix& W = fluxVarsCache.WTemp();

        // get the inside/outside temperatures
        const auto tInside = elemVolVars[scvf.insideScvIdx()].temperature();

        Scalar flux = 0.0;

        int indexFace = scvf.localFaceIdx();
        int indexLocal = 0;
        for (auto&& scvfIt : scvfs(fvGeometry))
        {
            Scalar tFace = globalFaceVars.faceVars(scvfIt.dofIndex()).facePriVars()[temperatureIdx];
            flux += W[indexFace][indexLocal] * (tInside - tFace);
            indexLocal++;
        }

        return insideLambda*flux;
    }

private:

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

        tij = scvf.area()*ti;

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

    // The flux variables cache has to be bound to an element prior to flux calculations
    // During the binding, the transmissibilities will be computed and stored using the method below.
    static DynamicMatrix calculateMatrix(const Problem& problem,
                                          const Element& element,
                                          const FVElementGeometry& fvGeometry,
                                          const ElementVolumeVariables& elemVolVars,
                                          const SubControlVolumeFace& scvf)
    {
        const auto numFaces = fvGeometry.numScvf();
        DynamicMatrix W(numFaces, numFaces, 0.0);
        DynamicMatrix R(numFaces, dim, 0.0);
        DynamicMatrix N(numFaces, dim, 0.0);
        DynamicMatrix Id(numFaces, numFaces, 0.0);
        DynamicMatrix C(numFaces, numFaces, 0.0);
        std::vector<Scalar> coNormalNorms(numFaces);

        for(int i=0; i<numFaces; i++)
            Id[i][i] = 1.0;

        for (auto&& scvf : scvfs(fvGeometry))
        {
            const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
            const auto& insideVolVars = elemVolVars[insideScv];

            int index = scvf.localFaceIdx();
            auto coNormal = scvf.unitOuterNormal();
            coNormalNorms[index] = coNormal.two_norm();
            N[index] = coNormal;
            //N[index] /= coNormalNorms[index];
            //N[index] *= scvf.area();
            R[index] =  scvf.ipGlobal();
            R[index] -= insideScv.center();
            C[index][index] = scvf.area();
            R[index] *= scvf.area();
        }

        DynamicMatrix LocalR = multiplyMatrices(getTransposed(R),R);
        LocalR.invert();
        //DynamicMatrix DTrans = multiplyMatrices(C,C);
        //DynamicMatrix DTrans_2 = multiplyMatrices(multiplyMatrices(R,LocalR),getTransposed(R));

        DynamicMatrix DTrans = Id;
        DTrans -= multiplyMatrices(multiplyMatrices(R,LocalR),getTransposed(R));

        DynamicMatrix LocalK = multiplyMatrices(getTransposed(N),R);
        LocalK.invert();
        W = multiplyMatrices(multiplyMatrices(N,LocalK),getTransposed(N));
        DynamicMatrix StabMatrix = DTrans; //multiplyMatrices(getTransposed(DTrans),DTrans);
        StabMatrix *= 0.5*trace(W);
        W += StabMatrix;

        DynamicMatrix sol = multiplyMatrices(C,multiplyMatrices(W,C));

        return sol;
    }
};

} // end namespace Dumux

#endif
