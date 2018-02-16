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
 *        volume and mass fluxes of fluid phases over a face of a finite volume by means
 *        of the Darcy approximation. Specializations are provided for the different discretization methods.
 */
#ifndef DUMUX_DISCRETIZATION_MIMETIC_DARCYS_LAW_HH
#define DUMUX_DISCRETIZATION_MIMETIC_DARCYS_LAW_HH

#include <memory>

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
NEW_PROP_TAG(ProblemEnableGravity);
NEW_PROP_TAG(GlobalFaceVars);
}

/*!
 * \ingroup DarcysLaw
 * \brief Specialization of Darcy's Law for the Mimetic method.
 */
template <class TypeTag>
class DarcysLawImplementation<TypeTag, DiscretizationMethods::Mimetic>
{
    using Implementation = DarcysLawImplementation<TypeTag, DiscretizationMethods::Mimetic>;
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using ElementFluxVarsCache = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);
    using FluxVariablesCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);

    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalFaceVars = typename GET_PROP_TYPE(TypeTag, GlobalFaceVars);
    using IndexType = typename GridView::IndexSet::IndexType;
    using Stencil = std::vector<IndexType>;

    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;

    using DimWorldMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

    using DynamicMatrix = Dune::DynamicMatrix<Scalar>;

    class MimeticDarcysLawCache
    {
    public:
        void updateAdvection(const Problem& problem,
                             const Element& element,
                             const FVElementGeometry& fvGeometry,
                             const ElementVolumeVariables& elemVolVars,
                             const SubControlVolumeFace &scvf)
        {
            W_ = Implementation::calculateMatrix(problem, element, fvGeometry, elemVolVars, scvf);
            //ti_ = Implementation::calculateTransmissibilities(problem, element, fvGeometry, elemVolVars, scvf);

        }

        const DynamicMatrix& W() const
        { return W_; }

        const Scalar& tij() const
        { return ti_; }

    private:
        DynamicMatrix W_;
        Scalar ti_;
    };

    //! Class that fills the cache corresponding to Mimetic Darcy's Law
    class MimeticDarcysLawCacheFiller
    {
    public:
        //! Function to fill a MimeticDarcysLawCache of a given scvf
        //! This interface has to be met by any advection-related cache filler class
        template<class FluxVariablesCacheFiller>
        static void fill(FluxVariablesCache& scvfFluxVarsCache,
                         const Problem& problem,
                         const Element& element,
                         const FVElementGeometry& fvGeometry,
                         const ElementVolumeVariables& elemVolVars,
                         const SubControlVolumeFace& scvf,
                         const FluxVariablesCacheFiller& fluxVarsCacheFiller)
        {
            scvfFluxVarsCache.updateAdvection(problem, element, fvGeometry, elemVolVars, scvf);
        }
    };

public:
    // state the discretization method this implementation belongs to
    static const DiscretizationMethods myDiscretizationMethod = DiscretizationMethods::Mimetic;

    // state the type for the corresponding cache and its filler
    using Cache = MimeticDarcysLawCache;
    using CacheFiller = MimeticDarcysLawCacheFiller;

    static Scalar flux(const Problem& problem,
                       const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const GlobalFaceVars& globalFaceVars,
                       const SubControlVolumeFace& scvf,
                       int phaseIdx,
                       const ElementFluxVarsCache& elemFluxVarsCache)
    {
        Scalar flux = 0.0;
        const auto& fluxVarsCache = elemFluxVarsCache[scvf];
        const auto facePriVars = globalFaceVars.faceVars(scvf.dofIndex()).facePriVars();
        //This index is used for the case that no capillary pressure is included
        int numVarsFace = facePriVars.size();
        int phaseIdxFace = std::min(phaseIdx,numVarsFace-1);
        flux = globalFaceVars.faceVars(scvf.dofIndex()).facePriVars()[phaseIdxFace];
        flux *= -scvf.area()*scvf.fluxMultiplier();

        return  flux;
    }

    // The flux variables cache has to be bound to an element prior to flux calculations
    // During the binding, the transmissibilities will be computed and stored using the method below.
    static DynamicMatrix calculateMatrix(const Problem& problem,
                                                      const Element& element,
                                                      const FVElementGeometry& fvGeometry,
                                                      const ElementVolumeVariables& elemVolVars,
                                                      const SubControlVolumeFace& scvf)
    {
        auto getPermeability = [&problem](const VolumeVariables& volVars,
                                          const GlobalPosition& scvfIpGlobal)
                               {
//                                    if (GET_PROP_VALUE(TypeTag, EvaluatePermeabilityAtScvfIP))
//                                        return problem.spatialParams().permeabilityAtPos(scvfIpGlobal);
//                                    else
                                        return volVars.permeability();
                               };

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
            auto coNormal = calculateCoNormal_(scvf, getPermeability(insideVolVars, scvf.ipGlobal()));
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

        W.invert();
        DynamicMatrix sol = multiplyMatrices(C,multiplyMatrices(W,C));

        for (auto&& scvf : scvfs(fvGeometry))
        {
            const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
            const auto& insideVolVars = elemVolVars[insideScv];

            int index = scvf.localFaceIdx();

//            const auto bcTypes = problem.boundaryTypes(element, scvf);
//            if (bcTypes.hasNeumann())
//            {
//                for (auto&& scvfIt : scvfs(fvGeometry))
//                {
//                    int indexIt = scvfIt.localFaceIdx();
//                    W[indexIt][index] = 0.0;
//                }
//            }
            //W[index] *= coNormalNorms[index];
        }


//        const auto numFaces = fvGeometry.numScvf();
//        DynamicMatrix W(numFaces, numFaces, 0.0);
//        for (auto&& scvf : scvfs(fvGeometry))
//        {
//            int index = scvf.localFaceIdx();
//            W[index][index] = 1.0/calculateTransmissibilities(problem, element, fvGeometry, elemVolVars, scvf);
//        }


        return W;
    }

    // The flux variables cache has to be bound to an element prior to flux calculations
    // During the binding, the transmissibilities will be computed and stored using the method below.
    static Scalar calculateTransmissibilities(const Problem& problem,
                                              const Element& element,
                                              const FVElementGeometry& fvGeometry,
                                              const ElementVolumeVariables& elemVolVars,
                                              const SubControlVolumeFace& scvf)
    {
        Scalar tij;

        const auto insideScvIdx = scvf.insideScvIdx();
        const auto& insideScv = fvGeometry.scv(insideScvIdx);
        const auto& insideVolVars = elemVolVars[insideScvIdx];
        // check if we evaluate the permeability in the volume (for discontinuous fields, default)
        // or at the scvf center for analytical permeability fields (e.g. convergence studies)
        auto getPermeability = [&problem](const VolumeVariables& volVars,
                                          const GlobalPosition& scvfIpGlobal)
                               {
//                                    if (GET_PROP_VALUE(TypeTag, EvaluatePermeabilityAtScvfIP))
//                                        return problem.spatialParams().permeabilityAtPos(scvfIpGlobal);
//                                    else
                                        return volVars.permeability();
                               };

        const Scalar ti = calculateOmega_(scvf, getPermeability(insideVolVars, scvf.ipGlobal()),
                                          insideScv, insideVolVars.extrusionFactor());

        // for the boundary (dirichlet) or at branching points we only need ti
//        if (scvf.boundary() || scvf.numOutsideScvs() > 1)
//        {
//            tij = scvf.area()*ti;
//        }
//
//        // otherwise we compute a tpfa harmonic mean
//        else
//        {
//            const auto outsideScvIdx = scvf.outsideScvIdx();
//            // as we assemble fluxes from the neighbor to our element the outside index
//            // refers to the scv of our element, so we use the scv method
//            const auto& outsideScv = fvGeometry.scv(outsideScvIdx);
//            const auto& outsideVolVars = elemVolVars[outsideScvIdx];
//
//            const Scalar tj = [&]()
//            {
//                    return -1.0*calculateOmega_(scvf, getPermeability(outsideVolVars, scvf.ipGlobal()),
//                                            outsideScv, outsideVolVars.extrusionFactor());
//
//            }();
//
//            // harmonic mean (check for division by zero!)
//            if (ti*tj <= 0.0)
//                tij = 0;
//            else
//                tij = scvf.area()*(ti * tj)/(ti + tj);
//        }

//        return scvf.area()*ti;
        return ti/scvf.area();
    }

private:
    static GlobalPosition calculateCoNormal_(const SubControlVolumeFace& scvf,
                                  const DimWorldMatrix &K)
    {
        GlobalPosition Knormal;
        K.mv(scvf.unitOuterNormal(), Knormal);

        return Knormal;
    }

    static GlobalPosition calculateCoNormal_(const SubControlVolumeFace& scvf,
                                  const Scalar K)
    {
        GlobalPosition Knormal = scvf.unitOuterNormal();
        Knormal *= K;

        return Knormal;
    }

    //! compute the transmissibility ti, overload for tensor permeabilites
    static Scalar calculateOmega_(const SubControlVolumeFace& scvf,
                                  const DimWorldMatrix &K,
                                  const SubControlVolume &scv,
                                  Scalar extrusionFactor)
    {
        GlobalPosition Knormal;
        K.mv(scvf.unitOuterNormal(), Knormal);

        auto distanceVector = scvf.ipGlobal();
        distanceVector -= scv.center();
        distanceVector /= distanceVector.two_norm2();

        Scalar omega = Knormal * distanceVector;
        omega *= extrusionFactor;

        return omega;
    }

    //! compute the transmissibility ti, overload for scalar permeabilites
    static Scalar calculateOmega_(const SubControlVolumeFace& scvf,
                                  const Scalar K,
                                  const SubControlVolume &scv,
                                  Scalar extrusionFactor)
    {
        auto distanceVector = scvf.ipGlobal();
        distanceVector -= scv.center();
        distanceVector /= distanceVector.two_norm2();

        Scalar omega = K * (distanceVector * scvf.unitOuterNormal());
        omega *= extrusionFactor;

        return omega;
    }
};

} // end namespace Dumux

#endif
