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

        // Get the inside and outside volume variables
        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& insideVolVars = elemVolVars[insideScv];
        const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];

        if (GET_PARAM_FROM_GROUP(TypeTag, bool, Problem, EnableGravity))
        {
            // do averaging for the density over all neighboring elements
            const auto rho = [&]()
            {
                // boundaries
                if (scvf.boundary())
                    return insideVolVars.density(phaseIdx);

                // inner faces with two neighboring elements
                else if (scvf.numOutsideScvs() == 1)
                    return (insideVolVars.density(phaseIdx) + outsideVolVars.density(phaseIdx))*0.5;

                // inner faces in networks (general case)
                else
                {
                    Scalar rho(insideVolVars.density(phaseIdx));
                    for (unsigned int i = 0; i < scvf.numOutsideScvs(); ++i)
                    {
                        const auto outsideScvIdx = scvf.outsideScvIdx(i);
                        const auto& outsideVolVars = elemVolVars[outsideScvIdx];
                        rho += outsideVolVars.density(phaseIdx);
                    }
                    return rho/(scvf.numOutsideScvs()+1);
                }
            }();

            // ask for the gravitational acceleration in the inside neighbor
            const auto xInside = insideScv.center();
            const auto gInside = problem.gravityAtPos(xInside);
            const auto hInside = insideVolVars.pressure(phaseIdx) - rho*(gInside*xInside);

            const DynamicMatrix& W = fluxVarsCache.W();

            int indexFace = scvf.localFaceIdx();
            int indexLocal = 0;
            for (auto&& scvfIt : scvfs(fvGeometry))
            {
                const auto xFace = scvfIt.ipGlobal();
                const auto gFace = problem.gravityAtPos(xFace);
                Scalar hFace = globalFaceVars.faceVars(scvfIt.dofIndex()).facePriVars()[0] - rho*(gFace*xFace);
                flux += W[indexFace][indexLocal] * (hInside - hFace);
                indexLocal++;
            }
//            const auto xFace = scvf.ipGlobal();
//            const auto gFace = problem.gravityAtPos(xFace);
//            Scalar hFace = globalFaceVars.faceVars(scvf.dofIndex()).facePriVars()[0] - rho*(gFace*xFace);
//            flux = fluxVarsCache.tij()* (hInside - hFace);
        }
        else // no gravity
        {
            const auto pInside = insideVolVars.pressure(phaseIdx);
            const DynamicMatrix& W = fluxVarsCache.W();

            int indexFace = scvf.localFaceIdx();
            int indexLocal = 0;
            for (auto&& scvfIt : scvfs(fvGeometry))
            {
                Scalar pFace = globalFaceVars.faceVars(scvfIt.dofIndex()).facePriVars()[0];
                flux += W[indexFace][indexLocal] * (pInside - pFace);
                indexLocal++;
            }

//            Scalar pFace = globalFaceVars.faceVars(scvf.dofIndex()).facePriVars()[0];
//            flux = fluxVarsCache.tij()* (pInside - pFace);

//            if(!scvf.boundary())
//                flux = fluxVarsCache.tij()* (pInside - outsideVolVars.pressure(phaseIdx));

//            Scalar pOutside = outsideVolVars.pressure(phaseIdx);
//            if(scvf.boundary())
//                flux = fluxVarsCache.tij()* (pInside - pOutside);
        }

        return flux;
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
        std::vector<std::vector<GlobalPosition>> yvals;
        yvals.resize(numFaces);

        for (auto&& scvf : scvfs(fvGeometry))
        {
            const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
            auto cellVolume = insideScv.volume();

            int indexI = scvf.localFaceIdx();
            for (auto&& scvfIt : scvfs(fvGeometry))
            {
                int indexJ = scvfIt.localFaceIdx();
                yvals[indexI].resize(numFaces);

                auto normalI = scvf.unitOuterNormal();
                auto normalJ = scvfIt.unitOuterNormal();

                auto areaI = scvf.area();
                auto areaJ = scvfIt.area();

                auto distVecI = scvf.ipGlobal();
                distVecI -= insideScv.center();

                auto distI = distVecI*normalI;

                if(indexI == indexJ)
                {
                    yvals[indexI][indexJ] = normalI;
                    yvals[indexI][indexJ] *= (areaI/cellVolume + std::sqrt(dimWorld)/distI*(1-areaI/cellVolume*distI));
                }
                else
                {
                    auto partI = normalJ;
                    partI *= areaJ/cellVolume;
                    yvals[indexI][indexJ] = normalI;
                    yvals[indexI][indexJ] *= -(std::sqrt(dimWorld) * areaJ/(distI*cellVolume)*(normalJ*distVecI));
                    yvals[indexI][indexJ] += partI;
                }
            }
        }

        DynamicMatrix W(numFaces, numFaces, 0.0);

        for (auto&& scvf : scvfs(fvGeometry))
        {
            const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
            const auto& insideVolVars = elemVolVars[insideScv];
            auto cellVolume = insideScv.volume();

            int indexI = scvf.localFaceIdx();
            for (auto&& scvfIt : scvfs(fvGeometry))
            {
                int indexJ = scvfIt.localFaceIdx();
                Scalar coeff = 0.0;

                for (auto&& scvfItIt : scvfs(fvGeometry))
                {
                    int indexK = scvfItIt.localFaceIdx();
                    auto K = getPermeability(insideVolVars, scvfItIt.ipGlobal());
                    GlobalPosition Ky;
                    K.mv(yvals[indexK][indexJ], Ky);
                    auto normalK = scvfItIt.unitOuterNormal();
                    auto areaK = scvfItIt.area();
                    auto distVecK = scvfItIt.ipGlobal();
                    distVecK -= insideScv.center();
                    auto distK = normalK*distVecK;
                    coeff += (yvals[indexK][indexI] * Ky)*(areaK*distK/dimWorld);
                }

                W[indexI][indexJ] = coeff;
            }
        }

        return W;
    }

    // The flux variables cache has to be bound to an element prior to flux calculations
    // During the binding, the transmissibilities will be computed and stored using the method below.
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
        if (scvf.boundary() || scvf.numOutsideScvs() > 1)
        {
            tij = scvf.area()*ti;
        }

        // otherwise we compute a tpfa harmonic mean
        else
        {
            const auto outsideScvIdx = scvf.outsideScvIdx();
            // as we assemble fluxes from the neighbor to our element the outside index
            // refers to the scv of our element, so we use the scv method
            const auto& outsideScv = fvGeometry.scv(outsideScvIdx);
            const auto& outsideVolVars = elemVolVars[outsideScvIdx];

            const Scalar tj = [&]()
            {
                    return -1.0*calculateOmega_(scvf, getPermeability(outsideVolVars, scvf.ipGlobal()),
                                            outsideScv, outsideVolVars.extrusionFactor());

            }();

            // harmonic mean (check for division by zero!)
            if (ti*tj <= 0.0)
                tij = 0;
            else
                tij = scvf.area()*(ti * tj)/(ti + tj);
        }

        return scvf.area()*ti;
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
