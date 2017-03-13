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
 * \brief Base class for the flux variables
 */
#ifndef DUMUX_1P_MIMETIC_IMPLICIT_FLUXVARIABLES_HH
#define DUMUX_1P_MIMETIC_IMPLICIT_FLUXVARIABLES_HH

#include <dumux/implicit/properties.hh>
#include <dumux/discretization/fluxvariablesbase.hh>
#include <dumux/common/math.hh>

namespace Dumux
{

namespace Properties
{
// forward declaration
NEW_PROP_TAG(EnableComponentTransport);
NEW_PROP_TAG(EnableEnergyBalance);
NEW_PROP_TAG(EnableInertiaTerms);
}

// forward declaration
template<class TypeTag, bool enableComponentTransport, bool enableEnergyBalance>
class OnePMimeticFluxVariablesImpl;

/*!
 * \ingroup ImplicitModel
 * \brief The flux variables class
 *        specializations are provided for combinations of physical processes
 * \note  Not all specializations are currently implemented
 */
template<class TypeTag>
using OnePMimeticFluxVariables = OnePMimeticFluxVariablesImpl<TypeTag, GET_PROP_VALUE(TypeTag, EnableComponentTransport),
                                                                 GET_PROP_VALUE(TypeTag, EnableEnergyBalance)>;

/*!
 * \ingroup Discretization
 * \brief Base class for the flux variables
 *        Actual flux variables inherit from this class
 */
// specialization for immiscible, isothermal flow
template<class TypeTag>
class OnePMimeticFluxVariablesImpl<TypeTag, false, false>
: public FluxVariablesBase<TypeTag>
{
    using ParentType = FluxVariablesBase<TypeTag>;
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Element = typename GridView::template Codim<0>::Entity;
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using GlobalFaceVars = typename GET_PROP_TYPE(TypeTag, GlobalFaceVars);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using FluxVariablesCache = typename GET_PROP_TYPE(TypeTag, FluxVariablesCache);
    using CellCenterPrimaryVariables = typename GET_PROP_TYPE(TypeTag, CellCenterPrimaryVariables);
    using FacePrimaryVariables = typename GET_PROP_TYPE(TypeTag, FacePrimaryVariables);
    using IndexType = typename GridView::IndexSet::IndexType;
    using Stencil = std::vector<IndexType>;
    using ElementFluxVariablesCache = typename GET_PROP_TYPE(TypeTag, ElementFluxVariablesCache);

    using DynamicMatrix = Dune::DynamicMatrix<Scalar>;

    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;

    using DimWorldMatrix = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

    using DofTypeIndices = typename GET_PROP(TypeTag, DofTypeIndices);
    typename DofTypeIndices::CellCenterIdx cellCenterIdx;
    typename DofTypeIndices::FaceIdx faceIdx;

public:

    CellCenterPrimaryVariables computeFluxForCellCenter(const Problem& problem,
                                  const Element &element,
                                  const FVElementGeometry& fvGeometry,
                                  const ElementVolumeVariables& elemVolVars,
                                  const GlobalFaceVars& globalFaceVars,
                                  const SubControlVolumeFace &scvf,
                                  const FluxVariablesCache& fluxVarsCache,
                                  const int phaseIdx)
    {
        const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
        const auto& insideVolVars = elemVolVars[insideScv];

        // check if we evaluate the permeability in the volume (for discontinuous fields, default)
        // or at the scvf center for analytical permeability fields (e.g. convergence studies)
//        auto getPermeability = [&problem](const VolumeVariables& volVars,
//                                          const GlobalPosition& scvfIpGlobal)
//                               {
//                                    if (GET_PROP_VALUE(TypeTag, EvaluatePermeabilityAtScvfIP))
//                                        return problem.spatialParams().permeabilityAtPos(scvfIpGlobal);
//                                    else
//                                        return volVars.permeability();
//                               };
//
//        auto perm = getPermeability(insideVolVars, scvf.ipGlobal());
//        const Scalar ti = calculateOmega_(scvf, getPermeability(insideVolVars, scvf.ipGlobal()),
//                                          insideScv, insideVolVars.extrusionFactor());

        // if we are on an inflow/outflow boundary, use the volVars of the element itself
        const auto& outsideVolVars = scvf.boundary() ?  insideVolVars : elemVolVars[scvf.outsideScvIdx()];

        Scalar density = (insideVolVars.density(phaseIdx) + outsideVolVars.density(phaseIdx))*0.5;

        const auto xInside = insideScv.center();
        const auto gInside = problem.gravityAtPos(xInside);
        Scalar potI = insideVolVars.pressure(phaseIdx) - density*(gInside*xInside);
//
//        const auto xFace = scvf.ipGlobal();
//        const auto gFace = problem.gravityAtPos(xFace);
//        Scalar facePot = globalFaceVars.faceVars(scvf.dofIndex()).facePriVars()[phaseIdx] - density*(gFace*xFace);
        CellCenterPrimaryVariables flux(0.0);
        //flux = ti*(potI - facePot);

        DynamicMatrix W = calcMatrix(problem, element, fvGeometry, elemVolVars);
        int indexFace = scvf.localFaceIdx();
        int indexLocal = 0;
        for (auto&& scvfIt : scvfs(fvGeometry))
        {
            const auto xFace = scvfIt.ipGlobal();
            const auto gFace = problem.gravityAtPos(xFace);
            Scalar facePot = globalFaceVars.faceVars(scvfIt.dofIndex()).facePriVars()[phaseIdx] - density*(gFace*xFace);
            flux += scvfIt.area() * W[indexFace][indexLocal] * (potI - facePot);
            indexLocal++;
        }

        return flux * scvf.area();
    }

    void computeCellCenterToCellCenterStencil(Stencil& stencil,
                                              const Problem& problem,
                                              const Element& element,
                                              const FVElementGeometry& fvGeometry,
                                              const SubControlVolumeFace& scvf)
    {
        // the first entry is always the cc dofIdx itself
        if(stencil.empty())
            stencil.push_back(scvf.insideScvIdx());
    }

    void computeCellCenterToFaceStencil(Stencil& stencil,
                                        const Problem& problem,
                                        const Element& element,
                                        const FVElementGeometry& fvGeometry,
                                        const SubControlVolumeFace& scvf)
    {
        stencil.push_back(scvf.dofIndex());
    }

    void computeFaceToCellCenterStencil(Stencil& stencil,
                                        const Problem& problem,
                                        const FVElementGeometry& fvGeometry,
                                        const SubControlVolumeFace& scvf)
    {
        const int eIdx = scvf.insideScvIdx();
        stencil.push_back(scvf.insideScvIdx());

        if(!scvf.boundary())
            stencil.push_back(scvf.outsideScvIdx());
    }

    void computeFaceToFaceStencil(Stencil& stencil,
                                  const Problem& problem,
                                  const FVElementGeometry& fvGeometry,
                                  const SubControlVolumeFace& scvf)
    {
        // the first entries are always the face dofIdx itself and the one of the opposing face
        if(stencil.empty())
        {
            stencil.push_back(scvf.dofIndex());
        }

        for (auto&& scvfIt : scvfs(fvGeometry))
        {
            if(scvfIt.dofIndex() != scvf.dofIndex())
                stencil.push_back(scvfIt.dofIndex());
        }
    }


    DynamicMatrix calcMatrix(const Problem& problem,
                             const Element& element,
                             const FVElementGeometry& fvGeometry,
                             const ElementVolumeVariables& elemVolVars)
    {
        auto getPermeability = [&problem](const VolumeVariables& volVars,
                                          const GlobalPosition& scvfIpGlobal)
                               {
                                    if (GET_PROP_VALUE(TypeTag, EvaluatePermeabilityAtScvfIP))
                                        return problem.spatialParams().permeabilityAtPos(scvfIpGlobal);
                                    else
                                        return volVars.permeability();
                               };
        const auto numFaces = fvGeometry.numScvf();
        DynamicMatrix W(numFaces, numFaces, 0.0);
        DynamicMatrix R(numFaces, dim, 0.0);
        DynamicMatrix N(numFaces, dim, 0.0);
        DynamicMatrix Id(numFaces, numFaces, 0.0);

        for(int i=0; i<numFaces; i++)
            Id[i][i] = 1.0;

        for (auto&& scvf : scvfs(fvGeometry))
        {
            const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());
            const auto& insideVolVars = elemVolVars[insideScv];

            int index = scvf.localFaceIdx();
            N[index] = calculateCoNormal_(scvf, getPermeability(insideVolVars, scvf.ipGlobal()));
            R[index] =  scvf.ipGlobal();
            R[index] -= insideScv.center();
            R[index] *= scvf.area();
        }

        DynamicMatrix LocalK = multiplyMatrices(getTransposed(N),R);
        LocalK.invert();
        DynamicMatrix DTrans = Id;
        DTrans -= multiplyMatrices(multiplyMatrices(R,LocalK),getTransposed(N));
        W = multiplyMatrices(multiplyMatrices(N,LocalK),getTransposed(N));
        DynamicMatrix StabMatrix = multiplyMatrices(getTransposed(DTrans),DTrans);
        StabMatrix *= 0.5*trace(W);
        W += StabMatrix;

        return W;
    }

private:
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

    static GlobalPosition calculateCoNormal_(const SubControlVolumeFace& scvf,
                                  const DimWorldMatrix &K)
    {
        GlobalPosition Knormal;
        K.mv(scvf.unitOuterNormal(), Knormal);

        return Knormal;
    }

    //! compute the transmissibility ti, overload for scalar permeabilites
    static GlobalPosition calculateCoNormal_(const SubControlVolumeFace& scvf,
                                  const Scalar K)
    {
        GlobalPosition Knormal = scvf.unitOuterNormal();
        Knormal *= K;

        return Knormal;
    }

};

} // end namespace

#endif
