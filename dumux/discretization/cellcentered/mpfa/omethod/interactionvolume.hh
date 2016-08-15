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
 * \brief Base classes for interaction volume of mpfa methods.
 */
#ifndef DUMUX_DISCRETIZATION_CC_MPFA_O_INTERACTIONVOLUME_HH
#define DUMUX_DISCRETIZATION_CC_MPFA_O_INTERACTIONVOLUME_HH

#include <dumux/common/math.hh>
#include <dumux/implicit/cellcentered/mpfa/properties.hh>
#include <dumux/discretization/cellcentered/mpfa/interactionvolumeseed.hh>
#include <dumux/discretization/cellcentered/mpfa/facetypes.hh>
#include <dumux/discretization/cellcentered/mpfa/methods.hh>

#include "localsubcontrolentityseeds.hh"
#include "localsubcontrolentities.hh"

namespace Dumux
{
/*!
 * \ingroup Mpfa
 * \brief Base class for the interaction volumes of the mpfa-o method
 */
template<class TypeTag>
class CCMpfaInteractionVolumeImplementation<TypeTag, MpfaMethods::oMethod>
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using MpfaHelper = typename GET_PROP_TYPE(TypeTag, MpfaHelper);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);

    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalIndexType = typename GridView::IndexSet::IndexType;
    using LocalIdxType = unsigned int;

    using LocalScvSeedType = CCMpfaOLocalScvSeed<GlobalIndexType, LocalIdxType>;
    using LocalScvfSeedType = CCMpfaOLocalScvfSeed<GlobalIndexType, LocalIdxType>;

    using LocalScvType = CCMpfaOLocalScv<TypeTag>;
    using LocalScvfType = CCMpfaOLocalScvf<TypeTag>;

    static const int dim = GridView::dimension;
    using DimVector = Dune::FieldVector<Scalar, dim>;
    using Tensor = Dune::FieldMatrix<Scalar, dim, dim>;

    using Stencil = std::vector<GlobalIndexType>;
    using LocalIndexSet = std::vector<LocalIdxType>;

    using Matrix = Dune::DynamicMatrix<Scalar>;
    using Vector = typename Matrix::row_type;

public:

    // the interaction volume seed type to be exported
    using LocalIndexType = LocalIdxType;
    using Seed = CCMpfaInteractionVolumeSeed<LocalScvSeedType, LocalScvfSeedType>;

    CCMpfaInteractionVolumeImplementation(const Seed& seed,
                                          const Problem& problem,
                                          const FVElementGeometry& fvGeometry,
                                          const ElementVolumeVariables& elemVolVars)
    : problemRef_(problem),
      fvGeometryRef_(fvGeometry),
      elemVolVarsRef_(elemVolVars),
      onBoundary_(seed.onBoundary())
    {
        // create local sub control entities from the seed
        createLocalEntities_(seed);

        // fill the dof stencil
        stencil_.reserve(localScvs_.size());
        volVarsPositions_.reserve(localScvs_.size());
        for (const auto& localScv : localScvs_)
        {
            stencil_.push_back(localScv.globalIndex());
            volVarsPositions_.push_back(localScv.center());
        }

        // fill vol vars stencil
        volVarsStencil_ = stencil_;
        for (const auto& localScvf : localScvfs_)
            if (localScvf.faceType() == MpfaFaceTypes::dirichlet || localScvf.faceType() == MpfaFaceTypes::interiorDirichlet)
            {
                volVarsStencil_.push_back(localScvf.outsideGlobalScvIndex());
                volVarsPositions_.push_back(localScvf.ip());
            }

        // set up local index sets of flux and dirichlet faces
        LocalIndexType localScvfIdx = 0;
        for (const auto& localScvf : localScvfs_)
        {
            auto faceType = localScvf.faceType();
            if (faceType != MpfaFaceTypes::dirichlet && faceType != MpfaFaceTypes::interiorDirichlet)
                fluxFaceIndexSet_.push_back(localScvfIdx++);
            else
                dirichletFaceIndexSet_.push_back(localScvfIdx++);
        }

        neumannFluxes_ = Vector(fluxFaceIndexSet_.size(), 0.0);
    }

    template<typename GetTensorFunction>
    void solveLocalSystem(const GetTensorFunction& getTensor)
    {
        const std::size_t numFluxFaces = fluxScvfIndexSet_().size();

        // if only dirichlet faces are present, assemble T_ directly
        if (numFluxFaces == 0)
            return assemblePureDirichletSystem_(getTensor);

        const std::size_t numFaces = localScvfs_.size();
        const std::size_t numPotentials = volVarsStencil().size();

        // the local matrices
        Matrix A(numFluxFaces, numFluxFaces, 0.0);
        Matrix B(numFluxFaces, numPotentials, 0.0);
        Matrix C(numFaces, numFluxFaces, 0.0);
        Matrix D(numFaces, numPotentials, 0.0);

        assembleLocalMatrices_(getTensor, A, B, C, D);

        // solve local system and store matrices
        Matrix copy(B);
        A.invert();
        AinvB_ = B.leftmultiply(A);
        CAinv_ = C.rightmultiply(A);
        T_ = multiplyMatrices(CAinv_, copy);
        T_ += D;
    }

    template<typename UpwindFactorFunction>
    void assembleNeumannFluxes(const UpwindFactorFunction& upwindFactor,
                               const unsigned int eqIdx)
    {

        LocalIndexType fluxFaceIdx = 0;
        for (const auto localFluxFaceIdx : fluxFaceIndexSet_)
        {
            const auto& localScvf = localScvf_(localFluxFaceIdx);
            const auto faceType = localScvf.faceType();
            if (faceType == MpfaFaceTypes::neumann || faceType == MpfaFaceTypes::interiorNeumann)
            {
                // boundary neumann fluxes stay zero when tpfa boundary handling is on
                if (GET_PROP_VALUE(TypeTag, UseTpfaBoundary) && faceType == MpfaFaceTypes::neumann)
                    continue;

                const auto& element = localElement_(localScvf.insideLocalScvIndex());
                const auto& globalScvf = fvGeometry_().scvf(localScvf.insideGlobalScvfIndex());
                auto neumannFlux = problem_().neumann(element, this->fvGeometry_(), this->elemVolVars_(), globalScvf)[eqIdx];
                neumannFlux *= globalScvf.area();

                const auto& insideScv = fvGeometry_().scv(globalScvf.insideScvIdx());
                const auto& volVars = elemVolVars_()[insideScv];
                neumannFlux /= upwindFactor(volVars);

                neumannFluxes_[fluxFaceIdx] = neumannFlux;
            }

            fluxFaceIdx++;
        }
    }

    const bool onBoundary() const
    { return onBoundary_; }

    const Stencil& stencil() const
    { return stencil_; }

    const Stencil& volVarsStencil() const
    { return volVarsStencil_; }

    const std::vector<DimVector>& volVarsPositions() const
    { return volVarsPositions_; }

    const Stencil globalScvfs() const
    {
        Stencil stencil;
        stencil.reserve(localScvfs_.size());

        for (const auto& localScvf : localScvfs_)
        {
            stencil.push_back(localScvf.insideGlobalScvfIndex());
            if (!localScvf.boundary())
                stencil.push_back(localScvf.outsideGlobalScvfIndex());
        }

        return stencil;
    }

    std::pair<LocalIndexType, bool> getLocalIndexPair(const SubControlVolumeFace& scvf) const
    {
        auto scvfGlobalIdx = scvf.index();

        LocalIndexType localIdx = 0;
        for (const auto& localScvf : localScvfs_)
        {
            if (localScvf.insideGlobalScvfIndex() == scvfGlobalIdx)
                return std::make_pair(localIdx, false);

            if (!localScvf.boundary() && localScvf.outsideGlobalScvfIndex() == scvfGlobalIdx)
                return std::make_pair(localIdx, true);

            localIdx++;
        }

        DUNE_THROW(Dune::InvalidStateException, "Could not find the local scv face in the interaction volume for the given scvf with index: " << scvf.index());
    }

    Vector getTransmissibilities(const std::pair<LocalIndexType, bool>& localIndexPair) const
    {
        auto tij = T_[localIndexPair.first];

        if (localIndexPair.second)
            tij *= -1.0;
        return tij;
    }

    Scalar getNeumannFlux(const std::pair<LocalIndexType, bool>& localIndexPair) const
    {
        if (fluxScvfIndexSet_().size() == 0)
            return 0.0;

        auto flux = CAinv_[localIndexPair.first] * neumannFluxes_;
        if (localIndexPair.second)
            return -1.0*flux;
        return flux;
    }

private:

    const LocalScvfType& localScvf_(const LocalIndexType localScvfIdx) const
    { return localScvfs_[localScvfIdx]; }

    const LocalScvType& localScv_(const LocalIndexType localScvIdx) const
    { return localScvs_[localScvIdx]; }

    const LocalIndexSet& fluxScvfIndexSet_() const
    { return fluxFaceIndexSet_; }

    const LocalIndexSet& dirichletScvfIndexSet_() const
    { return dirichletFaceIndexSet_; }

    const Element& localElement_(const LocalIndexType localScvIdx) const
    { return localElements_[localScvIdx]; }

    void createLocalEntities_(const Seed& seed)
    {
        auto numScvs = seed.scvSeeds().size();
        auto numScvfs = seed.scvfSeeds().size();

        localElements_.reserve(numScvs);
        localScvs_.reserve(numScvs);
        localScvfs_.reserve(numScvfs);

        for (const auto& scvSeed : seed.scvSeeds())
        {
            auto element = problem_().model().globalFvGeometry().element(scvSeed.globalIndex());
            localScvs_.emplace_back(LocalScvType(element, fvGeometry_(), scvSeed));
            localElements_.emplace_back(std::move(element));
        }

        for (const auto& scvfSeed : seed.scvfSeeds())
        {
            // we have to use the "inside" scv face here
            const auto& scvf = fvGeometry_().scvf(scvfSeed.insideGlobalScvfIndex());
            localScvfs_.emplace_back(LocalScvfType(scvfSeed, scvf));
        }
    }

    template<typename GetTensorFunction>
    void assembleLocalMatrices_(const GetTensorFunction& getTensor,
                                Matrix& A, Matrix& B, Matrix& C, Matrix& D)
    {
        const Scalar xi = GET_PROP_VALUE(TypeTag, Xi);
        const std::size_t numLocalScvs = localScvs_.size();

        // loop over the local faces
        LocalIndexType rowIdx = 0;
        for (const auto& localScvf : localScvfs_)
        {
            auto faceType = localScvf.faceType();
            bool hasUnknown = faceType != MpfaFaceTypes::dirichlet && faceType != MpfaFaceTypes::interiorDirichlet;
            LocalIndexType idxInFluxFaces = hasUnknown ? findLocalIndex_(fluxScvfIndexSet_(), rowIdx) : -1;

            // get diffusion tensor in "positive" sub volume
            const auto posLocalScvIdx = localScvf.insideLocalScvIndex();
            const auto& posLocalScv = localScv_(posLocalScvIdx);
            const auto& posGlobalScv = fvGeometry_().scv(posLocalScv.globalIndex());
            auto element = localElement_(posLocalScvIdx);
            auto tensor = getTensor(element, elemVolVars_()[posGlobalScv], posGlobalScv);

            // the omega factors of the "positive" sub volume
            auto posWijk = calculateOmegas_(posLocalScv, localScvf, tensor);
            posWijk *= problem_().boxExtrusionFactor(element, posGlobalScv);

            // Check the local directions of the positive sub volume
            for (int localDir = 0; localDir < dim; localDir++)
            {
                auto curLocalScvfIdx = posLocalScv.localScvfIndex(localDir);
                const auto& curLocalScvf = localScvf_(curLocalScvfIdx);
                auto curFaceType = curLocalScvf.faceType();
                bool curFaceHasUnknown = curFaceType != MpfaFaceTypes::dirichlet && curFaceType != MpfaFaceTypes::interiorDirichlet;

                // First, add the entries associated with face pressures (unkown or dirichlet)
                if (curFaceHasUnknown)
                {
                    // we need the index of the current local scvf in the flux face indices
                    auto curIdxInFluxFaces = findLocalIndex_(fluxScvfIndexSet_(), curLocalScvfIdx);

                    C[rowIdx][curIdxInFluxFaces] += posWijk[localDir];
                    if (hasUnknown)
                    {
                        if (faceType == MpfaFaceTypes::interiorNeumann)
                        {
                            A[idxInFluxFaces][curIdxInFluxFaces] += xi*posWijk[localDir];
                            // If facet coupling active, eventually add entry coming from the lowdim domain
                            if (GET_PROP_VALUE(TypeTag, FacetCoupling) && curIdxInFluxFaces == idxInFluxFaces)
                                DUNE_THROW(Dune::NotImplemented, "Facet coupling");
                        }
                        else
                            A[idxInFluxFaces][curIdxInFluxFaces] += posWijk[localDir];
                    }
                }
                else
                {
                    // the current face is a Dirichlet face and creates entries in D & eventually B
                    auto curIdxInDiriFaces = findLocalIndex_(dirichletScvfIndexSet_(), curLocalScvfIdx);
                    D[rowIdx][numLocalScvs + curIdxInDiriFaces] += posWijk[localDir];

                    if (hasUnknown)
                        B[idxInFluxFaces][numLocalScvs + curIdxInDiriFaces] -= posWijk[localDir];
                }

                // add entries related to pressures at the scv centers (dofs)
                D[rowIdx][posLocalScvIdx] -= posWijk[localDir];

                if (hasUnknown)
                {
                    if (faceType == MpfaFaceTypes::interiorNeumann)
                        B[idxInFluxFaces][posLocalScvIdx] += xi*posWijk[localDir];
                    else
                        B[idxInFluxFaces][posLocalScvIdx] += posWijk[localDir];
                }
            }

            // If not on a boundary or interior dirichlet face, add entries for the "negative" scv
            if (faceType == MpfaFaceTypes::interior || faceType == MpfaFaceTypes::interiorNeumann)
            {
                const auto negLocalScvIdx = localScvf.outsideLocalScvIndex();
                const auto& negLocalScv = localScv_(negLocalScvIdx);
                const auto& negGlobalScv = fvGeometry_().scv(negLocalScv.globalIndex());
                auto negElement = localElement_(negLocalScvIdx);;
                auto negTensor = getTensor(negElement, elemVolVars_()[negGlobalScv], negGlobalScv);

                // the omega factors of the "negative" sub volume
                auto negWijk = calculateOmegas_(negLocalScv, localScvf, negTensor);
                negWijk *= problem_().boxExtrusionFactor(negElement, negGlobalScv);

                // Check local directions of negative sub volume
                for (int localDir = 0; localDir < dim; localDir++)
                {
                    auto curLocalScvfIdx = negLocalScv.localScvfIndex(localDir);
                    const auto& curLocalScvf = localScvf_(curLocalScvfIdx);
                    auto curFaceType = curLocalScvf.faceType();
                    bool curFaceHasUnknown = curFaceType != MpfaFaceTypes::dirichlet && curFaceType != MpfaFaceTypes::interiorDirichlet;

                    if (curFaceHasUnknown)
                    {
                        // we need the index of the current local scvf in the flux face indices
                        auto curIdxInFluxFaces = findLocalIndex_(fluxScvfIndexSet_(), curLocalScvfIdx);

                        if (faceType == MpfaFaceTypes::interiorNeumann)
                            A[idxInFluxFaces][curIdxInFluxFaces] -= (1.0-xi)*negWijk[localDir];
                        else
                            A[idxInFluxFaces][curIdxInFluxFaces] -= negWijk[localDir];
                    }
                    else
                    {
                        // the current face is a Dirichlet face and creates entries in B
                        auto curIdxInDiriFaces = findLocalIndex_(dirichletScvfIndexSet_(), curLocalScvfIdx);
                        B[idxInFluxFaces][numLocalScvs + curIdxInDiriFaces] += negWijk[localDir];
                    }

                    // add entries to matrix B
                    if (faceType == MpfaFaceTypes::interiorNeumann)
                        B[idxInFluxFaces][negLocalScvIdx] -= (1.0-xi)*negWijk[localDir];
                    else
                        B[idxInFluxFaces][negLocalScvIdx] -= negWijk[localDir];
                }
            }
            // go to the next face
            rowIdx++;
        }
    }

    template<typename GetTensorFunction>
    void assemblePureDirichletSystem_(const GetTensorFunction& getTensor)
    {
        const std::size_t numLocalScvs = localScvs_.size();
        const std::size_t numFaces = localScvfs_.size();
        const std::size_t numPotentials = volVarsStencil().size();

        // resize matrices, only T_ will have entries
        T_.resize(numFaces, numPotentials, 0.0);
        AinvB_.resize(0, 0);
        CAinv_.resize(0, 0);

        // Loop over all the faces, in this case these are all dirichlet boundaries
        LocalIndexType rowIdx = 0;
        for (const auto& localScvf : localScvfs_)
        {
            // get diffusion tensor in "positive" sub volume
            const auto posLocalScvIdx = localScvf.insideLocalScvIndex();
            const auto& posLocalScv = localScv_(posLocalScvIdx);
            const auto& posGlobalScv = fvGeometry_().scv(posLocalScv.globalIndex());
            auto element = localElement_(posLocalScvIdx);
            auto tensor = getTensor(element, elemVolVars_()[posGlobalScv], posGlobalScv);

            // the omega factors of the "positive" sub volume
            auto posWijk = calculateOmegas_(posLocalScv, localScvf, tensor);
            posWijk *= problem_().boxExtrusionFactor(element, posGlobalScv);

            for (int localDir = 0; localDir < dim; localDir++)
            {
                auto curLocalScvfIdx = posLocalScv.localScvfIndex(localDir);
                auto curIdxInDiriFaces = findLocalIndex_(dirichletScvfIndexSet_(), curLocalScvfIdx);

                T_[rowIdx][numLocalScvs + curIdxInDiriFaces] += posWijk[localDir];
                T_[rowIdx][posLocalScvIdx] -= posWijk[localDir];
            }

            // go to the next face
            rowIdx++;
        }
    }

    DimVector calculateOmegas_(const LocalScvType& localScv,
                               const LocalScvfType& localScvf,
                               const Tensor& T) const
    {
        DimVector wijk;
        DimVector tmp;
        for (int dir = 0; dir < dim; ++dir)
        {
            T.mv(localScv.innerNormal(dir), tmp);
            wijk[dir] = tmp*localScvf.unitOuterNormal();
        }
        wijk *= localScvf.area();
        wijk /= localScv.detX();
        wijk *= -1.0;

        return wijk;
    }

    DimVector calculateOmegas_(const LocalScvType& localScv,
                               const LocalScvfType& localScvf,
                               const Scalar t) const
    {
        DimVector wijk;
        DimVector tmp(localScvf.unitOuterNormal());
        tmp *= t;

        for (int dir = 0; dir < dim; ++dir)
            wijk[dir] = tmp*localScv.innerNormal(dir);
        wijk *= localScvf.area();
        wijk /= localScv.detX();
        wijk *= -1.0;

        return wijk;
    }

    template<typename IdxType1, typename IdxType2>
    LocalIndexType findLocalIndex_(const std::vector<IdxType1>& vector, const IdxType2 globalIdx) const
    {
        auto it = std::find(vector.begin(), vector.end(), globalIdx);
        assert(it != vector.end() && "could not find local index in the vector for the given global index!");
        return std::distance(vector.begin(), it);
    }

    const Problem& problem_() const
    { return problemRef_; }

    const FVElementGeometry& fvGeometry_() const
    { return fvGeometryRef_; }

    const ElementVolumeVariables& elemVolVars_() const
    { return elemVolVarsRef_; }

    const Problem& problemRef_;
    const FVElementGeometry& fvGeometryRef_;
    const ElementVolumeVariables& elemVolVarsRef_;

    bool onBoundary_;
    LocalIndexType eqIdx_;

    std::vector<Element> localElements_;
    std::vector<LocalScvType> localScvs_;
    std::vector<LocalScvfType> localScvfs_;

    Stencil globalScvfIndices_;
    Stencil stencil_;
    Stencil volVarsStencil_;
    std::vector<DimVector> volVarsPositions_;

    LocalIndexSet fluxFaceIndexSet_;
    LocalIndexSet dirichletFaceIndexSet_;

    Matrix T_;
    Matrix AinvB_;
    Matrix CAinv_;

    Vector neumannFluxes_;
};

} // end namespace

#endif
