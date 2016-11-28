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
#include <dumux/discretization/cellcentered/mpfa/interactionvolumebase.hh>
#include <dumux/discretization/cellcentered/mpfa/facetypes.hh>
#include <dumux/discretization/cellcentered/mpfa/methods.hh>

#include "interactionvolumeseed.hh"
#include "localsubcontrolentities.hh"

namespace Dumux
{
//! Specialization of the interaction volume traits class for the mpfa-o method
template<class TypeTag>
class CCMpfaOInteractionVolumeTraits : public CCMpfaInteractionVolumeTraitsBase<TypeTag>
{
    using BaseTraits = CCMpfaInteractionVolumeTraitsBase<TypeTag>;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);

    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

public:
    using BoundaryInteractionVolume = CCMpfaInteractionVolumeImplementation<TypeTag, MpfaMethods::oMethod>;

    using PositionVector = std::vector<GlobalPosition>;
    using Matrix = Dune::DynamicMatrix<Scalar>;
    using Vector = typename Matrix::row_type;

    using LocalScvType = CCMpfaOLocalScv<TypeTag>;
    using LocalScvfType = CCMpfaOLocalScvf<TypeTag>;
    using typename BaseTraits::LocalIndexSet;
    using typename BaseTraits::GlobalIndexSet;
    using Seed = CCMpfaOInteractionVolumeSeed<GlobalIndexSet, LocalIndexSet, dim, dimWorld>;
};

//! Forward declaration of the mpfa-o interaction volume
template<class TypeTag, class Traits> class CCMpfaOInteractionVolume;

/*!
 * \ingroup Mpfa
 * \brief Base class for the interaction volumes of the mpfa-o method.
 *        We introduce one more level of inheritance here because the o-method with
 *        full pressure support uses the mpfa-o interaction volume but uses a different
 *        traits class.
 */
template<class TypeTag>
class CCMpfaInteractionVolumeImplementation<TypeTag, MpfaMethods::oMethod> : public CCMpfaOInteractionVolume<TypeTag, CCMpfaOInteractionVolumeTraits<TypeTag>>
{
    using Traits = CCMpfaOInteractionVolumeTraits<TypeTag>;
    using ParentType = CCMpfaOInteractionVolume<TypeTag, Traits>;

    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);

    using IVSeed = typename Traits::Seed;
public:
    CCMpfaInteractionVolumeImplementation(const IVSeed& seed,
                                          const Problem& problem,
                                          const FVElementGeometry& fvGeometry,
                                          const ElementVolumeVariables& elemVolVars)
    : ParentType(seed, problem, fvGeometry, elemVolVars)
    {}

};

template<class TypeTag, class Traits>
class CCMpfaOInteractionVolume : public CCMpfaInteractionVolumeBase<TypeTag, Traits>
{
    // The interaction volume implementation has to be friend,
    // because some methods use the mpfa-o interaction volume as base
    friend typename GET_PROP_TYPE(TypeTag, InteractionVolume);
    using ParentType = CCMpfaInteractionVolumeBase<TypeTag, Traits>;

    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using MpfaHelper = typename GET_PROP_TYPE(TypeTag, MpfaHelper);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);

    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;
    using Element = typename GridView::template Codim<0>::Entity;
    using DimVector = Dune::FieldVector<Scalar, dim>;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

    using DynamicVector = typename Traits::Vector;
    using DynamicMatrix = typename Traits::Matrix;
    using Tensor = typename Traits::Tensor;

    using LocalScvType = typename Traits::LocalScvType;
    using LocalScvfType = typename Traits::LocalScvfType;

public:
    using typename ParentType::LocalIndexType;
    using typename ParentType::LocalIndexSet;
    using typename ParentType::LocalFaceData;
    using typename ParentType::GlobalIndexSet;
    using typename ParentType::PositionVector;
    using typename ParentType::Seed;

    CCMpfaOInteractionVolume(const Seed& seed,
                             const Problem& problem,
                             const FVElementGeometry& fvGeometry,
                             const ElementVolumeVariables& elemVolVars)
    : problemPtr_(&problem),
      fvGeometryPtr_(&fvGeometry),
      elemVolVarsPtr_(&elemVolVars),
      onBoundary_(seed.onBoundary()),
      globalScvfIndices_(seed.globalScvfIndices())
    {
        // create local sub control entities from the seed
        createLocalEntities_(seed);

        // reserve memory
        const auto numLocalScvs = localScvs_.size();
        const auto numLocalScvfs = localScvfs_.size();
        const auto maxNumVolVars = numLocalScvs + numLocalScvfs;
        volVarsStencil_ = seed.globalScvIndices(); // boundary vol vars are placed at the end
        volVarsStencil_.reserve(maxNumVolVars);
        volVarsPositions_.reserve(maxNumVolVars);
        dirichletFaceIndexSet_.reserve(numLocalScvfs);
        fluxFaceIndexSet_.reserve(numLocalScvfs);

        // the positions where the vol vars are defined at (required for the gravitational acceleration)
        for (const auto& localScv : localScvs_)
            volVarsPositions_.push_back(localScv.center());

        // eventually add dirichlet vol var indices and set up local index sets of flux and dirichlet faces
        LocalIndexType localScvfIdx = 0;
        for (const auto& localScvf : localScvfs_)
        {
            // eventually add vol var index and corresponding position
            if (localScvf.faceType() == MpfaFaceTypes::dirichlet)
            {
                volVarsStencil_.push_back(localScvf.outsideGlobalScvIndex());
                volVarsPositions_.push_back(localScvf.ip());
                dirichletFaceIndexSet_.push_back(localScvfIdx++);
            }
            else
                fluxFaceIndexSet_.push_back(localScvfIdx++);
        }

        // initialize the neumann fluxes vector to zero
        neumannFluxes_ = DynamicVector(fluxFaceIndexSet_.size(), 0.0);
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
        DynamicMatrix A(numFluxFaces, numFluxFaces, 0.0);
        DynamicMatrix B(numFluxFaces, numPotentials, 0.0);
        DynamicMatrix C(numFaces, numFluxFaces, 0.0);
        DynamicMatrix D(numFaces, numPotentials, 0.0);

        assembleLocalMatrices_(getTensor, A, B, C, D);

        // solve local system and store matrices
        DynamicMatrix copy(B);
        A.invert();
        AinvB_ = B.leftmultiply(A);
        CAinv_ = C.rightmultiply(A);
        T_ = multiplyMatrices(CAinv_, copy);
        T_ += D;
    }

    template<typename UpwindFactorFunction>
    void assembleNeumannFluxes(const UpwindFactorFunction& upwindFactor, const unsigned int eqIdx)
    {
        if (!onBoundary() || GET_PROP_VALUE(TypeTag, UseTpfaBoundary))
            return;

        LocalIndexType fluxFaceIdx = 0;
        for (const auto localFluxFaceIdx : fluxFaceIndexSet_)
        {
            const auto& localScvf = localScvf_(localFluxFaceIdx);
            const auto faceType = localScvf.faceType();
            if (faceType == MpfaFaceTypes::neumann)
            {
                const auto& element = localElement_(localScvf.insideLocalScvIndex());
                const auto& globalScvf = fvGeometry_().scvf(localScvf.insideGlobalScvfIndex());
                auto neumannFlux = problem_().neumann(element, this->fvGeometry_(), this->elemVolVars_(), globalScvf)[eqIdx];
                neumannFlux *= globalScvf.area();

                // recover -k*gradh
                const auto& insideScv = fvGeometry_().scv(globalScvf.insideScvIdx());
                const auto& volVars = elemVolVars_()[insideScv];
                neumannFlux /= upwindFactor(volVars);

                neumannFluxes_[fluxFaceIdx] = neumannFlux;
            }

            fluxFaceIdx++;
        }
    }

    //! gets local data on a global scvf within the interaction volume
    //! specialization for dim = dimWorld
    template<int d = dim, int dw = dimWorld>
    typename std::enable_if< (d == dw), LocalFaceData >::type
    getLocalFaceData(const SubControlVolumeFace& scvf) const
    {
        auto scvfGlobalIdx = scvf.index();

        LocalIndexType localIdx = 0;
        for (const auto& localScvf : localScvfs_)
        {
            if (localScvf.insideGlobalScvfIndex() == scvfGlobalIdx)
                return LocalFaceData(localIdx, localScvf.insideLocalScvIndex(), false);

            if (!localScvf.boundary() && localScvf.outsideGlobalScvfIndex() == scvfGlobalIdx)
                return LocalFaceData(localIdx, localScvf.outsideLocalScvIndex(), true);

            localIdx++;
        }

        DUNE_THROW(Dune::InvalidStateException, "Could not find the local scv face in the interaction volume for the given scvf with index: " << scvf.index());
    }

    //! gets local data on a global scvf within the interaction volume
    //! specialization for dim < dimWorld
    template<int d = dim, int dw = dimWorld>
    typename std::enable_if< (d < dw), LocalFaceData >::type
    getLocalFaceData(const SubControlVolumeFace& scvf) const
    {
        const auto scvfGlobalIdx = scvf.index();

        LocalIndexType localFaceIdx = 0;
        for (const auto& localScvf : localScvfs_)
        {
            if (localScvf.insideGlobalScvfIndex() == scvfGlobalIdx)
                return LocalFaceData(localFaceIdx, localScvf.insideLocalScvIndex(), false);

            if (!localScvf.boundary() && MpfaHelper::contains(localScvf.outsideGlobalScvfIndices(), scvfGlobalIdx))
            {
                if (localScvf.outsideLocalScvIndices().size() == 1)
                    return LocalFaceData(localFaceIdx, localScvf.outsideLocalScvIndex(), true);

                // Now we have to find wich local scv is the one the global scvf is embedded in
                const auto scvGlobalIdx = scvf.insideScvIdx();
                for (auto localScvIdx : localScvf.outsideLocalScvIndices())
                    if (localScv_(localScvIdx).globalIndex() == scvGlobalIdx)
                        return LocalFaceData(localFaceIdx, localScvIdx, true);

                DUNE_THROW(Dune::InvalidStateException, "Could not find the local scv in the interaction volume for the given scvf with index: " << scvfGlobalIdx);
            }

            localFaceIdx++;
        }

        DUNE_THROW(Dune::InvalidStateException, "Could not find the local scv face in the interaction volume for the given scvf with index: " << scvfGlobalIdx);
    }

    //! gets the transmissibilities for a sub-control volume face within the interaction volume
    //! specialization for dim == dimWorld
    template<int d = dim, int dw = dimWorld>
    typename std::enable_if< (d == dw), DynamicVector >::type
    getTransmissibilities(const LocalFaceData& localFaceData) const
    {
        auto tij = T_[localFaceData.localScvfIndex];
        if (localFaceData.isOutside)
            tij *= -1.0;
        return tij;
    }

    //! gets the transmissibilities for a sub-control volume face within the interaction volume
    //! specialization for dim < dimWorld
    template<int d = dim, int dw = dimWorld>
    typename std::enable_if< (d < dw), DynamicVector >::type
    getTransmissibilities(const LocalFaceData& localFaceData) const
    {
        // if we are not on a branching point, do the usual
        if (localScvf_(localFaceData.localScvfIndex).outsideLocalScvIndices().size() == 1)
        {
            auto tij = T_[localFaceData.localScvfIndex];
            if (localFaceData.isOutside)
                tij *= -1.0;
            return tij;
        }

        // This means we are on a branching point. If we come from the inside, simply return tij
        if (!localFaceData.isOutside)
            return T_[localFaceData.localScvfIndex];

        // otherwise compute the outside transmissibilities
        DynamicVector tij(volVarsStencil().size(), 0.0);

        // get the local scv and iterate over local coordinates
        const std::size_t numLocalScvs = localScvs_.size();
        const auto& localScv = localScv_(localFaceData.localScvIndex);
        const auto& localScvf = localScvf_(localFaceData.localScvfIndex);

        auto idxInOutside = this->findLocalIndex(localScvf.outsideLocalScvIndices(), localFaceData.localScvIndex);
        const auto& wijk = wijk_[localFaceData.localScvfIndex][idxInOutside+1];
        for (int localDir = 0; localDir < dim; localDir++)
        {
            auto localScvfIdx = localScv.localScvfIndex(localDir);
            const auto& localScvf = localScvf_(localScvfIdx);

            // add entries from the face unknowns
            if (localScvf.faceType() != MpfaFaceTypes::dirichlet)
            {
                auto fluxFaceIndex = this->findLocalIndex(fluxScvfIndexSet_(), localScvfIdx);
                auto tmp = AinvB_[fluxFaceIndex];
                tmp *= wijk[localDir];

                tij += tmp;
            }
            else
            {
                auto idxInDiriFaces = this->findLocalIndex(dirichletScvfIndexSet_(), localScvfIdx);
                tij[numLocalScvs + idxInDiriFaces] += wijk[localDir];
            }

            // add entry from the scv unknown
            tij[localFaceData.localScvIndex] -= wijk[localDir];
        }

        return tij;
    }

    Scalar getNeumannFlux(const LocalFaceData& localFaceData) const
    {
        if (!onBoundary() || GET_PROP_VALUE(TypeTag, UseTpfaBoundary) || fluxScvfIndexSet_().size() == 0 )
            return 0.0;

        auto flux = CAinv_[localFaceData.localScvfIndex] * neumannFluxes_;
        if (localFaceData.isOutside)
            return -1.0*flux;
        return flux;
    }

    bool onBoundary() const
    { return onBoundary_; }

    const GlobalIndexSet& volVarsStencil() const
    { return volVarsStencil_; }

    const PositionVector& volVarsPositions() const
    { return volVarsPositions_; }

    const GlobalIndexSet& globalScvfs() const
    { return globalScvfIndices_; }

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
            localScvs_.emplace_back(LocalScvType(problem_(), element, fvGeometry_(), scvSeed));
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
                                DynamicMatrix& A,
                                DynamicMatrix& B,
                                DynamicMatrix& C,
                                DynamicMatrix& D)
    {
        const std::size_t numLocalScvs = localScvs_.size();

        // reserve space for the omegas
        wijk_.resize(localScvfs_.size());

        // loop over the local faces
        LocalIndexType rowIdx = 0;
        for (const auto& localScvf : localScvfs_)
        {
            auto faceType = localScvf.faceType();
            bool hasUnknown = faceType != MpfaFaceTypes::dirichlet;
            LocalIndexType idxInFluxFaces = hasUnknown ? this->findLocalIndex(fluxScvfIndexSet_(), rowIdx) : -1;

            // get diffusion tensor in "positive" sub volume
            const auto posLocalScvIdx = localScvf.insideLocalScvIndex();
            const auto& posLocalScv = localScv_(posLocalScvIdx);
            const auto& posGlobalScv = fvGeometry_().scv(posLocalScv.globalIndex());
            auto element = localElement_(posLocalScvIdx);
            auto tensor = getTensor(element, elemVolVars_()[posGlobalScv], posGlobalScv);

            // the omega factors of the "positive" sub volume
            auto posWijk = calculateOmegas_(posLocalScv, localScvf.unitOuterNormal(), localScvf.area(), tensor);
            posWijk *= problem_().boxExtrusionFactor(element, posGlobalScv);

            // Check the local directions of the positive sub volume
            for (int localDir = 0; localDir < dim; localDir++)
            {
                auto curLocalScvfIdx = posLocalScv.localScvfIndex(localDir);
                const auto& curLocalScvf = localScvf_(curLocalScvfIdx);

                // First, add the entries associated with face pressures (unkown or dirichlet)
                if (curLocalScvf.faceType() != MpfaFaceTypes::dirichlet)
                {
                    // we need the index of the current local scvf in the flux face indices
                    auto curIdxInFluxFaces = this->findLocalIndex(fluxScvfIndexSet_(), curLocalScvfIdx);

                    C[rowIdx][curIdxInFluxFaces] += posWijk[localDir];
                    if (hasUnknown)
                        A[idxInFluxFaces][curIdxInFluxFaces] += posWijk[localDir];
                }
                else
                {
                    // the current face is a Dirichlet face and creates entries in D & eventually B
                    auto curIdxInDiriFaces = this->findLocalIndex(dirichletScvfIndexSet_(), curLocalScvfIdx);

                    D[rowIdx][numLocalScvs + curIdxInDiriFaces] += posWijk[localDir];
                    if (hasUnknown)
                        B[idxInFluxFaces][numLocalScvs + curIdxInDiriFaces] -= posWijk[localDir];
                }

                // add entries related to pressures at the scv centers (dofs)
                D[rowIdx][posLocalScvIdx] -= posWijk[localDir];

                if (hasUnknown)
                    B[idxInFluxFaces][posLocalScvIdx] += posWijk[localDir];
            }

            // store the omegas
            wijk_[rowIdx].emplace_back(std::move(posWijk));

            // If face is not on a boundary or an interior dirichlet face, add entries for the "negative" scv
            if (faceType == MpfaFaceTypes::interior)
            {
                // loop over all the outside neighbors of this face and add entries
                unsigned int indexInOutsideData = 0;
                for (auto negLocalScvIdx : localScvf.outsideLocalScvIndices())
                {
                    const auto& negLocalScv = localScv_(negLocalScvIdx);
                    const auto& negGlobalScv = fvGeometry_().scv(negLocalScv.globalIndex());
                    auto negElement = localElement_(negLocalScvIdx);
                    auto negTensor = getTensor(negElement, elemVolVars_()[negGlobalScv], negGlobalScv);

                    // the omega factors of the "negative" sub volume
                    DimVector negWijk;

                    // if dim < dimWorld, use outside normal vector
                    if (dim < dimWorld)
                    {
                        // outside scvf
                        auto&& outsideScvf = fvGeometry_().scvf(localScvf.outsideGlobalScvfIndices()[indexInOutsideData]);
                        auto negNormal = outsideScvf.unitOuterNormal();
                        negNormal *= -1.0;
                        negWijk = calculateOmegas_(negLocalScv, negNormal, localScvf.area(), negTensor);
                    }
                    else
                        negWijk = calculateOmegas_(negLocalScv, localScvf.unitOuterNormal(), localScvf.area(), negTensor);

                    // scale by extrusion factpr
                    negWijk *= problem_().boxExtrusionFactor(negElement, negGlobalScv);

                    // Check local directions of negative sub volume
                    for (int localDir = 0; localDir < dim; localDir++)
                    {
                        auto curLocalScvfIdx = negLocalScv.localScvfIndex(localDir);
                        const auto& curLocalScvf = localScvf_(curLocalScvfIdx);

                        if (curLocalScvf.faceType() != MpfaFaceTypes::dirichlet)
                        {
                            // we need the index of the current local scvf in the flux face indices
                            auto curIdxInFluxFaces = this->findLocalIndex(fluxScvfIndexSet_(), curLocalScvfIdx);
                            A[idxInFluxFaces][curIdxInFluxFaces] -= negWijk[localDir];
                        }
                        else
                        {
                            // the current face is a Dirichlet face and creates entries in B
                            auto curIdxInDiriFaces = this->findLocalIndex(dirichletScvfIndexSet_(), curLocalScvfIdx);
                            B[idxInFluxFaces][numLocalScvs + curIdxInDiriFaces] += negWijk[localDir];
                        }

                        // add entries to matrix B
                        B[idxInFluxFaces][negLocalScvIdx] -= negWijk[localDir];
                    }

                    // store the omegas (negative, because of normal vector sign switch)
                    negWijk *= -1.0;
                    wijk_[rowIdx].emplace_back(std::move(negWijk));

                    // increment counter in outside data
                    indexInOutsideData++;
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

        // resize the omegas
        wijk_.resize(numFaces);

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
            auto posWijk = calculateOmegas_(posLocalScv, localScvf.unitOuterNormal(), localScvf.area(), tensor);
            posWijk *= problem_().boxExtrusionFactor(element, posGlobalScv);

            for (int localDir = 0; localDir < dim; localDir++)
            {
                auto curLocalScvfIdx = posLocalScv.localScvfIndex(localDir);
                auto curIdxInDiriFaces = this->findLocalIndex(dirichletScvfIndexSet_(), curLocalScvfIdx);

                T_[rowIdx][numLocalScvs + curIdxInDiriFaces] += posWijk[localDir];
                T_[rowIdx][posLocalScvIdx] -= posWijk[localDir];
            }

            // store the omegas
            wijk_[rowIdx].emplace_back(std::move(posWijk));

            // go to the next face
            rowIdx++;
        }
    }

    DimVector calculateOmegas_(const LocalScvType& localScv,
                               const GlobalPosition normal,
                               const Scalar area,
                               const Tensor& T) const
    {
        DimVector wijk;
        GlobalPosition tmp;
        for (int dir = 0; dir < dim; ++dir)
        {
            T.mv(localScv.innerNormal(dir), tmp);
            wijk[dir] = tmp*normal;
        }
        wijk *= area;
        wijk /= localScv.detX();
        wijk *= -1.0;

        return wijk;
    }

    DimVector calculateOmegas_(const LocalScvType& localScv,
                               const GlobalPosition normal,
                               const Scalar area,
                               const Scalar t) const
    {
        DimVector wijk;
        GlobalPosition tmp(normal);
        tmp *= t;

        for (int dir = 0; dir < dim; ++dir)
            wijk[dir] = tmp*localScv.innerNormal(dir);
        wijk *= area;
        wijk /= localScv.detX();
        wijk *= -1.0;

        return wijk;
    }

    const Problem& problem_() const
    { return *problemPtr_; }

    const FVElementGeometry& fvGeometry_() const
    { return *fvGeometryPtr_; }

    const ElementVolumeVariables& elemVolVars_() const
    { return *elemVolVarsPtr_; }

    const Problem* problemPtr_;
    const FVElementGeometry* fvGeometryPtr_;
    const ElementVolumeVariables* elemVolVarsPtr_;

    bool onBoundary_;

    std::vector<Element> localElements_;
    std::vector<LocalScvType> localScvs_;
    std::vector<LocalScvfType> localScvfs_;

    GlobalIndexSet globalScvfIndices_;
    GlobalIndexSet volVarsStencil_;
    PositionVector volVarsPositions_;

    LocalIndexSet fluxFaceIndexSet_;
    LocalIndexSet dirichletFaceIndexSet_;

    std::vector< std::vector< DimVector > > wijk_;

    DynamicMatrix T_;
    DynamicMatrix AinvB_;
    DynamicMatrix CAinv_;

    DynamicVector neumannFluxes_;
};

} // end namespace

#endif
