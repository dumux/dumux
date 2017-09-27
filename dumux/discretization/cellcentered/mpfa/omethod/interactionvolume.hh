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
#include <dumux/discretization/cellcentered/mpfa/methods.hh>
#include <dumux/discretization/cellcentered/mpfa/dualgridindexset.hh>

#include "localsubcontrolentities.hh"
#include "interactionvolumeindexset.hh"

namespace Dumux
{
//! Specialization of the interaction volume traits class for the mpfa-o method
template<class TypeTag>
class CCMpfaOInteractionVolumeTraits : public CCMpfaInteractionVolumeTraitsBase<TypeTag>
{
    using BaseTraits = CCMpfaInteractionVolumeTraitsBase<TypeTag>;
    using NodalIndexSet = typename CCMpfaDualGridIndexSet<TypeTag>::NodalIndexSet;

public:
    using SecondaryInteractionVolume = CCMpfaInteractionVolumeImplementation<TypeTag, MpfaMethods::oMethod>;

    using typename BaseTraits::DynamicLocalIndexContainer;
    using typename BaseTraits::DynamicGlobalIndexContainer;
    using IndexSet = CCMpfaOInteractionVolumeIndexSet<NodalIndexSet, DynamicGlobalIndexContainer, DynamicLocalIndexContainer>;

    // The matrix & vector types used in the interaction
    // volume are actually the dynamic types in the o-scheme
    using typename BaseTraits::DynamicVector;
    using typename BaseTraits::DynamicMatrix;
    using Vector = DynamicVector;
    using Matrix = DynamicMatrix;

    using LocalScvType = CCMpfaOInteractionVolumeLocalScv<TypeTag, IndexSet>;
    using LocalScvfType = CCMpfaOInteractionVolumeLocalScvf<TypeTag>;
};

//! Forward declaration of the mpfa-o interaction volume
template<class TypeTag, class Traits>
class CCMpfaOInteractionVolume;

/*!
 * \ingroup Mpfa
 * \brief Base class for the interaction volumes of the mpfa-o method.
 *        We introduce one more level of inheritance here because the o-method with
 *        full pressure support uses the mpfa-o interaction volume but uses a different
 *        traits class.
 */
template<class TypeTag>
class CCMpfaInteractionVolumeImplementation<TypeTag, MpfaMethods::oMethod>
          : public CCMpfaOInteractionVolume<TypeTag,
                                            CCMpfaOInteractionVolumeTraits<TypeTag>>
{
    using TraitsType = CCMpfaOInteractionVolumeTraits<TypeTag>;
public:
    // state the traits class type
    using Traits = TraitsType;
};

template<class TypeTag, class Traits>
class CCMpfaOInteractionVolume : public CCMpfaInteractionVolumeBase<TypeTag, Traits>
{
    using ParentType = CCMpfaInteractionVolumeBase<TypeTag, Traits>;
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using MpfaHelper = typename GET_PROP_TYPE(TypeTag, MpfaHelper);
    using PrimaryVariables = typename GET_PROP_TYPE(TypeTag, PrimaryVariables);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolumeFace = typename GET_PROP_TYPE(TypeTag, SubControlVolumeFace);
    using ElementVolumeVariables = typename GET_PROP_TYPE(TypeTag, ElementVolumeVariables);
    using DualGridNodalIndexSet = typename CCMpfaDualGridIndexSet<TypeTag>::NodalIndexSet;

    static const int dim = GridView::dimension;
    static const int dimWorld = GridView::dimensionworld;
    using Element = typename GridView::template Codim<0>::Entity;
    using DimVector = Dune::FieldVector<Scalar, dim>;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

    using Vector = typename Traits::DynamicVector;
    using Matrix = typename Traits::DynamicMatrix;
    using Tensor = typename Traits::Tensor;

    using LocalScvType = typename Traits::LocalScvType;
    using LocalScvfType = typename Traits::LocalScvfType;

    using IndexSet = typename Traits::IndexSet;
    using LocalIndexContainer = typename Traits::DynamicLocalIndexContainer;
    using LocalIndexType = typename LocalIndexContainer::value_type;
    using GlobalIndexContainer = typename Traits::DynamicGlobalIndexContainer;
    using DataHandle = typename Traits::DataHandle;
public:

    using typename ParentType::GlobalLocalFaceDataPair;
    using typename ParentType::LocalFaceData;
    using typename ParentType::DirichletDataContainer;

    //! Sets up the local scope for a given seed
    //! This function has to be called before using the IV!
    void bind(const IndexSet& indexSet,
              const Problem& problem,
              const FVElementGeometry& fvGeometry,
              const ElementVolumeVariables& elemVolVars,
              DataHandle& dataHandle)
    {
        problemPtr_ = &problem;
        fvGeometryPtr_ = &fvGeometry;
        elemVolVarsPtr_ = &elemVolVars;
        indexSetPtr_ = &indexSet;

        // set up the local scope of this interaction volume
        setLocalScope_(dataHandle);
    }

    //! Sets only the pointers to the local views
    //! Using the IV afterwards requires having called bind once before!
    //! Calling this with an fvGeometry or elemVolVars that differ from the
    //! ones that have been passed when calling bind leads to undefined behaviour
    void resetPointers(const FVElementGeometry& fvGeometry, const ElementVolumeVariables& elemVolVars)
    {
        fvGeometryPtr_ = &fvGeometry;
        elemVolVarsPtr_ = &elemVolVars;

        for (unsigned i = 0; i < globalLocalScvfPairedData_.size(); ++i)
            globalLocalScvfPairedData_[i].first = &fvGeometry.scvf(globalScvfIndices_[i]);
    }

    //! solves for the transmissibilities subject to a given tensor
    template<typename GetTensorFunction>
    void solveLocalSystem(const GetTensorFunction& getTensor, DataHandle& dataHandle)
    {
        // if only dirichlet faces are present, assemble T_ directly
        if (numUnknowns_ == 0)
            assemblePureDirichletSystem_(getTensor, dataHandle.T());
        else
        {
            // assemble
            assembleLocalMatrices_(getTensor);

            // solve
            A_.invert();

            // T = C*A^-1*B + D
            auto& T = dataHandle.T();
            T = multiplyMatrices(C_.rightmultiply(A_), B_);
            T += D_;

            // store A-1B only when gradient reconstruction is necessary
            if (requireABMatrix_())
                dataHandle.AB() = B_.leftmultiply(A_);
        }

        // // set vol vars stencil & positions pointer in handle
        dataHandle.setVolVarsStencilPointer(indexSet().nodalIndexSet().globalScvIndices());
        dataHandle.setDirichletDataPointer(dirichletData_);

        // on surface grids, additionally prepare the outside transmissibilities
        if (dim < dimWorld)
            computeOutsideTransmissibilities_(dataHandle);
    }

    //! Gets the transmissibilities for a sub-control volume face within the interaction volume.
    //! specialization for dim == dimWorld
    template<int d = dim, int dw = dimWorld>
    typename std::enable_if< (d == dw), const Vector& >::type
    getTransmissibilities(const SubControlVolumeFace& scvf, const LocalFaceData& localFaceData, const DataHandle& dataHandle) const
    { return dataHandle.T()[localFaceData.localScvfIndex]; }

    //! Gets the transmissibilities for a sub-control volume face within the interaction volume.
    //! specialization for dim < dimWorld.
    template<int d = dim, int dw = dimWorld>
    typename std::enable_if< (d < dw), const Vector& >::type
    getTransmissibilities(const SubControlVolumeFace& scvf, const LocalFaceData& localFaceData, const DataHandle& dataHandle) const
    {
        // If we come from the inside, simply return tij
        if (!localFaceData.isOutside)
            return dataHandle.T()[localFaceData.localScvfIndex];
        else
            return dataHandle.outsideTij()[this->findIndexInVector(outsideScvfIndices_, scvf.index())];
    }

    //! obtain the local data object for a given global scvf
    const LocalFaceData& getLocalFaceData(const SubControlVolumeFace& scvf) const
    { return globalLocalScvfPairedData_[this->findIndexInVector(globalScvfIndices_, scvf.index())].second; }

    //! returns the vector of data pairs storing the local face data for each global scvf
    const std::vector<GlobalLocalFaceDataPair>& globalLocalScvfPairedData() const
    { return globalLocalScvfPairedData_; }

    //! returns the grid element corresponding to a given local (in the iv) scv idx
    const Element& element(const LocalIndexType localScvIdx) const
    { return elements_[localScvIdx]; }

    //! returns the local scvf entity corresponding to a given local (in the iv) scvf idx
    const LocalScvfType& localScvf(const LocalIndexType localScvfIdx) const
    { return scvfs_[localScvfIdx]; }

    //! returns the local scv entity corresponding to a given local (in the iv) scv idx
    const LocalScvType& localScv(const LocalIndexType localScvIdx) const
    { return scvs_[localScvIdx]; }

    //! returns a reference to the problem to be solved
    const Problem& problem() const
    { return *problemPtr_; }

    //! returns a reference to the fvGeometry object
    const FVElementGeometry& fvGeometry() const
    { return *fvGeometryPtr_; }

    //! returns a reference to the element volume variables
    const ElementVolumeVariables& elemVolVars() const
    { return *elemVolVarsPtr_; }

    //! returns a reference to the corresponding iv index set
    const IndexSet& indexSet() const
    { return *indexSetPtr_; }

    const DirichletDataContainer& dirichletData() const { return dirichletData_; }

    //! returns the number of interaction volumes living around a vertex
    //! the mpfa-o scheme always constructs one iv per vertex
    static std::size_t numInteractionVolumesAtVertex(const DualGridNodalIndexSet& nodalIndexSet)
    { return 1; }

    //! adds the iv index sets living around a vertex to a given container
    //! and stores the the corresponding index in a map for each scvf
    template<class IvIndexSetContainer, class ScvfIndexMap>
    static void addInteractionVolumeIndexSets(IvIndexSetContainer& ivIndexSetContainer,
                                              ScvfIndexMap& scvfIndexMap,
                                              const DualGridNodalIndexSet& nodalIndexSet)
    {
        // the global index of the iv index set that is about to be created
        const auto curGlobalIndex = ivIndexSetContainer.size();

        // make the one index set for this node
        ivIndexSetContainer.emplace_back(nodalIndexSet);

        // store the index mapping
        for (const auto scvfIdx : nodalIndexSet.globalScvfIndices())
            scvfIndexMap[scvfIdx] = curGlobalIndex;
    }

private:
    //! returns a boolean whether or not the AB matrix has to be passed to the handles
    bool requireABMatrix_() const
    {
        static const bool requireAB = GET_PARAM_FROM_GROUP(TypeTag, bool, Vtk, AddVelocity) || dim < dimWorld;
        return requireAB;
    }

    //! clears all the containers
    void reset_()
    {
        elements_.clear();
        scvs_.clear();
        scvfs_.clear();
        globalLocalScvfPairedData_.clear();
        outsideScvfIndices_.clear();
        dirichletData_.clear();
    }

    //! sets up the local scvs and scvfs etc.. using the iv index set
    void setLocalScope_(DataHandle& dataHandle)
    {
        //! clear previous data
        reset_();

        //! number of interaction-volume-local faces
        numFaces_ = indexSet().numFaces();

        //! number of interaction-volume-local (and node-local) scvs
        const auto& scvIndices = indexSet().nodalIndexSet().globalScvIndices();
        const auto numLocalScvs = scvIndices.size();

        //! number of global scvfs appearing in this interaction volume
        const auto numGlobalScvfs = indexSet().nodalIndexSet().numScvfs();

        //! reserve memory for local entities
        elements_.reserve(numLocalScvs);
        scvs_.reserve(numLocalScvs);
        scvfs_.reserve(numFaces_);
        dirichletData_.reserve(numFaces_);
        globalLocalScvfPairedData_.reserve(numGlobalScvfs);
        globalScvfIndices_.reserve(numGlobalScvfs);

        // store outside scvf idx set on surface grids
        if (dim < dimWorld)
            outsideScvfIndices_.reserve(numGlobalScvfs);

        // set up quantities related to sub-control volumes
        for (LocalIndexType scvIdxLocal = 0; scvIdxLocal < numLocalScvs; scvIdxLocal++)
        {
            const auto scvIdxGlobal = scvIndices[scvIdxLocal];
            scvs_.emplace_back(fvGeometry(), fvGeometry().scv(scvIdxGlobal), scvIdxLocal, indexSet());
            elements_.emplace_back(fvGeometry().fvGridGeometry().element(scvIdxGlobal));
        }

        // keep track of the number of unknowns etc
        numUnknowns_ = 0;
        numPotentials_ = numLocalScvs;

        // set up quantitites related to sub-control volume faces
        for (LocalIndexType faceIdxLocal = 0; faceIdxLocal < numFaces_; ++faceIdxLocal)
        {
            const auto scvfIdxGlobal = indexSet().scvfIdxGlobal(faceIdxLocal);
            const auto& neighborScvIndicesLocal = indexSet().localNeighboringScvIndices(faceIdxLocal);
            const auto insideLocalScvIdx = neighborScvIndicesLocal[0];

            // we have to use the "inside" scv face here
            const auto& scvf = fvGeometry().scvf(scvfIdxGlobal);

            // create global/local face data for this face and the global/local map
            globalLocalScvfPairedData_.emplace_back(&scvf, LocalFaceData(faceIdxLocal, insideLocalScvIdx, false));
            globalScvfIndices_.push_back(scvf.index());

            // create iv-local scvf
            if (scvf.boundary())
            {
                const auto insideElement = elements_[insideLocalScvIdx];
                const auto bcTypes = problem().boundaryTypes(insideElement, scvf);

                if (bcTypes.hasOnlyDirichlet())
                {
                    scvfs_.emplace_back(scvf, neighborScvIndicesLocal, /*isDirichlet*/true, numPotentials_++);
                    dirichletData_.emplace_back(scvf.outsideScvIdx(), scvf.ipGlobal());
                }
                else
                    scvfs_.emplace_back(scvf, neighborScvIndicesLocal, /*isDirichlet*/false, numUnknowns_++);
            }
            else
            {
                scvfs_.emplace_back(scvf, neighborScvIndicesLocal, /*isDirichlet*/false, numUnknowns_++);

                // add outside faces to the global/local face data
                for (unsigned int i = 1; i < neighborScvIndicesLocal.size(); ++i)
                {
                    const auto outsideLocalScvIdx = neighborScvIndicesLocal[i];

                    // loop over scvfs in outside scv until we find the one coinciding with current scvf
                    for (int coord = 0; coord < dim; ++coord)
                    {
                        if (indexSet().localScvfIndexInScv(outsideLocalScvIdx, coord) == faceIdxLocal)
                        {
                            const auto nodeLocalScvfIdx = indexSet().nodalIndexSet().localScvfIndicesInScv(outsideLocalScvIdx)[coord];
                            const auto globalScvfIdx = indexSet().nodalIndexSet().scvfIdxGlobal(nodeLocalScvfIdx);
                            const auto& flipScvf = fvGeometry().scvf(globalScvfIdx);
                            globalLocalScvfPairedData_.emplace_back(&flipScvf, LocalFaceData(faceIdxLocal, outsideLocalScvIdx, true));
                            globalScvfIndices_.push_back(flipScvf.index());
                            if (dim < dimWorld)
                                outsideScvfIndices_.push_back(flipScvf.index());
                        }
                    }
                }
            }
        }

        // resize the matrices in the data handle
        dataHandle.resizeT(numFaces_, numPotentials_);
        if (requireABMatrix_())
            dataHandle.resizeAB(numUnknowns_, numPotentials_);

        // resize the local matrices
        A_.resize(numUnknowns_, numUnknowns_);
        B_.resize(numUnknowns_, numPotentials_);
        C_.resize(numFaces_, numUnknowns_);
        D_.resize(numFaces_, numPotentials_);

        // on surface grids, resize the vector containing the "outside" transmissibilities
        if (dim < dimWorld)
            dataHandle.resizeOutsideTij(outsideScvfIndices_.size(), numPotentials_);
    }

    //! Assembles the local matrices that define the local system of equations and flux expressions
    template<typename GetTensorFunction>
    void assembleLocalMatrices_(const GetTensorFunction& getTensor)
    {
        // reset matrices
        A_ = 0.0;
        B_ = 0.0;
        C_ = 0.0;
        D_ = 0.0;

        // reserve space for the omegas
        wijk_.resize(scvfs_.size());

        // loop over the local faces
        for (unsigned int faceIdx = 0; faceIdx < numFaces_; ++faceIdx)
        {
            const auto& curLocalScvf = localScvf(faceIdx);
            const auto& curGlobalScvf = fvGeometry().scvf(curLocalScvf.globalScvfIndex());
            const auto curIsDirichlet = curLocalScvf.isDirichlet();
            const auto curLocalDofIdx = curLocalScvf.localDofIndex();

            // get diffusion tensor in "positive" sub volume
            const auto& neighborScvIndices = curLocalScvf.neighboringLocalScvIndices();
            const auto posLocalScvIdx = neighborScvIndices[0];
            const auto& posLocalScv = localScv(posLocalScvIdx);
            const auto& posGlobalScv = fvGeometry().scv(posLocalScv.globalScvIndex());
            const auto& posVolVars = elemVolVars()[posGlobalScv];
            const auto& posElement = element(posLocalScvIdx);
            const auto tensor = getTensor(problem(), posElement, posVolVars, fvGeometry(), posGlobalScv);

            // the omega factors of the "positive" sub volume
            auto posWijk = calculateOmegas_(posLocalScv, curGlobalScvf.unitOuterNormal(), curGlobalScvf.area(), tensor);
            posWijk *= posVolVars.extrusionFactor();

            // go over the coordinate directions in the positive sub volume
            for (unsigned int localDir = 0; localDir < dim; localDir++)
            {
                const auto otherLocalScvfIdx = posLocalScv.scvfIdxLocal(localDir);
                const auto& otherLocalScvf = localScvf(otherLocalScvfIdx);
                const auto otherLocalDofIdx = otherLocalScvf.localDofIndex();

                // if we are not on a Dirichlet face, add entries associated with unknown face pressures
                // i.e. in matrix C and maybe A (if current face is not a Dirichlet face)
                if (!otherLocalScvf.isDirichlet())
                {
                    C_[faceIdx][otherLocalDofIdx] -= posWijk[localDir];
                    if (!curIsDirichlet)
                        A_[curLocalDofIdx][otherLocalDofIdx] -= posWijk[localDir];
                }
                // the current face is a Dirichlet face and creates entries in D & maybe B
                else
                {
                    D_[faceIdx][otherLocalDofIdx] -= posWijk[localDir];
                    if (!curIsDirichlet)
                        B_[curLocalDofIdx][otherLocalDofIdx] += posWijk[localDir];
                }

                // add entries related to pressures at the scv centers (dofs)
                const auto posScvLocalDofIdx = posLocalScv.localDofIndex();
                D_[faceIdx][posScvLocalDofIdx] += posWijk[localDir];

                if (!curIsDirichlet)
                    B_[curLocalDofIdx][posScvLocalDofIdx] -= posWijk[localDir];
            }

            // store the omegas
            wijk_[faceIdx].emplace_back(std::move(posWijk));

            // If we are on an interior face, add values from negative sub volume
            if (!curGlobalScvf.boundary())
            {
                // loop over all the outside neighbors of this face and add entries
                for (unsigned int idxInOutside = 0; idxInOutside < curGlobalScvf.numOutsideScvs(); ++idxInOutside)
                {
                    const auto negLocalScvIdx = neighborScvIndices[idxInOutside+1];
                    const auto& negLocalScv = localScv(negLocalScvIdx);
                    const auto& negGlobalScv = fvGeometry().scv(negLocalScv.globalScvIndex());
                    const auto& negVolVars = elemVolVars()[negGlobalScv];
                    const auto& negElement = element(negLocalScvIdx);
                    const auto negTensor = getTensor(problem(), negElement, negVolVars, fvGeometry(), negGlobalScv);

                    // the omega factors of the "negative" sub volume
                    DimVector negWijk;

                    // if dim < dimWorld, use outside normal vector
                    if (dim < dimWorld)
                    {
                        const auto& flipScvf = fvGeometry().flipScvf(curGlobalScvf.index(), idxInOutside);
                        auto negNormal = flipScvf.unitOuterNormal();
                        negNormal *= -1.0;
                        negWijk = calculateOmegas_(negLocalScv, negNormal, curGlobalScvf.area(), negTensor);
                    }
                    else
                        negWijk = calculateOmegas_(negLocalScv, curGlobalScvf.unitOuterNormal(), curGlobalScvf.area(), negTensor);

                    // scale by extrusion factpr
                    negWijk *= negVolVars.extrusionFactor();

                    // go over the coordinate directions in the positive sub volume
                    for (int localDir = 0; localDir < dim; localDir++)
                    {
                        const auto otherLocalScvfIdx = negLocalScv.scvfIdxLocal(localDir);
                        const auto& otherLocalScvf = localScvf(otherLocalScvfIdx);
                        const auto otherLocalDofIdx = otherLocalScvf.localDofIndex();

                        if (!otherLocalScvf.isDirichlet())
                            A_[curLocalDofIdx][otherLocalDofIdx] += negWijk[localDir];
                        else
                            B_[curLocalDofIdx][otherLocalDofIdx] -= negWijk[localDir];

                        // add entries to matrix B
                        B_[curLocalDofIdx][negLocalScv.localDofIndex()] += negWijk[localDir];
                    }

                    // store the omegas
                    wijk_[faceIdx].emplace_back(std::move(negWijk));
                }
            }
        }
    }

    //! for interaction volumes that have only dirichlet scvfs,
    //! the transmissibility matrix can be assembled directly
    template<typename GetTensorFunction>
    void assemblePureDirichletSystem_(const GetTensorFunction& getTensor, Matrix& T)
    {
        // reset the transmissibility matrix beforehand
        T = 0.0;

        // Loop over all the faces, in this case these are all dirichlet boundaries
        for (unsigned int faceIdx = 0; faceIdx < numFaces_; ++faceIdx)
        {
            const auto& curLocalScvf = localScvf(faceIdx);
            const auto& curGlobalScvf = fvGeometry().scvf(curLocalScvf.globalScvfIndex());

            // get diffusion tensor in "positive" sub volume
            const auto& neighborScvIndices = curLocalScvf.neighboringLocalScvIndices();
            const auto posLocalScvIdx = neighborScvIndices[0];
            const auto& posLocalScv = localScv(posLocalScvIdx);
            const auto& posGlobalScv = fvGeometry().scv(posLocalScv.globalScvIndex());
            const auto& posVolVars = elemVolVars()[posGlobalScv];
            const auto& posElement = element(posLocalScvIdx);
            const auto tensor = getTensor(problem(), posElement, posVolVars, fvGeometry(), posGlobalScv);

            // the omega factors of the "positive" sub volume
            auto posWijk = calculateOmegas_(posLocalScv, curGlobalScvf.unitOuterNormal(), curGlobalScvf.area(), tensor);
            posWijk *= posVolVars.extrusionFactor();

            const auto posScvLocalDofIdx = posLocalScv.localDofIndex();
            for (int localDir = 0; localDir < dim; localDir++)
            {
                const auto otherLocalScvfIdx = posLocalScv.scvfIdxLocal(localDir);
                const auto& otherLocalScvf = localScvf(otherLocalScvfIdx);
                const auto otherLocalDofIdx = otherLocalScvf.localDofIndex();
                T[faceIdx][otherLocalDofIdx] -= posWijk[localDir];
                T[faceIdx][posScvLocalDofIdx] += posWijk[localDir];
            }
        }
    }

    //! computes the transmissibilities associated with "outside" faces on surface grids
    void computeOutsideTransmissibilities_(DataHandle& dataHandle) const
    {
        assert(dim < dimWorld && "transmissibilities for outside scvfs can only be computed for dim < dimWorld!");

        for (const auto& globalLocalData : globalLocalScvfPairedData_)
        {
            const auto& localFaceData = globalLocalData.second;

            // continue only for "outside" faces
            if (!localFaceData.isOutside)
                continue;

            const auto localScvIdx = localFaceData.localScvIndex;
            const auto localScvfIdx = localFaceData.localScvfIndex;
            const auto& posLocalScv = localScv(localScvIdx);

            const auto idxInNeighbors = this->findIndexInVector(localScvf(localScvfIdx).neighboringLocalScvIndices(), localScvIdx);
            const auto& wijk = wijk_[localScvfIdx][idxInNeighbors];

            // store the calculated transmissibilities in the data handle
            const auto idxInOutsideFaces = this->findIndexInVector(outsideScvfIndices_, globalLocalData.first->index());
            auto& tij = dataHandle.outsideTij()[idxInOutsideFaces];
            tij = 0.0;

            for (int localDir = 0; localDir < dim; localDir++)
            {
                const auto curLocalScvfIdx = posLocalScv.scvfIdxLocal(localDir);
                const auto& curLocalScvf = localScvf(curLocalScvfIdx);
                if (!curLocalScvf.isDirichlet())
                {
                    auto tmp = dataHandle.AB()[curLocalScvf.localDofIndex()];
                    tmp *= wijk[localDir];
                    tij -= tmp;
                }
                else
                    tij[curLocalScvf.localDofIndex()] -= wijk[localDir];

                // add entry from the scv unknown
                tij[localScvIdx] += wijk[localDir];
            }
        }
    }

    // calculates n_i^T*K_j*nu_k
    DimVector calculateOmegas_(const LocalScvType& localScv,
                               const GlobalPosition normal,
                               const Scalar area,
                               const Tensor& T) const
    {
        // make sure we have positive definite diffsion tensors
        assert(this->tensorIsPositiveDefinite(T) && "only positive definite tensors can be handled by mpfa methods");

        DimVector wijk;
        GlobalPosition tmp;
        for (int dir = 0; dir < dim; ++dir)
        {
            T.mv(localScv.innerNormal(dir), tmp);
            wijk[dir] = tmp*normal;
        }
        wijk *= area;
        wijk /= localScv.detX();

        return wijk;
    }

    // calculates n_i^T*K_j*nu_k
    DimVector calculateOmegas_(const LocalScvType& localScv,
                               const GlobalPosition normal,
                               const Scalar area,
                               const Scalar t) const
    {
        // make sure we have positive diffusion coefficients
        assert(t > 0.0 && "non-positive diffusion coefficients cannot be handled by mpfa methods");

        DimVector wijk;
        GlobalPosition tmp(normal);
        tmp *= t;

        for (int dir = 0; dir < dim; ++dir)
            wijk[dir] = tmp*localScv.innerNormal(dir);
        wijk *= area;
        wijk /= localScv.detX();

        return wijk;
    }

    const Problem* problemPtr_;
    const FVElementGeometry* fvGeometryPtr_;
    const ElementVolumeVariables* elemVolVarsPtr_;
    const IndexSet* indexSetPtr_;

    // Variables defining the local scope
    std::vector<Element> elements_;
    std::vector<LocalScvType> scvs_;
    std::vector<LocalScvfType> scvfs_;
    std::vector<GlobalLocalFaceDataPair> globalLocalScvfPairedData_;
    DirichletDataContainer dirichletData_;
    GlobalIndexContainer globalScvfIndices_;
    GlobalIndexContainer outsideScvfIndices_;

    // sizes involved in the local matrices
    unsigned int numFaces_;
    unsigned int numUnknowns_;
    unsigned int numPotentials_;

    // The omega factors and the matrix A⁻1*B are stored
    // in order to recover the transmissibilities of outside faces on network grids
    std::vector< std::vector< DimVector > > wijk_;
    Matrix CAinv_;

    // Matrices involved in transmissibility calculations
    Matrix A_;
    Matrix B_;
    Matrix C_;
    Matrix D_;
};

} // end namespace

#endif
