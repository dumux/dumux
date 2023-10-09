// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup CCMpfaDiscretization
 * \brief Assembly functions for the local systems of equations
 *        involved in the transmissibility computations for the mpfa-o scheme.
 */
#ifndef DUMUX_DISCRETIZATION_CC_MPFA_O_IV_ASSEMBLER_HH
#define DUMUX_DISCRETIZATION_CC_MPFA_O_IV_ASSEMBLER_HH

#include <algorithm>

#include <dumux/discretization/extrusion.hh>
#include <dumux/discretization/cellcentered/mpfa/localassemblerhelper.hh>
#include <dumux/discretization/cellcentered/mpfa/computetransmissibility.hh>
#include <dumux/discretization/cellcentered/mpfa/dualgridindexset.hh>

namespace Dumux {

#ifndef DOXYGEN
namespace Detail {

/*!
 * \brief Assemble the matrices involved in the flux expressions
 *        across the scvfs inside an interaction volume as well as those involved
 *        in the interaction volume-local system of equations resulting from flux
 *        and solution continuity across the scvfs.
 *
 *        Flux expressions: \f$\mathbf{f} = \mathbf{C} \bar{\mathbf{u}} + \mathbf{D} \mathbf{u}\f$.
 *        Continuity equations: \f$\mathbf{A} \, \bar{\mathbf{u}} = \mathbf{B} \, \mathbf{u}\f$.
 *        From this we find: \f$\mathbf{f} = (\mathbf{C}\mathbf{A}^{-1}\mathbf{B} + \mathbf{D} ) \mathbf{u}\f$.
 *
 * \note  The matrices are expected to have been resized beforehand.
 *
 * \tparam IV The interaction volume type implementation
 * \tparam TensorFunc Lambda to obtain the tensor w.r.t. which the local system is to be solved
 *
 * \param A The A matrix of the iv-local equation system
 * \param B The B matrix of the iv-local equation system
 * \param C The C matrix of the iv-local flux expressions
 * \param D The D matrix of the iv-local flux expressions
 * \param iv The mpfa-o interaction volume
 * \param getT Lambda to evaluate the scv-wise tensors
 * \param wijZeroThresh Threshold below which transmissibilities are taken to be zero.
 *                      On the basis of this threshold, trivial (0 = 0) rows in the
 *                      A matrix are identified and modified accordingly in order to
 *                      avoid ending up with singular matrices. This can occur when the
 *                      tensor is zero in some cells.
 * \returns std::vector of pairs storing the local dof and the local face indices of those
 *          non-Dirichlet faces, for which the flux continuity condition results in a
 *          trivial 0 = 0 equation which could happen for zero tensors.
 */
template<class FVElementGeometry, class OmegaStorage, class IV, class TensorFunc, class Scalar>
auto assembleLocalMatrices(const FVElementGeometry& fvGeometry,
                           typename IV::Traits::MatVecTraits::AMatrix& A,
                           typename IV::Traits::MatVecTraits::BMatrix& B,
                           typename IV::Traits::MatVecTraits::CMatrix& C,
                           typename IV::Traits::MatVecTraits::DMatrix& D,
                           OmegaStorage& wijk,
                           IV& iv, const TensorFunc& getT,
                           Scalar wijZeroThresh)
{
    const auto computeOmega = [&] (const auto& scv, const auto& scvf, const auto& tensorVariant) {
        return std::visit([&] (const auto& tensor) {
            return computeMpfaTransmissibility(fvGeometry, scv, scvf, tensor, 1.0 /*todo: extrusion fac*/);
        }, tensorVariant);
    };

    static constexpr int dim = IV::Traits::GridView::dimension;
    static constexpr int dimWorld = IV::Traits::GridView::dimensionworld;

    std::vector< std::pair<unsigned int, unsigned int> > faceMarkers;
    InteractionVolumeAssemblerHelper::resizeVector(wijk, iv.numFaces());

    // if only Dirichlet faces are present in the iv,
    // the matrices A, B & C are undefined and D = T
    if (iv.numUnknowns() == 0)
    {
        // resize & reset D matrix
        InteractionVolumeAssemblerHelper::resizeMatrix(D, iv.numFaces(), iv.numKnowns()); D = 0.0;

        // Loop over all the faces, in this case these are all dirichlet boundaries
        for (unsigned int faceIdx = 0; faceIdx < iv.numFaces(); ++faceIdx)
        {
            const auto& curLocalScvf = iv.localScvf(faceIdx);
            const auto& curGlobalScvf = fvGeometry.scvf(curLocalScvf.gridScvfIndex());
            const auto& neighborScvIndices = curLocalScvf.neighboringLocalScvIndices();

            // get tensor in "positive" sub volume
            const auto& posLocalScv = iv.localScv(neighborScvIndices[0]);
            const auto& posScv = fvGeometry.scv(posLocalScv.gridScvIndex());
            const auto tensor = getT(posScv);

            // the omega factors of the "positive" sub volume
            InteractionVolumeAssemblerHelper::resizeVector(wijk[faceIdx], /*no outside scvs present*/1);
            wijk[faceIdx][0] = computeOmega(posLocalScv, curGlobalScvf, tensor);

            const auto posScvLocalDofIdx = posLocalScv.localDofIndex();
            for (unsigned int localDir = 0; localDir < dim; localDir++)
            {
                const auto& otherLocalScvf = iv.localScvf( posLocalScv.localScvfIndex(localDir) );
                const auto otherLocalDofIdx = otherLocalScvf.localDofIndex();
                D[faceIdx][otherLocalDofIdx] -= wijk[faceIdx][0][localDir];
                D[faceIdx][posScvLocalDofIdx] += wijk[faceIdx][0][localDir];
            }
        }
    }
    else
    {
        // resize & reset matrices
        InteractionVolumeAssemblerHelper::resizeMatrix(A, iv.numUnknowns(), iv.numUnknowns()); A = 0.0;
        InteractionVolumeAssemblerHelper::resizeMatrix(B, iv.numUnknowns(), iv.numKnowns());   B = 0.0;
        InteractionVolumeAssemblerHelper::resizeMatrix(C, iv.numFaces(), iv.numUnknowns());    C = 0.0;
        InteractionVolumeAssemblerHelper::resizeMatrix(D, iv.numFaces(), iv.numKnowns());      D = 0.0;

        for (unsigned int faceIdx = 0; faceIdx < iv.numFaces(); ++faceIdx)
        {
            const auto& curLocalScvf = iv.localScvf(faceIdx);
            const auto& curGlobalScvf = fvGeometry.scvf(curLocalScvf.gridScvfIndex());
            const auto curIsDirichlet = curLocalScvf.isDirichlet();
            const auto curLocalDofIdx = curLocalScvf.localDofIndex();

            // get tensor in "positive" sub volume
            const auto& neighborScvIndices = curLocalScvf.neighboringLocalScvIndices();
            const auto& posLocalScv = iv.localScv(neighborScvIndices[0]);
            const auto& posScv = fvGeometry.scv(posLocalScv.gridScvIndex());
            const auto tensor = getT(posScv);

            // the omega factors of the "positive" sub volume
            InteractionVolumeAssemblerHelper::resizeVector(wijk[faceIdx], neighborScvIndices.size());
            wijk[faceIdx][0] = computeOmega(posLocalScv, curGlobalScvf, tensor);

            using std::abs;
            bool insideZeroWij = false;

            // go over the coordinate directions in the positive sub volume
            for (unsigned int localDir = 0; localDir < dim; localDir++)
            {
                const auto& otherLocalScvf = iv.localScvf( posLocalScv.localScvfIndex(localDir) );
                const auto otherLocalDofIdx = otherLocalScvf.localDofIndex();

                if (otherLocalDofIdx == curLocalDofIdx)
                {
                    if (abs(wijk[faceIdx][0][localDir]) <= wijZeroThresh)
                    {
                        if (!curIsDirichlet)
                        {
                            insideZeroWij = true;
                            faceMarkers.emplace_back( std::make_pair(curLocalDofIdx, faceIdx) );
                        }
                    }
                }

                // if we are not on a Dirichlet face, add entries associated with unknown face pressures
                // i.e. in matrix C and maybe A (if current face is not a Dirichlet face)
                if (!otherLocalScvf.isDirichlet())
                {
                    C[faceIdx][otherLocalDofIdx] -= wijk[faceIdx][0][localDir];
                    if (!curIsDirichlet)
                        A[curLocalDofIdx][otherLocalDofIdx] -= wijk[faceIdx][0][localDir];
                }
                // the current face is a Dirichlet face and creates entries in D & maybe B
                else
                {
                    D[faceIdx][otherLocalDofIdx] -= wijk[faceIdx][0][localDir];
                    if (!curIsDirichlet)
                        B[curLocalDofIdx][otherLocalDofIdx] += wijk[faceIdx][0][localDir];
                }

                // add entries related to pressures at the scv centers (dofs)
                const auto posScvLocalDofIdx = posLocalScv.localDofIndex();
                D[faceIdx][posScvLocalDofIdx] += wijk[faceIdx][0][localDir];

                if (!curIsDirichlet)
                    B[curLocalDofIdx][posScvLocalDofIdx] -= wijk[faceIdx][0][localDir];
            }

            // If we are on an interior face, add values from negative sub volume
            if (!curGlobalScvf.boundary())
            {
                // loop over all the outside neighbors of this face and add entries
                for (unsigned int idxInOutside = 0; idxInOutside < curGlobalScvf.numOutsideScvs(); ++idxInOutside)
                {
                    const auto idxOnScvf = idxInOutside+1;
                    const auto& negLocalScv = iv.localScv( neighborScvIndices[idxOnScvf] );
                    const auto& negGlobalScv = fvGeometry.scv(negLocalScv.gridScvIndex());
                    const auto negTensor = getT(negGlobalScv);

                    // On surface grids, use outside face for "negative" transmissibility calculation
                    const auto& scvf = dim < dimWorld ? fvGeometry.flipScvf(curGlobalScvf.index(), idxInOutside)
                                                        : curGlobalScvf;
                    wijk[faceIdx][idxOnScvf] = computeOmega(negLocalScv, scvf, negTensor);

                    // flip sign on surface grids (since we used the "outside" normal)
                    if (dim < dimWorld)
                        wijk[faceIdx][idxOnScvf] *= -1.0;

                    // go over the coordinate directions in the negative sub volume
                    for (int localDir = 0; localDir < dim; localDir++)
                    {
                        const auto otherLocalScvfIdx = negLocalScv.localScvfIndex(localDir);
                        const auto& otherLocalScvf = iv.localScvf(otherLocalScvfIdx);
                        const auto otherLocalDofIdx = otherLocalScvf.localDofIndex();

                        // check for zero transmissibilities (skip if inside has been zero already)
                        if (otherLocalDofIdx == curLocalDofIdx && !insideZeroWij)
                            if (abs(wijk[faceIdx][idxOnScvf][localDir]) <= wijZeroThresh)
                                faceMarkers.emplace_back( std::make_pair(curLocalDofIdx, faceIdx) );

                        if (!otherLocalScvf.isDirichlet())
                            A[curLocalDofIdx][otherLocalDofIdx] += wijk[faceIdx][idxOnScvf][localDir];
                        else
                            B[curLocalDofIdx][otherLocalDofIdx] -= wijk[faceIdx][idxOnScvf][localDir];

                        // add entries to matrix B
                        B[curLocalDofIdx][negLocalScv.localDofIndex()] += wijk[faceIdx][idxOnScvf][localDir];
                    }
                }
            }
        }
    }

    return faceMarkers;
}

/*!
 * \brief Solves a previously assembled iv-local system of equations
 *        and stores the resulting transmissibilities in the provided
 *        containers within the interaction volume data handle.
 *
 * \param fvGeometry The bound element finite volume geometry
 * \param handle The data handle in which the matrices are stored
 * \param iv The interaction volume
 */
template< class FVElementGeometry, class DataHandle, class IV >
static void solveLocalSystem(const FVElementGeometry& fvGeometry,
                             DataHandle& handle,
                             IV& iv)
{
    assert(iv.numUnknowns() > 0);

    // T = C*(A^-1)*B + D

    handle.AInverse().invert();
    handle.CAInverse().rightmultiply(handle.AInverse());
    handle.T() += multiplyMatrices(handle.CAInverse(), handle.AInverseB());
    handle.AInverseB().leftmultiply(handle.AInverse());

    // On surface grids, compute the "outside" transmissibilities
    using GridView = typename IV::Traits::GridView;
    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;
    if (dim < dimWorld)
    {
        // bring outside tij container to the right size
        auto& tijOut = handle.tijOutside();
        tijOut.resize(iv.numFaces());
        for (unsigned int fIdx = 0; fIdx < iv.numFaces(); ++fIdx)
        {
            const auto& curGlobalScvf = fvGeometry.scvf(iv.localScvf(fIdx).gridScvfIndex());
            const auto numOutsideFaces = curGlobalScvf.boundary() ? 0 : curGlobalScvf.numOutsideScvs();
            // resize each face entry to the right number of outside faces
            tijOut[fIdx].resize(numOutsideFaces);
            std::for_each(tijOut[fIdx].begin(),
                          tijOut[fIdx].end(),
                          [&](auto& v) { InteractionVolumeAssemblerHelper::resizeVector(v, iv.numKnowns()); });
        }

        // compute outside transmissibilities
        for (const auto& localFaceData : iv.localFaceData())
        {
            // continue only for "outside" faces
            if (!localFaceData.isOutsideFace())
                continue;

            const auto scvIdx = localFaceData.ivLocalInsideScvIndex();
            const auto scvfIdx = localFaceData.ivLocalScvfIndex();
            const auto idxInOut = localFaceData.scvfLocalOutsideScvfIndex();

            const auto& wijk = handle.omegas()[scvfIdx][idxInOut+1];
            auto& tij = tijOut[scvfIdx][idxInOut];

            tij = 0.0;
            for (unsigned int dir = 0; dir < dim; dir++)
            {
                // the scvf corresponding to this local direction in the scv
                const auto& scvf = iv.localScvf(iv.localScv(scvIdx).localScvfIndex(dir));

                // on interior faces the coefficients of the AB matrix come into play
                if (!scvf.isDirichlet())
                {
                    auto tmp = handle.AInverseB()[scvf.localDofIndex()];
                    tmp *= wijk[dir];
                    tij -= tmp;
                }
                else
                    tij[scvf.localDofIndex()] -= wijk[dir];

                // add entry from the scv unknown
                tij[scvIdx] += wijk[dir];
            }
        }
    }
}

} // namespace Detail
#endif // DOXYGEN


/*!
 * \brief Assembles the matrices involved in the flux
 *        expressions and the local system of equations
 *        within an interaction volume in an mpfa-o type way.
 *
 * \param handle The data handle in which the matrices are stored
 * \param iv The interaction volume
 * \param getT Lambda to evaluate the scv-wise tensors
 * \param wijZeroThresh Threshold below which transmissibilities are taken to be zero.
 *                      On the basis of this threshold, trivial (0 = 0) rows in the
 *                      A matrix are identified and modified accordingly in order to
 *                      avoid ending up with singular matrices. This can occur when the
 *                      tensor is zero in some cells.
 */
template<class DataHandle, class IV, class TensorFunc, class FVElementGeometry, class Scalar = double>
void assembleMatrices(DataHandle& handle,
                      IV& iv,
                      const TensorFunc& getT,
                      const FVElementGeometry& fvGeometry,
                      Scalar wijZeroThresh = 0.0)
{
    const auto zeroRows = Detail::assembleLocalMatrices(
        fvGeometry,
        handle.AInverse(),
        handle.AInverseB(),
        handle.CAInverse(),
        handle.T(),
        handle.omegas(),
        iv,
        getT,
        wijZeroThresh
    );

    // maybe solve the local system
    if (iv.numUnknowns() > 0)
    {
        // for now, B is carried in what will be AB after solve
        auto& B = handle.AInverseB();
        auto& A = handle.AInverse();

        // If the tensor is zero in some cells, we might have zero rows in the matrix.
        // In this case we set the diagonal entry of A to 1.0 and the row of B to zero.
        // This will ultimatively lead to zero transmissibilities for this face.
        for (const auto& zeroRowIndices : zeroRows)
        {
            const auto zeroRowDofIdx = zeroRowIndices.first;
            for (auto& row : A)
                row[zeroRowDofIdx] = 0.0;
            A[zeroRowDofIdx] = 0.0;
            A[zeroRowDofIdx][zeroRowDofIdx] = 1.0;
            B[zeroRowDofIdx] = 0.0;
        }

        Detail::solveLocalSystem(fvGeometry, handle, iv);

        // make sure the inverse of A now carries zeroes in the zero rows, as
        // well as CA and T. AB will have the zeroes already because we set the rows
        // of B to zero above.
        for (const auto& zeroRowIndices : zeroRows)
        {
            const auto faceIdx = zeroRowIndices.second;
            A[zeroRowIndices.first] = 0.0;
            handle.CAInverse()[faceIdx] = 0.0;
            handle.T()[faceIdx] = 0.0;

            // reset outside transmissibilities on surface grids
            static constexpr int dim = IV::Traits::GridView::dimension;
            static constexpr int dimWorld = IV::Traits::GridView::dimensionworld;
            if constexpr (dim < dimWorld)
                std::for_each( handle.tijOutside()[faceIdx].begin(),
                                handle.tijOutside()[faceIdx].end(),
                                [] (auto& outsideTij) { outsideTij = 0.0; } );
        }
    }
}

/*!
 * \brief Assembles the force contributions per face, involved
 *        in the flux expressions and the local system of equations
 *        within an interaction volume in an mpfa-o type way.
 *
 * \param handle The data handle in which the matrices are stored
 * \param iv The interaction volume
 * \param getT Lambda to evaluate the scv-wise tensors
 * \param getF Lambda to evaluate the scv-wise forces
 */
template<class DataHandle, class IV, class TensorFunc, class ForceFunc, class FVElementGeometry>
void assembleForces(DataHandle& handle,
                    IV& iv,
                    const TensorFunc& getT,
                    const ForceFunc& getF,
                    const FVElementGeometry& fvGeometry)
{
    using GridView = typename IV::Traits::GridView;
    using Extrusion = Extrusion_t<typename FVElementGeometry::GridGeometry>;

    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;
    static constexpr bool isSurfaceGrid = dim < dimWorld;

    InteractionVolumeAssemblerHelper::resizeVector(handle.forces(), iv.numFaces());
    InteractionVolumeAssemblerHelper::resizeVector(handle.deltaForces(), iv.numUnknowns());
    if constexpr (isSurfaceGrid)
        InteractionVolumeAssemblerHelper::resizeVector(handle.outsideForces(), iv.numFaces());

    //! For each face, we...
    //! - compute the term \f$ \alpha := \mathbf{n}^T \mathbf{T} \mathbf{f} \f$ in each neighboring cell
    //! - compute \f$ \delta \alpha^ = \sum{\alpha_{outside, i}} - \alpha_{inside} \f$
    using NodalIndexSet = CCMpfaDualGridNodalIndexSet<GridView>;
    using LocalIndexType = typename NodalIndexSet::LocalIndexType;


    for (LocalIndexType faceIdx = 0; faceIdx < iv.numFaces(); ++faceIdx)
    {
        const auto& localScvf = iv.localScvf(faceIdx);
        const auto& gridScvf = fvGeometry.scvf(localScvf.gridScvfIndex());
        const auto insideLocalScvIdx = localScvf.neighboringLocalScvIndices()[0];
        const auto& insideScv = fvGeometry.scv(iv.localScv(insideLocalScvIdx).gridScvIndex());
        std::visit([&] (const auto& tensor) {
            handle.forces()[faceIdx] = 1.0/*todo:extrusion*/
                                    *vtmv(gridScvf.unitOuterNormal(), tensor, getF(insideScv))
                                    *Extrusion::area(fvGeometry, gridScvf);
        }, getT(insideScv));

        const auto localDofIdx = localScvf.localDofIndex();
        if (!localScvf.isDirichlet())
            handle.deltaForces()[localDofIdx] = 0.0;
        if (gridScvf.boundary())
            continue;

        const auto numOutsideFaces = gridScvf.numOutsideScvs();
        handle.deltaForces()[localDofIdx] = handle.forces()[faceIdx];

        if constexpr (isSurfaceGrid)
        {
            InteractionVolumeAssemblerHelper::resizeVector(handle.outsideForces()[faceIdx], numOutsideFaces);
            std::fill(handle.outsideForces().begin(), handle.outsideForces().end(), 0.0);
        }

        for (unsigned int idxInOutside = 0; idxInOutside < numOutsideFaces; ++idxInOutside)
        {
            const auto outLocalScvIdx = localScvf.neighboringLocalScvIndices()[idxInOutside+1];
            const auto& outGlobalScv = fvGeometry.scv(iv.localScv(outLocalScvIdx).gridScvIndex());

            // todo:extrusion
            std::visit([&] (const auto& outTensor) {
                if constexpr (isSurfaceGrid)
                {
                    const auto& flipScvf = fvGeometry.flipScvf(gridScvf.index(), idxInOutside);
                    handle.outsideForces()[faceIdx][idxInOutside] = 1.0*vtmv(
                        flipScvf.unitOuterNormal(), outTensor, getF(outGlobalScv)
                    )*Extrusion::area(fvGeometry, flipScvf);
                    handle.deltaForces()[localDofIdx] -= handle.outsideForces()[faceIdx][idxInOutside];
                }
                else
                    handle.deltaForces()[localDofIdx] -= -1.0
                                                         *vtmv(gridScvf.unitOuterNormal(), outTensor, getF(outGlobalScv))
                                                         *Extrusion::area(fvGeometry, gridScvf);
            }, getT(outGlobalScv));
        }
    }
}

} // end namespace Dumux

#endif
