// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
 * \ingroup CCMpfaDiscretization
 * \brief Class for the assembly of the local systems of equations
 *        involved in the transmissibility computaion in an mpfa-o type manner.
 */
#ifndef DUMUX_DISCRETIZATION_CC_MPFA_O_LOCAL_ASSEMBLER_HH
#define DUMUX_DISCRETIZATION_CC_MPFA_O_LOCAL_ASSEMBLER_HH

#include <algorithm>

#include <dumux/discretization/cellcentered/mpfa/localassemblerbase.hh>
#include <dumux/discretization/cellcentered/mpfa/localassemblerhelper.hh>
#include <dumux/discretization/cellcentered/mpfa/computetransmissibility.hh>

namespace Dumux {

/*!
 * \ingroup CCMpfaDiscretization
 * \brief Specialization of the interaction volume-local
 *        assembler class for the schemes using an mpfa-o type assembly.
 *
 * \tparam P The problem type
 * \tparam EG The element finite volume geometry
 * \tparam EV The element volume variables type
 */
template< class P, class EG, class EV >
class MpfaOInteractionVolumeAssembler
: public InteractionVolumeAssemblerBase< P, EG, EV >
{
    using ParentType = InteractionVolumeAssemblerBase< P, EG, EV >;
    using Helper = InteractionVolumeAssemblerHelper;

    template< class IV >
    using Scalar = typename IV::Traits::MatVecTraits::FaceVector::value_type;

public:
    //! Pull up constructor of the base class
    using ParentType::ParentType;

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
    template< class DataHandle, class IV, class TensorFunc >
    void assembleMatrices(DataHandle& handle, IV& iv, const TensorFunc& getT, Scalar<IV> wijZeroThresh = 0.0)
    {
        const auto zeroRows = assembleLocalMatrices_(handle.A(), handle.AB(), handle.CA(), handle.T(), handle.omegas(), iv, getT, wijZeroThresh);

        // maybe solve the local system
        if (iv.numUnknowns() > 0)
        {
            // for now, B is carried in what will be AB after solve
            auto& B = handle.AB();
            auto& A = handle.A();

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

            Helper::solveLocalSystem(this->fvGeometry(), handle, iv);

            // make sure the inverse of A now carries zeroes in the zero rows, as
            // well as CA and T. AB will have the zeroes already because we set the rows
            // of B to zero above.
            for (const auto& zeroRowIndices : zeroRows)
            {
                const auto faceIdx = zeroRowIndices.second;
                A[zeroRowIndices.first] = 0.0;
                handle.CA()[faceIdx] = 0.0;
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
     * \brief Assembles the vector of primary (cell) unknowns and (maybe)
     *        Dirichlet boundary conditions within an interaction volume.
     *
     * \param handle The data handle in which the vector is stored
     * \param iv The mpfa-o interaction volume
     * \param getU Lambda to obtain the desired cell/Dirichlet value from vol vars
     */
    template< class DataHandle, class IV, class GetU >
    void assembleU(DataHandle& handle, const IV& iv, const GetU& getU)
    {
        auto& u = handle.uj();
        Helper::resizeVector(u, iv.numKnowns());

        // put the cell unknowns first, then Dirichlet values
        typename IV::Traits::IndexSet::LocalIndexType i = 0;
        for (; i < iv.numScvs(); i++)
            u[i] = getU( this->elemVolVars()[iv.localScv(i).gridScvIndex()] );
        for (const auto& data : iv.dirichletData())
            u[i++] = getU( this->elemVolVars()[data.volVarIndex()] );
    }

private:
    /*!
     * \brief Assemble the matrices involved in the flux expressions
     *        across the scvfs inside an interaction volume as well as those involved
     *        in the interaction volume-local system of equations resulting from flux
     *        and solution continuity across the scvfs.
     *
     *        Flux expressions: \f$\mathbf{f} = \mathbf{C} \bar{\mathbf{u}} + \mathbf{D} \mathbf{u}\f$.
     *
     *        Continuity equations: \f$\mathbf{A} \, \bar{\mathbf{u}} = \mathbf{B} \, \mathbf{u}\f$.
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
    template< class IV, class TensorFunc >
    auto assembleLocalMatrices_(typename IV::Traits::MatVecTraits::AMatrix& A,
                                typename IV::Traits::MatVecTraits::BMatrix& B,
                                typename IV::Traits::MatVecTraits::CMatrix& C,
                                typename IV::Traits::MatVecTraits::DMatrix& D,
                                typename IV::Traits::MatVecTraits::OmegaStorage& wijk,
                                IV& iv, const TensorFunc& getT,
                                Scalar<IV> wijZeroThresh)
    {
        using LocalIndexType = typename IV::Traits::IndexSet::LocalIndexType;
        static constexpr int dim = IV::Traits::GridView::dimension;
        static constexpr int dimWorld = IV::Traits::GridView::dimensionworld;

        std::vector< std::pair<LocalIndexType, LocalIndexType> > faceMarkers;
        Helper::resizeVector(wijk, iv.numFaces());

        // if only Dirichlet faces are present in the iv,
        // the matrices A, B & C are undefined and D = T
        if (iv.numUnknowns() == 0)
        {
            // resize & reset D matrix
            Helper::resizeMatrix(D, iv.numFaces(), iv.numKnowns()); D = 0.0;

            // Loop over all the faces, in this case these are all dirichlet boundaries
            for (LocalIndexType faceIdx = 0; faceIdx < iv.numFaces(); ++faceIdx)
            {
                const auto& curLocalScvf = iv.localScvf(faceIdx);
                const auto& curGlobalScvf = this->fvGeometry().scvf(curLocalScvf.gridScvfIndex());
                const auto& neighborScvIndices = curLocalScvf.neighboringLocalScvIndices();

                // get tensor in "positive" sub volume
                const auto& posLocalScv = iv.localScv(neighborScvIndices[0]);
                const auto& posGlobalScv = this->fvGeometry().scv(posLocalScv.gridScvIndex());
                const auto& posVolVars = this->elemVolVars()[posGlobalScv];
                const auto tensor = getT(posVolVars);

                // the omega factors of the "positive" sub volume
                Helper::resizeVector(wijk[faceIdx], /*no outside scvs present*/1);
                wijk[faceIdx][0] = computeMpfaTransmissibility<EG>(posLocalScv, curGlobalScvf, tensor, posVolVars.extrusionFactor());

                const auto posScvLocalDofIdx = posLocalScv.localDofIndex();
                for (LocalIndexType localDir = 0; localDir < dim; localDir++)
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
            Helper::resizeMatrix(A, iv.numUnknowns(), iv.numUnknowns()); A = 0.0;
            Helper::resizeMatrix(B, iv.numUnknowns(), iv.numKnowns());   B = 0.0;
            Helper::resizeMatrix(C, iv.numFaces(), iv.numUnknowns());    C = 0.0;
            Helper::resizeMatrix(D, iv.numFaces(), iv.numKnowns());      D = 0.0;

            for (LocalIndexType faceIdx = 0; faceIdx < iv.numFaces(); ++faceIdx)
            {
                const auto& curLocalScvf = iv.localScvf(faceIdx);
                const auto& curGlobalScvf = this->fvGeometry().scvf(curLocalScvf.gridScvfIndex());
                const auto curIsDirichlet = curLocalScvf.isDirichlet();
                const auto curLocalDofIdx = curLocalScvf.localDofIndex();

                // get tensor in "positive" sub volume
                const auto& neighborScvIndices = curLocalScvf.neighboringLocalScvIndices();
                const auto& posLocalScv = iv.localScv(neighborScvIndices[0]);
                const auto& posGlobalScv = this->fvGeometry().scv(posLocalScv.gridScvIndex());
                const auto& posVolVars = this->elemVolVars()[posGlobalScv];
                const auto tensor = getT(posVolVars);

                // the omega factors of the "positive" sub volume
                Helper::resizeVector(wijk[faceIdx], neighborScvIndices.size());
                wijk[faceIdx][0] = computeMpfaTransmissibility<EG>(posLocalScv, curGlobalScvf, tensor, posVolVars.extrusionFactor());

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
                        const auto& negGlobalScv = this->fvGeometry().scv(negLocalScv.gridScvIndex());
                        const auto& negVolVars = this->elemVolVars()[negGlobalScv];
                        const auto negTensor = getT(negVolVars);

                        // On surface grids, use outside face for "negative" transmissibility calculation
                        const auto& scvf = dim < dimWorld ? this->fvGeometry().flipScvf(curGlobalScvf.index(), idxInOutside)
                                                          : curGlobalScvf;
                        wijk[faceIdx][idxOnScvf] = computeMpfaTransmissibility<EG>(negLocalScv, scvf, negTensor, negVolVars.extrusionFactor());

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
};

} // end namespace

#endif
