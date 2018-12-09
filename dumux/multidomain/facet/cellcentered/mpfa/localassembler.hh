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
 * \ingroup CCMpfaDiscretization
 * \brief Class for the assembly of the local systems of equations
 *        involved in the transmissibility computaion in the mpfa-o
 *        scheme with domain coupling across the element facets.
 */
#ifndef DUMUX_DISCRETIZATION_CC_MPFA_O_FACET_COUPLING_LOCAL_ASSEMBLER_HH
#define DUMUX_DISCRETIZATION_CC_MPFA_O_FACET_COUPLING_LOCAL_ASSEMBLER_HH

#include <dumux/common/math.hh>
#include <dumux/common/matrixvectorhelper.hh>

#include <dumux/discretization/cellcentered/mpfa/methods.hh>
#include <dumux/discretization/cellcentered/mpfa/localassembler.hh>
#include <dumux/discretization/cellcentered/mpfa/computetransmissibility.hh>

namespace Dumux
{

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
class InteractionVolumeAssemblerImpl< P, EG, EV, MpfaMethods::oMethodFacetCoupling >
      : public InteractionVolumeAssemblerBase< P, EG, EV >
{
    using ParentType = InteractionVolumeAssemblerBase< P, EG, EV >;

public:
    //! Use the constructor of the base class
    using ParentType::ParentType;

    /*!
     * \brief Assembles the transmissibility matrix within an
     *        interaction volume of the mpfa-o scheme in the
     *        context of .
     *
     * \tparam IV The interaction volume type implementation
     * \tparam TensorFunc Lambda to obtain the tensor w.r.t.
     *                   which the local system is to be solved
     *
     * \param T The transmissibility matrix to be assembled
     * \param iv The interaction volume
     * \param getT Lambda to evaluate the scv-wise tensors
     */
    template< class IV, class TensorFunc >
    void assemble(typename IV::Traits::MatVecTraits::TMatrix& T, IV& iv, const TensorFunc& getT)
    {
        // assemble D into T directly
        assembleLocalMatrices_(iv.A(), iv.B(), iv.C(), T, iv, getT);

        // maybe solve the local system
        if (iv.numUnknowns() > 0)
        {
            // T = C*A^-1*B + D
            iv.A().invert();
            iv.C().rightmultiply(iv.A());
            T += multiplyMatrices(iv.C(), iv.B());
        }
    }

    /*!
     * \brief Assembles the vector of primary (cell) unknowns and (maybe)
     *        Dirichlet boundary conditions within an interaction volume.
     *
     * \tparam IV The interaction volume type implementation
     * \tparam GetU Lambda to obtain the cell unknowns from grid indices
     *
     * \param u The vector to be filled with the cell unknowns
     * \param iv The mpfa-o interaction volume
     * \param getU Lambda to obtain the desired cell/Dirichlet value from grid index
     */
    template< class IV, class GetU >
    void assemble(typename IV::Traits::MatVecTraits::CellVector& u, const IV& iv, const GetU& getU)
    {
        // put the cell pressures first
        using LocalIndexType = typename IV::Traits::IndexSet::LocalIndexType;
        for (LocalIndexType i = 0; i < iv.numScvs(); ++i)
            u[i] = getU( iv.localScv(i).globalScvIndex() );

        // Dirichlet BCs come afterwards
        LocalIndexType i = iv.numScvs();
        for (const auto& data : iv.dirichletData())
            u[i++] = getU( data.volVarIndex() );
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
     * \tparam TensorFunc Lambda to obtain the tensor w.r.t.
     *                    which the local system is to be solved
     *
     * \param A The A matrix of the iv-local equation system
     * \param B The B matrix of the iv-local equation system
     * \param C The C matrix of the iv-local flux expressions
     * \param D The D matrix of the iv-local flux expressions
     * \param iv The mpfa-o interaction volume
     * \param getT Lambda to evaluate the scv-wise tensors
     */
    template< class IV, class TensorFunc >
    void assembleLocalMatrices_(typename IV::Traits::MatVecTraits::AMatrix& A,
                                typename IV::Traits::MatVecTraits::BMatrix& B,
                                typename IV::Traits::MatVecTraits::CMatrix& C,
                                typename IV::Traits::MatVecTraits::DMatrix& D,
                                IV& iv, const TensorFunc& getT)
    {
        using LocalIndexType = typename IV::Traits::IndexSet::LocalIndexType;
        static constexpr int dim = IV::Traits::GridView::dimension;
        static constexpr int dimWorld = IV::Traits::GridView::dimensionworld;

        // Matrix D is assumed to have the right size already
        assert(D.rows() == iv.numFaces() && D.cols() == iv.numKnowns());

        // if only Dirichlet faces are present in the iv,
        // the matrices A, B & C are undefined and D = T
        if (iv.numUnknowns() == 0)
        {
            // reset matrix beforehand
            D = 0.0;

            // Loop over all the faces, in this case these are all dirichlet boundaries
            for (LocalIndexType faceIdx = 0; faceIdx < iv.numFaces(); ++faceIdx)
            {
                const auto& curLocalScvf = iv.localScvf(faceIdx);
                const auto& curGlobalScvf = this->fvGeometry().scvf(curLocalScvf.globalScvfIndex());
                const auto& neighborScvIndices = curLocalScvf.neighboringLocalScvIndices();

                // get tensor in "positive" sub volume
                const auto& posLocalScv = iv.localScv(neighborScvIndices[0]);
                const auto& posGlobalScv = this->fvGeometry().scv(posLocalScv.globalScvIndex());
                const auto& posVolVars = this->elemVolVars()[posGlobalScv];
                const auto& posElement = iv.element(neighborScvIndices[0]);
                const auto tensor = getT(this->problem(), posElement, posVolVars, this->fvGeometry(), posGlobalScv);

                // the omega factors of the "positive" sub volume
                const auto wijk = computeMpfaTransmissibility(posLocalScv, curGlobalScvf, tensor, posVolVars.extrusionFactor());

                const auto posScvLocalDofIdx = posLocalScv.localDofIndex();
                for (LocalIndexType localDir = 0; localDir < dim; localDir++)
                {
                    const auto& otherLocalScvf = iv.localScvf( posLocalScv.scvfIdxLocal(localDir) );
                    const auto otherLocalDofIdx = otherLocalScvf.localDofIndex();
                    D[faceIdx][otherLocalDofIdx] -= wijk[localDir];
                    D[faceIdx][posScvLocalDofIdx] += wijk[localDir];
                }
            }
        }
        else
        {
            // we require the matrices A,B,C to have the correct size already
            assert(A.rows() == iv.numUnknowns() && A.cols() == iv.numUnknowns());
            assert(B.rows() == iv.numUnknowns() && B.cols() == iv.numKnowns());
            assert(C.rows() == iv.numFaces() && C.cols() == iv.numUnknowns());

            // reset matrices
            A = 0.0;
            B = 0.0;
            C = 0.0;
            D = 0.0;

            auto& wijk = iv.omegas();
            for (LocalIndexType faceIdx = 0; faceIdx < iv.numFaces(); ++faceIdx)
            {
                const auto& curLocalScvf = iv.localScvf(faceIdx);
                const auto& curGlobalScvf = this->fvGeometry().scvf(curLocalScvf.globalScvfIndex());
                const auto curIsDirichlet = curLocalScvf.isDirichlet();
                const auto curLocalDofIdx = curLocalScvf.localDofIndex();

                // get tensor in "positive" sub volume
                const auto& neighborScvIndices = curLocalScvf.neighboringLocalScvIndices();
                const auto& posLocalScv = iv.localScv(neighborScvIndices[0]);
                const auto& posGlobalScv = this->fvGeometry().scv(posLocalScv.globalScvIndex());
                const auto& posVolVars = this->elemVolVars()[posGlobalScv];
                const auto& posElement = iv.element(neighborScvIndices[0]);
                const auto tensor = getT(this->problem(), posElement, posVolVars, this->fvGeometry(), posGlobalScv);

                // the omega factors of the "positive" sub volume
                wijk[faceIdx][0] = computeMpfaTransmissibility(posLocalScv, curGlobalScvf, tensor, posVolVars.extrusionFactor());

                // go over the coordinate directions in the positive sub volume
                for (unsigned int localDir = 0; localDir < dim; localDir++)
                {
                    const auto& otherLocalScvf = iv.localScvf( posLocalScv.scvfIdxLocal(localDir) );
                    const auto otherLocalDofIdx = otherLocalScvf.localDofIndex();

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
                        const auto& negGlobalScv = this->fvGeometry().scv(negLocalScv.globalScvIndex());
                        const auto& negVolVars = this->elemVolVars()[negGlobalScv];
                        const auto& negElement = iv.element( neighborScvIndices[idxOnScvf] );
                        const auto negTensor = getT(this->problem(), negElement, negVolVars, this->fvGeometry(), negGlobalScv);

                        // On surface grids, use outside face for "negative" transmissibility calculation
                        const auto& scvf = dim < dimWorld ? this->fvGeometry().flipScvf(curGlobalScvf.index(), idxInOutside)
                                                          : curGlobalScvf;
                        wijk[faceIdx][idxOnScvf] = computeMpfaTransmissibility(negLocalScv, scvf, negTensor, negVolVars.extrusionFactor());

                        // flip sign on surface grids (since we used the "outside" normal)
                        if (dim < dimWorld)
                            wijk[faceIdx][idxOnScvf] *= -1.0;

                        // go over the coordinate directions in the positive sub volume
                        for (int localDir = 0; localDir < dim; localDir++)
                        {
                            const auto otherLocalScvfIdx = negLocalScv.scvfIdxLocal(localDir);
                            const auto& otherLocalScvf = iv.localScvf(otherLocalScvfIdx);
                            const auto otherLocalDofIdx = otherLocalScvf.localDofIndex();

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
    }
};

} // end namespace

#endif
