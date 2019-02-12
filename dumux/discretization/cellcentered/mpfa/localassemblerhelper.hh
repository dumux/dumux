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
 * \brief A class that contains helper functions as well as functionality
 *        which is common to different mpfa schemes and which solely
 *        operate on the interaction volume.
 */
#ifndef DUMUX_DISCRETIZATION_CC_MPFA_LOCAL_ASSEMBLER_HELPER_HH
#define DUMUX_DISCRETIZATION_CC_MPFA_LOCAL_ASSEMBLER_HELPER_HH

#include <algorithm>
#include <vector>
#include <cassert>
#include <utility>
#include <type_traits>

#include <dumux/common/typetraits/isvalid.hh>

namespace Dumux {

/*!
 * \ingroup CCMpfaDiscretization
 * \brief A class that contains helper functions as well as functionality
 *        which is common to different mpfa schemes and which solely
 *        operate on the interaction volume.
 */
class InteractionVolumeAssemblerHelper
{
    // Helper structs to detect if matrix has resize function
    struct HasMatrixResize
    {
        template<class M>
        auto operator()(const M& m) -> decltype(std::declval<M>().resize(0, 0))
        {}
    };

    // Helper structs to detect if vector has resize function
    struct HasVectorResize
    {
        template<class V>
        auto operator()(const V& v) -> decltype(std::declval<V>().resize(0))
        {}
    };

    template<class Matrix>
    static constexpr auto matrixHasResizeFunction()
    { return decltype( isValid(HasMatrixResize())(std::declval<Matrix>()) )::value; }

    template<class Vector>
    static constexpr auto vectorHasResizeFunction()
    { return decltype( isValid(HasVectorResize())(std::declval<Vector>()) )::value; }

public:
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
        handle.A().invert();
        handle.CA().rightmultiply(handle.A());
        handle.T() += multiplyMatrices(handle.CA(), handle.AB());
        handle.AB().leftmultiply(handle.A());

        // On surface grids, compute the "outside" transmissibilities
        using GridView = typename IV::Traits::GridView;
        static constexpr int dim = GridView::dimension;
        static constexpr int dimWorld = GridView::dimensionworld;
        if (dim < dimWorld)
        {
            // bring outside tij container to the right size
            auto& tijOut = handle.tijOutside();
            tijOut.resize(iv.numFaces());
            using LocalIndexType = typename IV::Traits::IndexSet::LocalIndexType;
            for (LocalIndexType fIdx = 0; fIdx < iv.numFaces(); ++fIdx)
            {
                const auto& curGlobalScvf = fvGeometry.scvf(iv.localScvf(fIdx).gridScvfIndex());
                const auto numOutsideFaces = curGlobalScvf.boundary() ? 0 : curGlobalScvf.numOutsideScvs();
                // resize each face entry to the right number of outside faces
                tijOut[fIdx].resize(numOutsideFaces);
                std::for_each(tijOut[fIdx].begin(),
                              tijOut[fIdx].end(),
                              [&](auto& v) { resizeVector(v, iv.numKnowns()); });
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
                        auto tmp = handle.AB()[scvf.localDofIndex()];
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

    /*!
     * \brief Assembles the vector of face unknowns within an interaction volume.
     * \note  This requires the data handle to be fully assembled already.
     *
     * \param handle The data handle in which the vector is stored
     * \param iv The interaction volume
     */
    template< class DataHandle, class IV >
    static typename IV::Traits::MatVecTraits::FaceVector
    assembleFaceUnkowns(const DataHandle& handle, const IV& iv)
    {
        typename IV::Traits::MatVecTraits::FaceVector u;
        resizeVector(u, iv.numFaces());

        handle.AB().mv(handle.uj(), u);

        // maybe add gravity terms
        if (handle.deltaG().size() == iv.numUnknowns())
            handle.A().umv(handle.deltaG(), u);

        return u;
    }

    /*!
     * \brief Assembles the solution gradients in the
     *        sub-control volumes within an interaction volume.
     * \note  This requires the data handle to be fully assembled already.
     *
     * \param handle The data handle in which the vector is stored
     * \param iv The interaction volume
     */
    template< class DataHandle, class IV >
    static std::vector< typename IV::Traits::LocalScvType::GlobalCoordinate >
    assembleScvGradients(const DataHandle& handle, const IV& iv)
    {
        const auto u = assembleFaceUnkowns(handle, iv);

        using LocalScv = typename IV::Traits::LocalScvType;
        using Gradient = typename LocalScv::GlobalCoordinate;

        std::vector<Gradient> result; result.reserve(iv.numScvs());
        for (unsigned int scvIdx = 0; scvIdx < iv.numScvs(); ++scvIdx)
        {
            const auto& scv = iv.localScv(scvIdx);

            Gradient gradU(0.0);
            for (unsigned int dir = 0; dir < LocalScv::myDimension; ++dir)
            {
                auto nu = scv.nu(dir);

                // obtain face pressure
                const auto& scvf = iv.localScvf( scv.localScvfIndex(dir) );
                const auto faceU = !scvf.isDirichlet() ? u[scvf.localDofIndex()]
                                                       : handle.uj()[scvf.localDofIndex()];

                nu *= faceU - handle.uj()[scv.localDofIndex()];
                gradU += nu;
            }

            gradU /= scv.detX();
            result.emplace_back( std::move(gradU) );
        }

        return result;
    }

    //! resizes a matrix to the given sizes (specialization for dynamic matrix type)
    template< class Matrix,
              class size_type,
              std::enable_if_t<matrixHasResizeFunction<Matrix>(), int> = 0 >
    static void resizeMatrix(Matrix& M, size_type rows, size_type cols)
    {
        M.resize(rows, cols);
    }

    //! resizes a matrix to the given sizes (specialization for static matrix type - do nothing)
    template< class Matrix,
              class size_type,
              std::enable_if_t<!matrixHasResizeFunction<Matrix>(), int> = 0 >
    static void resizeMatrix(Matrix& M, size_type rows, size_type cols)
    {}

    //! resizes a vector to the given size (specialization for dynamic matrix type)
    template< class Vector,
              class size_type,
              std::enable_if_t<vectorHasResizeFunction<Vector>(), int> = 0 >
    static void resizeVector(Vector& v, size_type size)
    {
        v.resize(size);
    }

    //! resizes a vector to the given size (specialization for static vector type - do nothing)
    template< class Vector,
              class size_type,
              std::enable_if_t<!vectorHasResizeFunction<Vector>(), int> = 0 >
    static void resizeVector(Vector& v, size_type rows)
    {}
};

} // end namespace Dumux

#endif
