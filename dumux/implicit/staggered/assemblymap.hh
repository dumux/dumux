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
 * \brief Stores the face indices corresponding to the neighbors of an element
 *        that contribute to the derivative calculation. This is used for
 *        finite-volume schemes with symmetric sparsity pattern in the global matrix.
 */
#ifndef DUMUX_STAGGERED_ASSEMBLY_MAP_HH
#define DUMUX_STAGGERED_ASSEMBLY_MAP_HH

#include <dune/istl/bcrsmatrix.hh>

#include <dumux/implicit/properties.hh>

namespace Dumux
{

template<class TypeTag>
class StaggeredAssemblyMap
{
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using FluxVariables = typename GET_PROP_TYPE(TypeTag, FluxVariables);

    using IndexType = typename GridView::IndexSet::IndexType;

    using DofTypeIndices = typename GET_PROP(TypeTag, DofTypeIndices);
    typename DofTypeIndices::CellCenterIdx cellCenterIdx;
    typename DofTypeIndices::FaceIdx faceIdx;

    using CellCenterToCellCenterMap = std::vector<std::vector<IndexType>>;
    using CellCenterToFaceMap = std::vector<std::vector<IndexType>>;
    using FaceToCellCenterMap = std::vector<std::vector<IndexType>>;
    using FaceToFaceMap = std::vector<std::vector<IndexType>>;

public:

    /*!
     * \brief Initialize the AssemblyMap object.
     *
     * \param problem The problem which we want to simulate.
     */
    void init(const Problem& problem)
    {

        const auto numDofsCC = problem.model().numCellCenterDofs();
        const auto numDofsFace = problem.model().numFaceDofs();
        cellCenterToCellCenterMap_.resize(numDofsCC);
        cellCenterToFaceMap_.resize(numDofsCC);
        faceToCellCenterMap_.resize(numDofsFace);
        faceToFaceMap_.resize(numDofsFace);

        cellCenterToCellCenterMap_[0] = std::vector<IndexType>{1,2,3};
        cellCenterToFaceMap_[0] = std::vector<IndexType>{4,5,6};
        faceToCellCenterMap_[0] = std::vector<IndexType>{7,8,9};
        faceToFaceMap_[0] = std::vector<IndexType>{10,11,12};
    }

    template< std::size_t index1, std::size_t index2 >
    const std::vector<IndexType>& operator() (const std::integral_constant< std::size_t, index1 > indexVariable1,
                                              const std::integral_constant< std::size_t, index2 > indexVariable2,
                                              const IndexType globalI)
    {
        DUNE_UNUSED_PARAMETER(indexVariable1);
        DUNE_UNUSED_PARAMETER(indexVariable2);

        if(index1 == cellCenterIdx && index2 == cellCenterIdx)
            return cellCenterToCellCenterMap_[globalI];
        else if(index1 == cellCenterIdx && index2 == faceIdx)
            return cellCenterToFaceMap_[globalI];
        else if(index1 == faceIdx && index2 == cellCenterIdx)
            return cellCenterToFaceMap_[globalI];
        else if(index1 == faceIdx && index2 == faceIdx)
            return faceToFaceMap_[globalI];
        else
            DUNE_THROW(Dune::InvalidStateException, "Invalid indices");
    }



private:
    CellCenterToCellCenterMap cellCenterToCellCenterMap_;
    CellCenterToFaceMap cellCenterToFaceMap_;
    FaceToCellCenterMap faceToCellCenterMap_;
    FaceToFaceMap faceToFaceMap_;
};

} // end namespace Dumux

#endif
