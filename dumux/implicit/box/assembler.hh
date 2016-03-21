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
 * \brief An assembler for the global Jacobian matrix for fully implicit models.
 */
#ifndef DUMUX_BOX_ASSEMBLER_HH
#define DUMUX_BOX_ASSEMBLER_HH

#include <dumux/implicit/properties.hh>
#include <dumux/implicit/assembler.hh>

namespace Dumux {

/*!
 * \ingroup ImplicitModel
 * \brief An assembler for the global Jacobian matrix for fully implicit models.
 */
template<class TypeTag>
class BoxAssembler : public ImplicitAssembler<TypeTag>
{
    typedef ImplicitAssembler<TypeTag> ParentType;
    friend class ImplicitAssembler<TypeTag>;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GridView::IndexSet::IndexType IndexType;

    // assemble an interior element
    void assembleElement_(const Element &element)
    {

        this->model_().localJacobian().assemble(element, this->matrix());

        for(int vIdxLocal = 0; vIdxLocal < element.subEntities(dim); ++vIdxLocal)
        {
            auto globalI = vertexMapper_().subIndex(element, dim, vIdxLocal);

            // update the right hand side
            this->residual_[globalI] += this->model_().localJacobian().residual(vIdxLocal);

            for (int j = 0; j < residual_[globalI].dimension; ++j)
                assert(std::isfinite(residual_[globalI][j]));
        }
    }

    void setRowSizes_()
    {
        for (const auto& vertex : vertices(this->gridView_()))
        {
            // the global index of the element at hand
            const auto globalI = this->vertexMapper_().index(vertex);
            const auto& stencil = this->model_().stencils(vertex).vertexStencil();

            this->matrix().setrowsize(globalI, stencil.size());
        }
        this->matrix().endrowsizes();

    }

    void addIndices_()
    {
        for (const auto& vertex : vertices(this->gridView_()))
        {
            // the global index of the element at hand
            const auto globalI = this->vertexMapper_().index(vertex);
            const auto& stencil = this->model_().stencils(vertex).vertexStencil();

            for (auto&& globalJ : stencil)
                this->matrix().addindex(globalI, globalJ);

            this->matrix().endindices();
        }
    }
};

} // namespace Dumux

#endif
