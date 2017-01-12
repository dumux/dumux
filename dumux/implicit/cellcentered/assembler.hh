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
#ifndef DUMUX_CC_ASSEMBLER_HH
#define DUMUX_CC_ASSEMBLER_HH

#include <dune/istl/matrixindexset.hh>

#include <dumux/implicit/properties.hh>
#include <dumux/implicit/assembler.hh>

namespace Dumux {

/*!
 * \ingroup ImplicitModel
 * \brief An assembler for the global Jacobian matrix for fully implicit models.
 */
template<class TypeTag>
class CCAssembler : public ImplicitAssembler<TypeTag>
{
    friend class ImplicitAssembler<TypeTag>;
    using JacobianMatrix = typename GET_PROP_TYPE(TypeTag, JacobianMatrix);

    void createMatrix_()
    {
        const auto numDofs = this->problem_().model().numDofs();

        // allocate raw matrix
        this->matrix_ = std::make_shared<JacobianMatrix>(numDofs, numDofs, JacobianMatrix::random);

        // get occupation pattern of the matrix
        Dune::MatrixIndexSet occupationPattern;
        occupationPattern.resize(numDofs, numDofs);

        const auto& assemblyMap = this->problem_().model().localJacobian().assemblyMap();
        for (unsigned int globalI = 0; globalI < numDofs; ++globalI)
        {
            occupationPattern.add(globalI, globalI);
            for (const auto& dataJ : assemblyMap[globalI])
                occupationPattern.add(dataJ.globalJ, globalI);
        }

        // export pattern to matrix
        occupationPattern.exportIdx(*this->matrix_);
    }
};

} // namespace Dumux

#endif
