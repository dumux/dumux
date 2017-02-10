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
 *
 * \brief This class provides the infrastructure to write the
 *        convergence behaviour of the newton method into a VTK file.
 */
#ifndef DUMUX_STAGGERED_NEWTON_CONVERGENCE_WRITER_HH
#define DUMUX_STAGGERED_NEWTON_CONVERGENCE_WRITER_HH

#include <dune/grid/io/file/vtk/vtksequencewriter.hh>

#include <dumux/common/basicproperties.hh>

#include "newtoncontroller.hh"

namespace Dumux
{

namespace Properties
{
NEW_PROP_TAG(NewtonController);
NEW_PROP_TAG(SolutionVector);
}

/*!
 * \ingroup Newton
 * \brief Writes the intermediate solutions during
 *        the Newton scheme
 */
template <class TypeTag>
class StaggeredNewtonConvergenceWriter
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using NewtonController = typename GET_PROP_TYPE(TypeTag, NewtonController);
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);

    static constexpr int numEq = GET_PROP_VALUE(TypeTag, NumEq);
    static constexpr bool isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox);

public:
    StaggeredNewtonConvergenceWriter(NewtonController &ctl, const GridView& gridView)
    : ctl_(ctl),
      writer_(gridView, "convergence", "", "")
    {
        DUNE_THROW(Dune::NotImplemented, "StaggeredNewtonConvergenceWriter not implemented! Do it yourself!");
//         timeStepIndex_ = 0;
//         iteration_ = 0;
//
//         const auto numDofs = ctl_.method().model().numDofs();
//         for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
//         {
//             def_[eqIdx].resize(numDofs);
//             delta_[eqIdx].resize(numDofs);
//             x_[eqIdx].resize(numDofs);
//         }
//
//         if (isBox)
//         {
//             for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
//             {
//                 writer_.addVertexData(x_[eqIdx], "x_" + std::to_string(eqIdx));
//                 writer_.addVertexData(delta_[eqIdx], "delta_" + std::to_string(eqIdx));
//                 writer_.addVertexData(def_[eqIdx], "defect_" + std::to_string(eqIdx));
//             }
//         }
//         else
//         {
//             for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
//             {
//                 writer_.addCellData(x_[eqIdx], "x_" + std::to_string(eqIdx));
//                 writer_.addCellData(delta_[eqIdx], "delta_" + std::to_string(eqIdx));
//                 writer_.addCellData(def_[eqIdx], "defect_" + std::to_string(eqIdx));
//             }
//         }
    }

    void advanceTimeStep()
    {
        ++timeStepIndex_;
        iteration_ = 0;
    }

    void advanceIteration()
    {
        ++iteration_;
    }

    void write(const SolutionVector &uLastIter,
               const SolutionVector &deltaU)
    {
//         SolutionVector residual(uLastIter);
//         ctl_.method().model().globalResidual(residual, uLastIter);
//
//         for (unsigned int dofIdxGlobal = 0; dofIdxGlobal < deltaU.size(); dofIdxGlobal++)
//         {
// //             for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
// //             {
// //                 x_[eqIdx][dofIdxGlobal] = uLastIter[dofIdxGlobal][eqIdx];
// //                 delta_[eqIdx][dofIdxGlobal] = - deltaU[dofIdxGlobal][eqIdx];
// //                 def_[eqIdx][dofIdxGlobal] = residual[dofIdxGlobal][eqIdx];
// //             }
//         }
//
//         writer_.write(timeStepIndex_ + iteration_ / 100.0);
    }

private:
    int timeStepIndex_;
    int iteration_;

    std::array<std::vector<Scalar>, numEq> def_;
    std::array<std::vector<Scalar>, numEq> delta_;
    std::array<std::vector<Scalar>, numEq> x_;

    NewtonController &ctl_;
    Dune::VTKSequenceWriter<GridView> writer_;
};

} // namespace Dumux

#endif
