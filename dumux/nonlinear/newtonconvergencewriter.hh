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
 * \ingroup Nonlinear
 * \brief This class provides the infrastructure to write the
 *        convergence behaviour of the newton method into a VTK file.
 */
#ifndef DUMUX_NEWTON_CONVERGENCE_WRITER_HH
#define DUMUX_NEWTON_CONVERGENCE_WRITER_HH

#include <string>

#include <dune/grid/io/file/vtk/vtksequencewriter.hh>
#include <dumux/discretization/method.hh>

namespace Dumux {

//! provide an interface as a form of type erasure
//! this is the minimal requirements a convergence write passed to a newton method has to fulfill
template <class SolutionVector>
struct ConvergenceWriterInterface
{
    virtual ~ConvergenceWriterInterface() = default;

    virtual void write(const SolutionVector &uLastIter, const SolutionVector &deltaU, const SolutionVector &residual) {}
};

/*!
 * \ingroup Nonlinear
 * \brief Writes the intermediate solutions for every Newton iteration
 * \note To use this create a shared_ptr to an instance of this class in the main file
 *       and pass it to newton.solve(x, convergencewriter). You can use the reset method
 *       to write out multiple Newton solves with a unique id, if you don't call use all
 *       Newton iterations just come after each other in the pvd file.
 */
template <class GridGeometry, class SolutionVector>
class NewtonConvergenceWriter : public ConvergenceWriterInterface<SolutionVector>
{
    using GridView = typename GridGeometry::GridView;
    static constexpr auto numEq = SolutionVector::block_type::dimension;
    using Scalar = typename SolutionVector::block_type::value_type;

    static_assert(GridGeometry::discMethod != DiscretizationMethod::staggered,
                  "This convergence writer does not work for the staggered method, use the StaggeredNewtonConvergenceWriter instead");
public:
    /*!
     * \brief Constructor
     * \param gridGeometry The finite-volume grid geometry
     * \param name Base name of the vtk output
     */
    NewtonConvergenceWriter(const GridGeometry& gridGeometry,
                            const std::string& name = "newton_convergence")
    : gridGeometry_(gridGeometry)
    , writer_(gridGeometry.gridView(), name, "", "")
    {
        resize();

        if (GridGeometry::discMethod == DiscretizationMethod::box)
        {
            for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
            {
                writer_.addVertexData(x_[eqIdx], "x_" + std::to_string(eqIdx));
                writer_.addVertexData(delta_[eqIdx], "delta_" + std::to_string(eqIdx));
                writer_.addVertexData(def_[eqIdx], "defect_" + std::to_string(eqIdx));
            }
        }
        else
        {
            for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
            {
                writer_.addCellData(x_[eqIdx], "x_" + std::to_string(eqIdx));
                writer_.addCellData(delta_[eqIdx], "delta_" + std::to_string(eqIdx));
                writer_.addCellData(def_[eqIdx], "defect_" + std::to_string(eqIdx));
            }
        }
    }

    //! Resizes the output fields. This has to be called whenever the grid changes
    void resize()
    {
        const auto numDofs = gridGeometry_.numDofs();

        // resize the output fields
        for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
        {
            def_[eqIdx].resize(numDofs);
            delta_[eqIdx].resize(numDofs);
            x_[eqIdx].resize(numDofs);
        }
    }

    //! Reset the convergence writer for a possible next Newton step
    //! You may set a different id in case you don't want the output to be overwritten by the next step
    void reset(std::size_t newId = 0UL)
    { id_ = newId; iteration_ = 0UL; }

    void write(const SolutionVector& uLastIter,
               const SolutionVector& deltaU,
               const SolutionVector& residual) override
    {
        assert(uLastIter.size() == deltaU.size() && uLastIter.size() == residual.size());

        for (std::size_t dofIdxGlobal = 0; dofIdxGlobal < deltaU.size(); ++dofIdxGlobal)
        {
            for (int eqIdx = 0; eqIdx < numEq; ++eqIdx)
            {
                x_[eqIdx][dofIdxGlobal] = uLastIter[dofIdxGlobal][eqIdx];
                delta_[eqIdx][dofIdxGlobal] = - deltaU[dofIdxGlobal][eqIdx];
                def_[eqIdx][dofIdxGlobal] = residual[dofIdxGlobal][eqIdx];
            }
        }

        writer_.write(static_cast<double>(id_) + static_cast<double>(iteration_)/1000);
        ++iteration_;
    }

private:
    std::size_t id_ = 0UL;
    std::size_t iteration_ = 0UL;

    const GridGeometry& gridGeometry_;

    Dune::VTKSequenceWriter<GridView> writer_;

    std::array<std::vector<Scalar>, numEq> def_;
    std::array<std::vector<Scalar>, numEq> delta_;
    std::array<std::vector<Scalar>, numEq> x_;
};

} // end namespace Dumux

#endif
