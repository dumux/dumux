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
 *        convergence behaviour of the newton method for the
 *        staggered discretization scheme into a VTK file.
 */
#ifndef DUMUX_STAGGERED_NEWTON_CONVERGENCE_WRITER_HH
#define DUMUX_STAGGERED_NEWTON_CONVERGENCE_WRITER_HH

#include <string>

#include <dune/grid/io/file/vtk/vtksequencewriter.hh>
#include <dumux/io/vtksequencewriter.hh>
#include <dumux/io/pointcloudvtkwriter.hh>
#include "newtonconvergencewriter.hh"

namespace Dumux {

/*!
 * \ingroup Nonlinear
 * \brief Writes the intermediate solutions for every Newton iteration (for staggered grid scheme)
 * \note To use this create a shared_ptr to an instance of this class in the main file
 *       and pass it to newton.solve(x, convergencewriter). You can use the reset method
 *       to write out multiple Newton solves with a unique id, if you don't call use all
 *       Newton iterations just come after each other in the pvd file.
 */
template <class GridGeometry, class SolutionVector>
class StaggeredNewtonConvergenceWriter : public ConvergenceWriterInterface<SolutionVector>
{
    using GridView = typename GridGeometry::GridView;

    using CellCenterSolutionVector = typename std::decay_t<decltype(std::declval<SolutionVector>()[GridGeometry::cellCenterIdx()])>;
    using FaceSolutionVector = typename std::decay_t<decltype(std::declval<SolutionVector>()[GridGeometry::faceIdx()])>;

    using Scalar = typename CellCenterSolutionVector::block_type::value_type;

    static constexpr auto numEqCellCenter = CellCenterSolutionVector::block_type::dimension;
    static constexpr auto numEqFace = FaceSolutionVector::block_type::dimension;

    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    static_assert(GridGeometry::discMethod == DiscretizationMethod::staggered,
                  "This convergence writer does only work for the staggered method, use the NewtonConvergenceWriter instead");
public:
    /*!
     * \brief Constructor
     * \param gridGeometry The finite volume geometry on the grid view
     * \param name Base name of the vtk output
     */
    StaggeredNewtonConvergenceWriter(const GridGeometry& gridGeometry,
                                     const std::string& name = "newton_convergence")
    : gridGeometry_(gridGeometry)
    , ccWriter_(gridGeometry.gridView(), name, "", "")
    , faceWriter_(std::make_shared<PointCloudVtkWriter<Scalar, GlobalPosition>>(coordinates_))
    , faceSequenceWriter_(faceWriter_, name + "-face", "","",
                          gridGeometry.gridView().comm().rank(),
                          gridGeometry.gridView().comm().size())
    {
        resize();

        for (int eqIdx = 0; eqIdx < numEqCellCenter; ++eqIdx)
        {
            ccWriter_.addCellData(xCellCenter_[eqIdx], "x_" + std::to_string(eqIdx));
            ccWriter_.addCellData(deltaCellCenter_[eqIdx], "delta_" + std::to_string(eqIdx));
            ccWriter_.addCellData(defCellCenter_[eqIdx], "defect_" + std::to_string(eqIdx));
        }
    }

    //! Resizes the output fields. This has to be called whenever the grid changes
    void resize()
    {
        const auto numCellCenterDofs = gridGeometry_.numCellCenterDofs();
        const auto numFaceDofs = gridGeometry_.numFaceDofs();

        // resize the cell center output fields
        for (int eqIdx = 0; eqIdx < numEqCellCenter; ++eqIdx)
        {
            defCellCenter_[eqIdx].resize(numCellCenterDofs);
            deltaCellCenter_[eqIdx].resize(numCellCenterDofs);
            xCellCenter_[eqIdx].resize(numCellCenterDofs);
        }

        // resize the face output fields
        for (int eqIdx = 0; eqIdx < numEqFace; ++eqIdx)
        {
            defFace_[eqIdx].resize(numFaceDofs);
            deltaFace_[eqIdx].resize(numFaceDofs);
            xFace_[eqIdx].resize(numFaceDofs);
        }

        coordinates_.resize(numFaceDofs);
        for (auto&& facet : facets(gridGeometry_.gridView()))
        {
            const auto dofIdxGlobal = gridGeometry_.gridView().indexSet().index(facet);
            coordinates_[dofIdxGlobal] = facet.geometry().center();
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

        for (std::size_t dofIdxGlobal = 0; dofIdxGlobal < deltaU[GridGeometry::cellCenterIdx()].size(); ++dofIdxGlobal)
        {
            for (int eqIdx = 0; eqIdx < numEqCellCenter; ++eqIdx)
            {
                xCellCenter_[eqIdx][dofIdxGlobal] = uLastIter[GridGeometry::cellCenterIdx()][dofIdxGlobal][eqIdx];
                deltaCellCenter_[eqIdx][dofIdxGlobal] = - deltaU[GridGeometry::cellCenterIdx()][dofIdxGlobal][eqIdx];
                defCellCenter_[eqIdx][dofIdxGlobal] = residual[GridGeometry::cellCenterIdx()][dofIdxGlobal][eqIdx];
            }
        }

        for (int eqIdx = 0; eqIdx < numEqFace; ++eqIdx)
        {
            faceWriter_->addPointData(xFace_[eqIdx], "x_" + std::to_string(eqIdx));
            faceWriter_->addPointData(deltaFace_[eqIdx], "delta_" + std::to_string(eqIdx));
            faceWriter_->addPointData(defFace_[eqIdx], "defect_" + std::to_string(eqIdx));
        }

        for (std::size_t dofIdxGlobal = 0; dofIdxGlobal < deltaU[GridGeometry::faceIdx()].size(); ++dofIdxGlobal)
        {
            for (int eqIdx = 0; eqIdx < numEqFace; ++eqIdx)
            {
                xFace_[eqIdx][dofIdxGlobal] = uLastIter[GridGeometry::faceIdx()][dofIdxGlobal][eqIdx];
                deltaFace_[eqIdx][dofIdxGlobal] = - deltaU[GridGeometry::faceIdx()][dofIdxGlobal][eqIdx];
                defFace_[eqIdx][dofIdxGlobal] = residual[GridGeometry::faceIdx()][dofIdxGlobal][eqIdx];
            }
        }

        ccWriter_.write(static_cast<double>(id_) + static_cast<double>(iteration_)/1000);
        faceSequenceWriter_.write(static_cast<double>(id_) + static_cast<double>(iteration_)/1000);
        ++iteration_;
    }

private:
    std::size_t id_ = 0UL;
    std::size_t iteration_ = 0UL;

    const GridGeometry& gridGeometry_;

    Dune::VTKSequenceWriter<GridView> ccWriter_;

    std::vector<GlobalPosition> coordinates_;
    std::shared_ptr<PointCloudVtkWriter<Scalar, GlobalPosition>> faceWriter_;
    VTKSequenceWriter<PointCloudVtkWriter<Scalar, GlobalPosition>> faceSequenceWriter_;

    std::array<std::vector<Scalar>, numEqCellCenter> defCellCenter_;
    std::array<std::vector<Scalar>, numEqCellCenter> deltaCellCenter_;
    std::array<std::vector<Scalar>, numEqCellCenter> xCellCenter_;

    std::array<std::vector<Scalar>, numEqFace> defFace_;
    std::array<std::vector<Scalar>, numEqFace> deltaFace_;
    std::array<std::vector<Scalar>, numEqFace> xFace_;
};

} // end namespace Dumux

#endif
