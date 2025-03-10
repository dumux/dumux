// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Newton
 * \brief This class provides the infrastructure to write the
 *        convergence behaviour of the newton method into a VTK file.
 */
#ifndef DUMUX_NEWTON_CONVERGENCE_WRITER_HH
#define DUMUX_NEWTON_CONVERGENCE_WRITER_HH

#include <string>

#include <dune/grid/io/file/vtk/vtksequencewriter.hh>
#include <dumux/discretization/method.hh>

namespace Dumux {

/*!
 * \ingroup Newton
 * \brief A convergence writer interface
 * Provide an interface that show the minimal requirements a convergence write passed to a newton method has to fulfill
 * \note This is used together with a Newton solver, see documentation of the Newton solver for
 *       more information on how to use this class.
 */
template <class SolutionVector, class ResidualVector>
struct ConvergenceWriterInterface
{
    virtual ~ConvergenceWriterInterface() = default;

    virtual void write(const SolutionVector &uLastIter, const ResidualVector &deltaU, const ResidualVector &residual) {}
};

/*!
 * \ingroup Newton
 * \brief Writes the intermediate solutions for every Newton iteration
 * \note This is used together with a Newton solver, see documentation of the Newton solver for
 *       more information on how to use this class.
 */
template <class GridGeometry, class SolutionVector, class ResidualVector>
class NewtonConvergenceWriter : public ConvergenceWriterInterface<SolutionVector, ResidualVector>
{
    using GridView = typename GridGeometry::GridView;
    static constexpr auto numEq = SolutionVector::block_type::dimension;
    using Scalar = typename SolutionVector::block_type::value_type;

    static_assert(GridGeometry::discMethod != DiscretizationMethods::staggered,
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

        if (GridGeometry::discMethod == DiscretizationMethods::box)
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
               const ResidualVector& deltaU,
               const ResidualVector& residual) override
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
