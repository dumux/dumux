// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup InputOutput
 * \brief Minimal Dune-Functions-compatible wrapper around a CVFE element solution
 *        for use with Dumux::IO::GridWriter at higher interpolation orders.
 *
 * Usage example:
 * \code
 *   GridWriter writer{Format::vtu, gridView, order<3>};
 *   writer.setPointField("velocity", Dumux::IO::cvfeGridFunction(gridGeometry, x));
 * \endcode
 *
 * The returned object satisfies `GridFormat::Dune::Concepts::Function<GridView>`:
 *   - `localFunction(f)` returns a local function object
 *   - local function has `bind(element)` and `operator()(LocalCoordinate)`
 */
#ifndef DUMUX_IO_CVFE_GRID_FUNCTION_HH
#define DUMUX_IO_CVFE_GRID_FUNCTION_HH

#include <optional>

#include <dumux/discretization/elementsolution.hh>
#include <dumux/discretization/evalsolution.hh>

namespace Dumux::IO {

/*!
 * \brief Thin dune-functions-compatible wrapper for a CVFE grid function.
 *
 * Satisfies `GridFormat::Dune::Concepts::Function<GridView>`:
 * - `localFunction(f)` returns a `LocalFunction` object
 * - `LocalFunction::bind(element)` binds to an element
 * - `LocalFunction::operator()(localCoord)` evaluates the FE solution
 *
 * \tparam GridGeometry  The finite-volume grid geometry (e.g. PQ3FVGridGeometry)
 * \tparam SolutionVector  The global solution vector type
 */
template<class GridGeometry, class SolutionVector>
class CVFEGridFunction
{
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;

public:
    /*!
     * \brief The local function type satisfying the dune-functions local function concept.
     */
    class LocalFunction
    {
    public:
        LocalFunction(const GridGeometry& gg, const SolutionVector& sol)
        : gridGeometry_(&gg), sol_(&sol)
        {}

        void bind(const Element& element)
        { element_.emplace(element); }

        void unbind()
        { element_.reset(); }

        //! Evaluate FE solution at local coordinate \p localPos in the bound element.
        auto operator()(const typename Element::Geometry::LocalCoordinate& localPos) const
        {
            const auto elemSol = elementSolution(*element_, *sol_, *gridGeometry_);
            return evalSolutionAtLocalPos(
                *element_, element_->geometry(),
                *gridGeometry_, elemSol, localPos
            );
        }

    private:
        const GridGeometry* gridGeometry_;
        const SolutionVector* sol_;
        std::optional<Element> element_;
    };

    CVFEGridFunction(const GridGeometry& gg, const SolutionVector& sol)
    : gridGeometry_(&gg), sol_(&sol)
    {}

    //! Returns a local function. Found via ADL as required by gridformat.
    friend LocalFunction localFunction(const CVFEGridFunction& f)
    { return LocalFunction(*f.gridGeometry_, *f.sol_); }

private:
    const GridGeometry* gridGeometry_;
    const SolutionVector* sol_;
};

/*!
 * \brief Convenience factory function for `CVFEGridFunction`.
 * \code
 *   writer.setPointField("velocity", Dumux::IO::cvfeGridFunction(gg, x));
 * \endcode
 */
template<class GridGeometry, class SolutionVector>
auto cvfeGridFunction(const GridGeometry& gg, const SolutionVector& sol)
{ return CVFEGridFunction<GridGeometry, SolutionVector>{gg, sol}; }

} // namespace Dumux::IO

#endif
