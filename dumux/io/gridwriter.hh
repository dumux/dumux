// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup InputOutput
 * \brief Generic writer for a variety of grid file formats.
 */
#ifndef DUMUX_IO_GRID_WRITER_HH
#define DUMUX_IO_GRID_WRITER_HH

#include <config.h>

#if DUMUX_HAVE_GRIDFORMAT

#include <type_traits>

#include <gridformat/gridformat.hpp>
#include <gridformat/traits/dune.hpp>

namespace Dumux::IO {


#ifndef DOXYGEN
namespace Detail {

template<typename Grid, typename Format, typename Comm, typename... Args>
auto makeParallelWriter(const Grid& grid, const Format& fmt, const Dune::Communication<Comm>& comm, Args&&... args)
{
    return comm.size() > 1
        ? GridFormat::Writer<Grid>{fmt, grid, static_cast<Comm>(comm), std::forward<Args>(args)...}
        : GridFormat::Writer<Grid>{fmt, grid, std::forward<Args>(args)...};
}

template<typename Grid, typename Format, typename... Args>
auto makeParallelWriter(const Grid& grid, const Format& fmt, const Dune::Communication<Dune::No_Comm>&, Args&&... args)
{ return GridFormat::Writer<Grid>{fmt, grid, std::forward<Args>(args)...}; }

template<typename Grid, typename Format, typename... Args>
auto makeWriter(const Grid& grid, const Format& fmt, Args&&... args)
{
    const auto& comm = GridFormat::Dune::Traits::GridView<Grid>::get(grid).comm();
    return makeParallelWriter(grid, fmt, comm, std::forward<Args>(args)...);
}

} // namespace Detail
#endif // DOXYGEN


namespace VTK { using namespace GridFormat::VTK; }
namespace Format { using namespace GridFormat::Formats; }
namespace Encoding { using namespace GridFormat::Encoding; }
namespace Compression { using namespace GridFormat::Compression; using GridFormat::none; }
namespace Precision {
    using GridFormat::float32;
    using GridFormat::float64;

    using GridFormat::uint64;
    using GridFormat::uint32;
    using GridFormat::uint16;
    using GridFormat::uint8;

    using GridFormat::int64;
    using GridFormat::int32;
    using GridFormat::int16;
    using GridFormat::int8;
} // namespace Precision


/*!
 * \ingroup InputOutput
 * \brief Represents the interpolation order with which fields are written out.
 */
template<int order>
struct Order { static_assert(order > 0, "order must be > 0"); };

template<int o>
inline constexpr auto order = Order<o>{};

/*!
 * \ingroup InputOutput
 * \brief Generic writer for a variety of grid file formats.
 *        Supports higher-order output of fields provided as Dune::Function objects.
 *        To create a writer for higher-order output, you may write
 *        \code
 *            GridWriter writer{gridView, order<2>};
 *        \endcode
 */
template<GridFormat::Concepts::Grid GridView, int order = 1>
class GridWriter
{
    using Grid = std::conditional_t<
        (order > 1),
        GridFormat::Dune::LagrangePolynomialGrid<GridView>,
        GridView
    >;
    using Cell = GridFormat::Cell<Grid>;
    using Vertex = typename GridView::template Codim<GridView::dimension>::Entity;
    using Element = typename GridView::template Codim<0>::Entity;
    using Coordinate = typename Element::Geometry::GlobalCoordinate;
    using Writer = GridFormat::Writer<Grid>;

 public:
    /*!
     * \brief Constructor for non-transient file formats.
     * \note This does not compile if the chosen format is a transient file format.
     */
    template<typename Format>
    explicit GridWriter(const Format& fmt,
                        const GridView& gridView,
                        const Order<order>& = {})
    : gridView_{gridView}
    , grid_{makeGrid_(gridView)}
    , writer_{Detail::makeWriter(grid_, fmt)}
    { writer_.set_meta_data("rank", gridView_.comm().rank()); }

    /*!
     * \brief Constructor for transient file formats, i.e. time series.
     * \note This does not compile if the chosen format is not a transient file format.
     */
    template<typename Format>
    explicit GridWriter(const Format& fmt,
                        const GridView& gridView,
                        const std::string& filename,
                        const Order<order>& = {})
    : gridView_{gridView}
    , grid_{makeGrid_(gridView)}
    , writer_{Detail::makeWriter(grid_, fmt, filename)}
    { writer_.set_meta_data("rank", gridView_.comm().rank()); }

    /*!
     * \brief Write the registered fields into the file with the given name.
     * \note This function throws if the writer was constructed for time series output.
     */
    std::string write(const std::string& name) const
    { return writer_.write(name); }

    /*!
     * \brief Write a step in a time series.
     * \note This function throws if the writer was not constructed for time series output.
     */
    template<std::floating_point T>
    std::string write(T time) const
    { return writer_.write(time); }

    //! Set a cell field via a lambda invoked with grid elements
    template<GridFormat::Concepts::CellFunction<GridView> F,
             GridFormat::Concepts::Scalar T = GridFormat::FieldScalar<std::invoke_result_t<F, Element>>>
    void setCellField(const std::string& name, F&& f, const GridFormat::Precision<T>& prec = {})
    { writer_.set_cell_field(name, std::move(f), prec); }

    //! Set a Dune::Function as cell field
    template<GridFormat::Dune::Concepts::Function<GridView> F>
    void setCellField(const std::string& name, F&& f)
    { GridFormat::Dune::set_cell_function(std::forward<F>(f), writer_, name); }

    //! Set a Dune::Function as cell field with custom precision
    template<typename F, GridFormat::Concepts::Scalar T>
    void setCellField(const std::string& name, F&& f, const GridFormat::Precision<T>& prec)
    { GridFormat::Dune::set_cell_function(std::forward<F>(f), writer_, name, prec); }

    //! Set a point field via a lambda invoked with grid vertices
    template<GridFormat::Concepts::PointFunction<GridView> F,
             GridFormat::Concepts::Scalar T = GridFormat::FieldScalar<std::invoke_result_t<F, Vertex>>>
    void setPointField(const std::string& name, F&& f, const GridFormat::Precision<T>& prec = {})
    {
        static_assert(order == 1, "Point lambdas can only be used for order == 1. Use Dune::Functions instead.");
        writer_.set_point_field(name, std::move(f), prec);
    }

    //! Set a Dune::Function as point field
    template<GridFormat::Dune::Concepts::Function<GridView> F>
    void setPointField(const std::string& name, F&& f)
    { GridFormat::Dune::set_point_function(std::forward<F>(f), writer_, name); }

    //! Set a Dune::Function as point data with custom precision
    template<GridFormat::Dune::Concepts::Function<GridView> F, GridFormat::Concepts::Scalar T>
    void setPointField(const std::string& name, F&& f, const GridFormat::Precision<T>& prec)
    { GridFormat::Dune::set_point_function(std::forward<F>(f), writer_, name, prec); }

    //! Clear all data
    void clear()
    { writer_.clear(); }

    //! Update (must be called when the grid has changed)
    void update()
    {
        if constexpr (order > 1)
            grid_.update(gridView_);
    }

 private:
    Grid makeGrid_(const GridView& gv) const
    {
        if constexpr (order > 1)
            return Grid{gv, order};
        else
            return gv;
    }

    GridView gridView_;
    Grid grid_;
    Writer writer_;
};

} // namespace Dumux::IO

#else // DUMUX_HAVE_GRIDFORMAT

namespace Dumux::IO {

template<class... Args>
class GridWriter
{
public:
    template<class... _Args>
    GridWriter(_Args&&...)
    {
        static_assert(
            false,
            "GridWriter only available when the GridFormat library is available. "
            "Use `git submodule update --init` to pull it."
        );
    }
};

} // namespace Dumux::IO

#endif // DUMUX_HAVE_GRIDFORMAT
#endif // DUMUX_IO_GRID_HH
