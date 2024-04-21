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

#include <ranges>
#include <execution>
#include <type_traits>

#include <gridformat/gridformat.hpp>
#include <gridformat/traits/dune.hpp>

#include <dumux/discretization/method.hh>

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

/*!
 * \ingroup InputOutput
 * \brief Generic output module for entity & volume variable fields. Supports a variety of grid file formats.
 */
template<typename GridVariables, typename SolutionVector>
class OutputModule : private GridWriter<typename GridVariables::GridGeometry::GridView, 1> {
    using ParentType = GridWriter<typename GridVariables::GridGeometry::GridView, 1>;

    static constexpr bool isCVFE = DiscretizationMethods::isCVFE<typename GridVariables::GridGeometry::DiscretizationMethod>;
    static constexpr int dimWorld = GridVariables::GridGeometry::GridView::dimensionworld;
    using Scalar = typename GridVariables::Scalar;
    using Vector = Dune::FieldVector<Scalar, dimWorld>;
    using Tensor = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;
    using VolVar = typename GridVariables::VolumeVariables;

    class VolVarFieldStorage;

 public:
    using VolumeVariables = VolVar;

    static constexpr auto defaultFileFormat = IO::Format::pvd_with(
        IO::Format::vtu.with({
            .encoder = IO::Encoding::ascii,
            .compressor = IO::Compression::none,
            .data_format = VTK::DataFormat::inlined
        })
    );

    /*!
     * \brief Constructor for output to the default file format.
     */
    explicit OutputModule(const GridVariables& gridVariables,
                          const SolutionVector& sol,
                          const std::string& filename)
    : ParentType{defaultFileFormat, gridVariables.gridGeometry().gridView(), filename, order<1>}
    , gridVariables_{gridVariables}
    , solutionVector_{sol}
    {}

    /*!
     * \brief Constructor for stationary file formats.
     * \note This does not compile if the chosen format is not a stationary file format.
     */
    template<typename Format>
    explicit OutputModule(const Format& fmt,
                          const GridVariables& gridVariables,
                          const SolutionVector& sol)
    : ParentType{fmt, gridVariables.gridGeometry().gridView(), order<1>}
    , gridVariables_{gridVariables}
    , solutionVector_{sol}
    {}

    /*!
     * \brief Constructor for transient file formats, i.e. time series.
     * \note This does not compile if the chosen format is not a transient file format.
     */
    template<typename Format>
    explicit OutputModule(const Format& fmt,
                          const GridVariables& gridVariables,
                          const SolutionVector& sol,
                          const std::string& filename)
    : ParentType{fmt, gridVariables.gridGeometry().gridView(), filename, order<1>}
    , gridVariables_{gridVariables}
    , solutionVector_{sol}
    {}

    /*!
     * \brief Register a volume variable to be added to the output.
     */
    template<std::invocable<const VolumeVariables&> VolVarFunction>
    void addVolumeVariable(VolVarFunction&& f, const std::string& name)
    {
        using ResultType = std::invoke_result_t<std::remove_cvref_t<VolVarFunction>, const VolumeVariables&>;
        if constexpr (GridFormat::Concepts::Scalar<ResultType>)
            setVolVarField_<ResultType>(name, volVarFields_.registerScalarField(name, [_f=std::move(f)] (const auto& vv) {
                return static_cast<Scalar>(_f(vv));
            }));
        else if constexpr (GridFormat::mdrange_dimension<ResultType> == 1)
            setVolVarField_<GridFormat::MDRangeScalar<ResultType>>(name, volVarFields_.registerVectorField(name, [_f=std::move(f)] (const auto& vv) {
                return VolVarFieldStorage::toStorageVector(_f(vv));
            }));
        else if constexpr (GridFormat::mdrange_dimension<ResultType> == 2)
            setVolVarField_<GridFormat::MDRangeScalar<ResultType>>(name, volVarFields_.registerTensorField(name, [_f=std::move(f)] (const auto& vv) {
                return VolVarFieldStorage::toStorageTensor(_f(vv));
            }));
        else
        {
            static_assert(
                Dune::AlwaysFalse<VolVarFunction>::value,
                "Could not identify the given volume variable as scalar, vector or tensor."
            );
        }
    }

    /*!
     * \brief Register a dof field (automatically selects point or cell field depending on the scheme).
     */
    template<typename DofFunction>
    void addField(DofFunction&& f, const std::string& name)
    {
        if constexpr (isCVFE)
            addPointField(std::forward<DofFunction>(f), name);
        else
            addCellField(std::forward<DofFunction>(f), name);
    }

    /*!
     * \brief Register a point field.
     */
    template<typename DofFunction>
    void addPointField(DofFunction&& f, const std::string& name)
    { this->setPointField(name, std::forward<DofFunction>(f)); }

    /*!
     * \brief Register a point field.
     */
    template<typename DofFunction>
    void addCellField(DofFunction&& f, const std::string& name)
    { this->setCellField(name, std::forward<DofFunction>(f)); }

    /*!
     * \brief Write the registered fields into the file with the given name.
     */
    std::string write(const std::string& name)
    {
        volVarFields_.updateFieldData(gridVariables_, solutionVector_);
        auto filename = ParentType::write(name);
        volVarFields_.clearFieldData();
        return filename;
    }

    /*!
     * \brief Write a step in a time series.
     */
    template<std::floating_point T>
    std::string write(T time)
    {
        volVarFields_.updateFieldData(gridVariables_, solutionVector_);
        auto filename = ParentType::write(time);
        volVarFields_.clearFieldData();
        return filename;
     }

    //! clear all registered data
    void clear() {
        ParentType::clear();
        volVarFields_.clear();
    }

 private:
    template<typename ResultType, typename Id>
    void setVolVarField_(const std::string& name, Id&& volVarFieldId)
    {
        auto dofEntityField = [&, _id=std::move(volVarFieldId)] (const auto& entity) {
            return volVarFields_.getValue(_id, gridVariables_.gridGeometry().dofMapper().index(entity));
        };
        if constexpr (isCVFE)
            this->setPointField(name, std::move(dofEntityField), GridFormat::Precision<ResultType>{});
        else
            this->setCellField(name, std::move(dofEntityField), GridFormat::Precision<ResultType>{});
    }

    const GridVariables& gridVariables_;
    const SolutionVector& solutionVector_;
    VolVarFieldStorage volVarFields_;
};

// Class to store vol var fields; implementation detail of the OutputModule class
template<typename GridVariables, typename SolutionVector>
class OutputModule<GridVariables, SolutionVector>::VolVarFieldStorage
{
    enum class FieldType
    { scalar, vector, tensor };

    struct FieldInfo
    {
        std::string name;
        FieldType type;
        std::size_t index;
    };

    template<typename T>
    struct FieldStorage
    {
        std::vector<T> data;
        std::function<T(const VolVar&)> getter;
    };

public:
    template<FieldType ft>
    struct FieldId { std::size_t index; };

    template<std::ranges::range R>
    static constexpr auto toStorageVector(R&& in)
    {
        Vector result;
        std::ranges::copy(in, result.begin());
        return result;
    }

    template<GridFormat::Concepts::MDRange<2> R>
    static constexpr auto toStorageTensor(R&& in)
    {
        Tensor result;
        std::ranges::for_each(in, [&, i=0] (const auto& row) mutable {
            std::ranges::copy(row, result[i++].begin());
        });
        return result;
    }

    template<FieldType ft>
    const auto& getValue(const FieldId<ft>& id, std::size_t idx) const
    {
        if constexpr (ft == FieldType::scalar)
            return scalarFieldStorage_.at(id.index).data.at(idx);
        else if constexpr (ft == FieldType::vector)
            return vectorFieldStorage_.at(id.index).data.at(idx);
        else
            return tensorFieldStorage_.at(id.index).data.at(idx);
    }

    auto registerScalarField(std::string name, std::function<Scalar(const VolVar&)> f)
    { return register_<FieldType::scalar>(std::move(name), scalarFieldStorage_, std::move(f)); }

    auto registerVectorField(std::string name, std::function<Vector(const VolVar&)> f)
    { return register_<FieldType::vector>(std::move(name), vectorFieldStorage_, std::move(f)); }

    auto registerTensorField(std::string name, std::function<Tensor(const VolVar&)> f)
    { return register_<FieldType::tensor>(std::move(name), tensorFieldStorage_, std::move(f)); }

    void updateFieldData(const GridVariables& gridVars, const SolutionVector& x)
    {
        resizeFieldData_(gridVars.gridGeometry().numDofs());
        const auto range = GridFormat::cells(gridVars.gridGeometry().gridView());
        std::for_each(
            std::execution::par_unseq,
            std::ranges::begin(range),
            std::ranges::end(range),
            [&] (const auto& element) {
                auto fvGeometry = localView(gridVars.gridGeometry()).bindElement(element);
                auto elemVolVars = localView(gridVars.curGridVolVars()).bindElement(element, fvGeometry, x);
                for (const auto& scv : scvs(fvGeometry))
                {
                    const auto& volVars = elemVolVars[scv];
                    for (auto& s : scalarFieldStorage_) { s.data.at(scv.dofIndex()) = s.getter(volVars); }
                    for (auto& s : vectorFieldStorage_) { s.data.at(scv.dofIndex()) = s.getter(volVars); }
                    for (auto& s : tensorFieldStorage_) { s.data.at(scv.dofIndex()) = s.getter(volVars); }
                }
        });
    }

    void clearFieldData()
    {
        for (auto& s : scalarFieldStorage_) { s.data.clear(); }
        for (auto& s : vectorFieldStorage_) { s.data.clear(); }
        for (auto& s : tensorFieldStorage_) { s.data.clear(); }
    }

    void clear()
    {
        fields_.clear();
        scalarFieldStorage_.clear();
        vectorFieldStorage_.clear();
        tensorFieldStorage_.clear();
    }

private:
    template<FieldType ft, typename T>
    auto register_(std::string&& name,
                   std::vector<FieldStorage<T>>& storage,
                   std::function<T(const VolVar&)>&& f)
    {
        if (exists_<ft>(name))
            DUNE_THROW(Dune::InvalidStateException, "Volume variables field '" << name << "' is already defined.");

        FieldId<ft> id{storage.size()};
        fields_.emplace_back(FieldInfo{std::move(name), ft, id.index});
        storage.push_back({{}, std::move(f)});
        return id;
    }

    template<FieldType ft>
    bool exists_(const std::string& name) const
    {
        return std::ranges::any_of(fields_, [&] (const FieldInfo& info) {
            return info.type == ft && info.name == name;
        });
    }

    void resizeFieldData_(std::size_t size)
    {
        std::ranges::for_each(scalarFieldStorage_, [&] (auto& s) { s.data.resize(size); });
        std::ranges::for_each(vectorFieldStorage_, [&] (auto& s) { s.data.resize(size); });
        std::ranges::for_each(tensorFieldStorage_, [&] (auto& s) { s.data.resize(size); });
    }

    std::vector<FieldInfo> fields_;
    std::vector<FieldStorage<Scalar>> scalarFieldStorage_;
    std::vector<FieldStorage<Vector>> vectorFieldStorage_;
    std::vector<FieldStorage<Tensor>> tensorFieldStorage_;
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

template<class... Args>
class OutputModule
{
public:
    template<class... _Args>
    OutputModule(_Args&&...)
    {
        static_assert(
            false,
            "OutputModule only available when the GridFormat library is available. "
            "Use `git submodule update --init` to pull it."
        );
    }
};

} // namespace Dumux::IO

#endif // DUMUX_HAVE_GRIDFORMAT
#endif // DUMUX_IO_GRID_HH
