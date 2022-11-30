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
 * \ingroup InputOutput
 * \brief Generic writer for time series files.
 */
#ifndef DUMUX_IO_TIME_SERIES_WRITER_HH
#define DUMUX_IO_TIME_SERIES_WRITER_HH

#include <config.h>

#if HAVE_GRIDFORMAT
#include <string>
#include <utility>
#include <concepts>
#include <type_traits>
#include <array>
#include <ranges>
#include <vector>

#include <gridformat/grid/adapters/dune.hpp>
#include <gridformat/gridformat.hpp>

#include <dune/common/timer.hh>
#include <dumux/discretization/method.hh>

namespace Dumux {
namespace TimeSeriesWriterDetail {

#if HAVE_MPI
inline constexpr bool haveMPI_ = true;
#else
inline constexpr bool haveMPI_ = false;
#endif

// TODO: This could be put into some common place or we expose dofCodim
//       as a static constexpr int of the GGs; However, a clear dof to
//       codim mapping probably doesn't exist for all schemes (e.g. higher-order FEM)
//       For such schemes, setDofField should maybe not be accessible. Alternatively,
//       we could use higher-order cell types, which would require some more involved
//       way of defining the points of the grids and the connectivity...
template<typename GG> struct DofCodim;
template<typename GG> requires(GG::discMethod == DiscretizationMethods::cctpfa)
struct DofCodim<GG> : public std::integral_constant<int, 0> {};
template<typename GG> requires(GG::discMethod == DiscretizationMethods::ccmpfa)
struct DofCodim<GG> : public std::integral_constant<int, 0> {};
template<typename GG> requires(GG::discMethod == DiscretizationMethods::box)
struct DofCodim<GG> : public std::integral_constant<int, GG::GridView::dimension> {};

template<typename Format> struct IsTimeSeriesFormat : public std::false_type {};
template<typename F> struct IsTimeSeriesFormat<GridFormat::FileFormat::PVD<F>> : public std::true_type {};
template<typename F> struct IsTimeSeriesFormat<GridFormat::FileFormat::TimeSeries<F>> : public std::true_type {};
template<typename Format>
inline constexpr bool isTimeSeriesFormat = IsTimeSeriesFormat<Format>::value;

template<typename T>
struct To;
template<GridFormat::Concepts::Scalar T>
struct To<T>
{
    static constexpr T from(const std::convertible_to<T> auto in)
    { return static_cast<T>(in); }
};
template<GridFormat::Concepts::Scalar T, std::size_t n>
struct To<std::array<T, n>>
{
    static constexpr auto from(const GridFormat::Concepts::MDRange<1> auto& vector)
    { return GridFormat::Ranges::to_array<n, T>(vector); }
};
template<GridFormat::Concepts::Scalar T, std::size_t m, std::size_t n>
struct To<std::array<std::array<T, n>, m>>
{
    static constexpr auto from(const GridFormat::Concepts::MDRange<2> auto& tensor)
    {
        int dir = 0;
        std::array<std::array<T, n>, m> result;
        std::ranges::for_each(tensor, [&] (const auto& tensorRow) {
            result[dir++] = GridFormat::Ranges::to_array<n, T>(tensorRow);
        });
        return result;
    }
};

} // namespace TimeSeriesWriterDetail


namespace IO::Format { using namespace GridFormat::FileFormat; }
namespace IO::Encoding { using namespace GridFormat::Encoding; }
namespace IO::Compression { using namespace GridFormat::Compression; using GridFormat::none; }

/*!
 * \ingroup InputOutput
 * \brief Generic writer for time series files.
 *        Per default, this writes .pvd files, which for now
 *        also is the only available option. XDMF-support is
 *        currently ongoing.
 */
template<typename GridVariables, typename SolutionVector>
class TimeSeriesWriter
{
    using GridGeometry = typename GridVariables::GridGeometry;
    using GridView = typename GridGeometry::GridView;
    static constexpr int dimWorld = GridView::dimensionworld;
    static constexpr int dofCodim = TimeSeriesWriterDetail::DofCodim<GridGeometry>::value;

    template<typename F>
    static constexpr bool isPointFunction = GridFormat::Concepts::PointFunction<F, GridView>;
    template<typename F>
    static constexpr bool isCellFunction = GridFormat::Concepts::CellFunction<F, GridView>;

    using VolumeVariables = typename GridVariables::GridVolumeVariables::VolumeVariables;
    using Writer = GridFormat::Writer<GridView>;

    using Scalar = typename GridVariables::Scalar;
    using Vector = std::array<Scalar, dimWorld>;
    using Tensor = std::array<Vector, dimWorld>;

    class VolVarFunctionStorage;

 public:
    TimeSeriesWriter(const GridVariables& gridVars,
                     const SolutionVector& x,
                     const std::string& filename)
    : TimeSeriesWriter(gridVars, x, filename, GridFormat::pvd(GridFormat::default_for<GridView>()))
    {}

    template<typename Format> requires(!TimeSeriesWriterDetail::isTimeSeriesFormat<Format>)
    TimeSeriesWriter(const GridVariables& gridVars,
                     const SolutionVector& x,
                     const std::string& filename,
                     const Format& f)
    : TimeSeriesWriter(gridVars, x, filename, GridFormat::pvd(f))
    {}

    template<typename Format> requires(TimeSeriesWriterDetail::isTimeSeriesFormat<Format>)
    TimeSeriesWriter(const GridVariables& gridVars,
                     const SolutionVector& x,
                     const std::string& filename,
                     const Format& f)
    : gridVariables_(gridVars)
    , solutionVector_(x)
    {
        if (gridView_().comm().size() == 1)
            writer_ = std::make_unique<Writer>(f, gridView_(), filename);
        else if constexpr (TimeSeriesWriterDetail::haveMPI_)
            writer_ = std::make_unique<Writer>(f, gridView_(), getComm_(), filename);
        else
            DUNE_THROW(Dune::NotImplemented, "Communicators other than MPI");
    }

    template<GridFormat::Concepts::PointFunction<GridView> F>
    void setPointField(const std::string& name, F&& field)
    { writer_->set_point_field(name, std::forward<F>(field)); }

    template<GridFormat::Concepts::CellFunction<GridView> F>
    void setCellField(const std::string& name, F&& field)
    { writer_->set_cell_field(name, std::forward<F>(field)); }

    template<std::invocable<const VolumeVariables&> F>
    void setVolVarField(const std::string& name, F&& field)
    { volVarFunctions_.add(name, std::forward<F>(field)); }

    template<typename F>
    void setDofField(const std::string& name, F&& field)
    {
        if constexpr (dofCodim == 0) {
            static_assert(isCellFunction<F>, "The GridGeometry requires dof functions to be callable with elements.");
            setCellField(name, std::forward<F>(field));
        }
        else if constexpr (dofCodim == GridView::dimension) {
            static_assert(isPointFunction<F>, "The GridGeometry requires dof functions to be callable with vertices.");
            setPointField(name, std::forward<F>(field));
        }
        else
            DUNE_THROW(Dune::NotImplemented, "Field output for edges or faces");
    }

    void write(std::floating_point auto time_value)
    {
        std::unordered_map<std::string, std::vector<Scalar>> scalars;
        std::unordered_map<std::string, std::vector<Vector>> vectors;
        std::unordered_map<std::string, std::vector<Tensor>> tensors;
        fillVolVarData_(scalars, vectors, tensors);

        for (const auto& [name, values] : scalars) setDofValues_(name, values);
        for (const auto& [name, values] : vectors) setDofValues_(name, values);
        for (const auto& [name, values] : tensors) setDofValues_(name, values);
        setCellField("rank", [&] (const auto&) { return gridView_().comm().rank(); });

        Dune::Timer timer;
        const std::string filename = writer_->write(time_value);
        timer.stop();

        if (gridView_().comm().rank() == 0)
            std::cout << "Writing '" << filename << "' took " << timer.elapsed() << " seconds" << std::endl;
    }

 private:
    const GridView& gridView_() const { return gridGeometry_().gridView(); }
    const GridGeometry& gridGeometry_() const { return gridVariables_.gridGeometry(); }
    auto getComm_() requires(TimeSeriesWriterDetail::haveMPI_)
    {
        MPI_Comm comm = gridView_().comm();
        return comm;
    }

    void fillVolVarData_(std::unordered_map<std::string, std::vector<Scalar>>& scalars,
                         std::unordered_map<std::string, std::vector<Vector>>& vectors,
                         std::unordered_map<std::string, std::vector<Tensor>>& tensors)
    {
        volVarFunctions_.setSize(scalars, gridGeometry_().numDofs());
        volVarFunctions_.setSize(vectors, gridGeometry_().numDofs());
        volVarFunctions_.setSize(tensors, gridGeometry_().numDofs());
        for (const auto& element : elements(gridView_(), Dune::Partitions::interior))
        {
            auto fvGeometry = localView(gridGeometry_()).bindElement(element);
            auto elemVolVars = localView(gridVariables_.curGridVolVars());
            elemVolVars.bindElement(element, fvGeometry, solutionVector_);
            for (const auto& scv : scvs(fvGeometry))
            {
                const auto& volVars = elemVolVars[scv];
                volVarFunctions_.setValues(scalars, volVars, scv.dofIndex());
                volVarFunctions_.setValues(vectors, volVars, scv.dofIndex());
                volVarFunctions_.setValues(tensors, volVars, scv.dofIndex());
            }
        }
    }

    template<typename Values>
    void setDofValues_(const std::string& name, const Values& values)
    {
        if constexpr (dofCodim == 0)
            writer_->set_cell_field(name, [&] (const auto& e) {
                return values[gridGeometry_().elementMapper().index(e)];
            });
        else if constexpr (dofCodim == GridView::dimension)
            writer_->set_point_field(name, [&] (const auto& v) {
                return values[gridGeometry_().vertexMapper().index(v)];
            });
        else
            DUNE_THROW(Dune::NotImplemented, "Edge/Face output");
    }

    const GridVariables& gridVariables_;
    const SolutionVector& solutionVector_;

    std::unique_ptr<Writer> writer_{nullptr};
    VolVarFunctionStorage volVarFunctions_;
};

template<typename GridVariables, typename SolutionVector>
class TimeSeriesWriter<GridVariables, SolutionVector>::VolVarFunctionStorage
{
    template<typename T> struct StoredType;
    template<GridFormat::Concepts::Scalar T> struct StoredType<T> : public std::type_identity<Scalar> {};
    template<GridFormat::Concepts::MDRange<2> T> struct StoredType<T> : public std::type_identity<Vector> {};
    template<GridFormat::Concepts::MDRange<3> T> struct StoredType<T> : public std::type_identity<Tensor> {};
    template<typename T, int size> struct StoredType<Dune::FieldVector<T, size>> : public std::type_identity<Vector> {};
    template<typename T, int m, int n> struct StoredType<Dune::FieldMatrix<T, m, n>> : public std::type_identity<Tensor> {};

    template<typename T>
    class Storage
    {
        using Function = std::function<T(const VolumeVariables&)>;

    public:
        template<std::invocable<const VolumeVariables&> F>
            requires(!std::is_lvalue_reference_v<F>)
        void set(const std::string& name, F&& f)
        {
            functions_[name] = [_f = std::move(f)] (const VolumeVariables& v) {
                return TimeSeriesWriterDetail::To<T>::from(std::invoke(_f, v));
            };
        }

        auto size() const { return functions_.size(); }
        auto begin() const { return functions_.begin(); }
        auto end() const { return functions_.end(); }

    private:
        std::unordered_map<std::string, Function> functions_;
    };

public:
    template<std::invocable<const VolumeVariables&> F>
    void add(const std::string& name, F&& field)
    {
        using T = std::invoke_result_t<F, const VolumeVariables&>;
        using S = typename StoredType<T>::type;
        std::get<Storage<S>>(storage_).set(name, std::forward<F>(field));
    }

    template<typename T>
    void setSize(std::unordered_map<std::string, std::vector<T>>& out, std::size_t size)
    {
        for (const auto& [name, _] : std::get<Storage<T>>(storage_))
            out[name].resize(size);
    }

    template<typename T>
    void setValues(std::unordered_map<std::string, std::vector<T>>& out,
                    const VolumeVariables& volVars,
                    std::size_t index) const
    {
        for (const auto& [name, function] : std::get<Storage<T>>(storage_))
            out[name][index] = function(volVars);
    }

private:
    template<int i, typename Action, typename... Args>
    void applyToAll_(const Action& action, std::tuple<Args...>& tuple) const
    {
        action(std::get<i>(tuple));
        if constexpr (i < sizeof...(Args) - 1)
            applyToAll_<i+1>(action, tuple);
    }

    std::tuple<Storage<Scalar>, Storage<Vector>, Storage<Tensor>> storage_;
};

} // end namespace Dumux

#else // HAVE_GRIDFORMAT

template<typename GridVariables>
class TimeSeriesWriter
{ static_assert(false, "GridFormat library required for TimeSeriesWriter"); };

#endif // HAVE_GRIDFORMAT
#endif
