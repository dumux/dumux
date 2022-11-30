#include <config.h>
#include <cmath>
#include <memory>

#include <dune/common/fmatrix.hh>
#include <dune/grid/yaspgrid.hh>

#include <dumux/discretization/cellcentered/tpfa/fvgridgeometry.hh>
#include <dumux/io/timeserieswriter.hh>

template<typename Position>
double testFunction(const Position& pos)
{ return std::sin(10.0*pos[0])*std::cos(10.0*pos[1]); }

template<typename GG>
class MockGridVariables
{
    static constexpr int dim = GG::GridView::dimension;

    using Element = typename GG::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using FVGeometry = typename GG::LocalView;
    using SubControlVolume = typename GG::SubControlVolume;

public:
    using GridGeometry = GG;
    using Scalar = double;

    struct GridVolumeVariables
    {
        struct VolumeVariables
        {
            GlobalPosition position;

            std::uint8_t integral() const
            { return 1; }

            double scalar() const
            { return testFunction(position); }

            GlobalPosition vector() const
            {
                GlobalPosition result;
                std::fill(result.begin(), result.end(), testFunction(position));
                return result;
            }

            auto tensor() const
            {
                Dune::FieldMatrix<double, dim, dim> result;
                std::fill(result.begin(), result.end(), vector());
                return result;
            }
        };

        struct LocalView
        {
            template<typename SolutionVector>
            void bindElement(const Element& e,
                             const FVGeometry& fvg,
                             const SolutionVector&)
            {
                volVars_.clear();
                volVars_.resize(fvg.numScv());
                for (const auto& scv : scvs(fvg))
                    volVars_[scv.localDofIndex()] = VolumeVariables{scv.dofPosition()};
            }

            const VolumeVariables& operator[](const SubControlVolume& scv) const
            { return volVars_[scv.localDofIndex()]; }

        private:
            std::vector<VolumeVariables> volVars_;
        };

        friend LocalView localView(const GridVolumeVariables&)
        { return LocalView{}; }
    };

    explicit MockGridVariables(std::shared_ptr<GridGeometry> gg)
    : gridGeometry_(gg)
    {}

    const GridGeometry& gridGeometry() const
    { return *gridGeometry_; }

    GridVolumeVariables curGridVolVars() const
    { return {}; }

private:
    std::shared_ptr<GridGeometry> gridGeometry_;
};


using Grid = Dune::YaspGrid<2>;
using GridView = typename Grid::LeafGridView;


template<typename Writer, typename GridGeometry>
void writeWithFields(Writer&& writer, const GridGeometry& gridGeometry)
{
    using GridView = typename GridGeometry::GridView;
    std::vector<double> pdata(gridGeometry.gridView().size(GridView::dimension));
    for (const auto& v : vertices(gridGeometry.gridView()))
        pdata[gridGeometry.vertexMapper().index(v)] = testFunction(v.geometry().center());

    std::vector<double> cdata(gridGeometry.gridView().size(0));
    for (const auto& e : elements(gridGeometry.gridView()))
        cdata[gridGeometry.elementMapper().index(e)] = testFunction(e.geometry().center());

    writer.setPointField("pfunc", [] (const auto& vertex) {
        return testFunction(vertex.geometry().center());
    });
    writer.setPointField("pdata_from_vector", [&] (const auto& v) {
        return pdata[gridGeometry.vertexMapper().index(v)];
    });

    writer.setCellField("cfunc", [] (const auto& element) {
        return testFunction(element.geometry().center());
    });
    writer.setCellField("cdata_from_vector", [&] (const auto& e) {
        return cdata[gridGeometry.elementMapper().index(e)];
    });

    writer.setDofField("dfunc", [] (const auto& dof) {
        return testFunction(dof.geometry().center());
    });

    writer.setVolVarField("vv_scalar", [] (const auto& vv) { return vv.scalar(); });
    writer.setVolVarField("vv_vector", [] (const auto& vv) { return vv.vector(); });
    writer.setVolVarField("vv_tensor", [] (const auto& vv) { return vv.tensor(); });
    writer.setVolVarField("vv_integral", [] (const auto& vv) { return vv.integral(); });

    writer.write(1.0);
}

template<typename Format>
void write(const std::string& filename, std::optional<Format> format = {})
{
    Grid grid{{1.0, 1.0}, {10, 10}};
    grid.loadBalance();

    using GridGeometry = Dumux::CCTpfaFVGridGeometry<GridView>;
    auto gridGeometry = std::make_shared<GridGeometry>(grid.leafGridView());
    MockGridVariables gridVariables{gridGeometry};
    double dummySolutionVector;

    if (format)
        writeWithFields(
            Dumux::TimeSeriesWriter{gridVariables, dummySolutionVector, filename, *format},
            *gridGeometry
        );
    else
        writeWithFields(
            Dumux::TimeSeriesWriter{gridVariables, dummySolutionVector, filename},
            *gridGeometry
        );
}

template<typename Format>
void write(const std::string& filename, const Format& format)
{ write(filename, std::optional{format}); }


int main(int argc, char** argv) {
    Dune::MPIHelper::instance(argc, argv);

    // use default setting (compression - if available - & appended base64 encoding)
    // template arg is needed s.t. optional can be resolved & writer is created w/o format argument
    write<Dumux::IO::Format::VTU>("time_series_writer_default");

    // use default format, but set custom options
    write(
        "time_series_default_with_opts",
        Dumux::IO::Format::default_for<GridView>().with({.encoder = Dumux::IO::Encoding::ascii})
    );

    // explicitly use vtu with ascii encoding
    write(
        "time_series_writer_vtu_ascii",
        Dumux::IO::Format::vtu({.encoder = Dumux::IO::Encoding::ascii})
    );

    // use vti (available because we use YaspGrid) with raw encoding
    write(
        "time_series_writer_vti_raw",
        Dumux::IO::Format::vti({.encoder = Dumux::IO::Encoding::raw})
    );

    // use vtr (available because we use YaspGrid) without compression
    write(
        "time_series_writer_vtr",
        Dumux::IO::Format::vtr({.compressor = Dumux::IO::Compression::none})
    );

    // use vtp and fully specify both encoding and compression
    write(
        "time_series_writer_vtp",
        Dumux::IO::Format::vtp({
            .encoder = Dumux::IO::Encoding::base64,
            .compressor = Dumux::IO::Compression::none
        })
    );

    // we can also write a time-series directly without a .pvd file. This
    // simply writes a file per time step, adding the time information in the metadata
    write(
        "time_series_writer_vtu_series",
        Dumux::IO::Format::time_series(Dumux::IO::Format::vtu)
    );


#if HAVE_VTK_HDF
    // let's write into the VTK-HDF format
    write(
        "time_series_writer_hdf",
        Dumux::IO::Format::time_series(Dumux::IO::Format::vtk_hdf)
    );
#endif

    return 0;
}
