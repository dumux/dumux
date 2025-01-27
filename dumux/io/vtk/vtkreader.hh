// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup InputOutput
 * \brief A vtk file reader using tinyxml2 as xml backend
 */
#ifndef DUMUX_IO_VTK_VTKREADER_HH
#define DUMUX_IO_VTK_VTKREADER_HH

#include <iostream>
#include <iterator>
#include <algorithm>
#include <memory>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <filesystem>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/exceptions.hh>
#include <dune/grid/common/capabilities.hh>
#include <dune/grid/io/file/vtk/common.hh>
#include <dune/grid/common/gridfactory.hh>

#include <dumux/io/container.hh>

#if DUMUX_HAVE_GRIDFORMAT

#include <gridformat/gridformat.hpp>
#include <gridformat/traits/dune.hpp>
#include <gridformat/reader.hpp>

namespace Dumux::Detail::VTKReader {

// factory adapter that fulfills the GridFormat::Concepts::GridFactory concept
template<typename Grid>
class GridFactoryAdapter : public GridFormat::Dune::GridFactoryAdapter<Grid>
{
    using ParentType = GridFormat::Dune::GridFactoryAdapter<Grid>;
public:
    GridFactoryAdapter(Dune::GridFactory<Grid>& factory)
    : ParentType(factory)
    {}

    void insert_cell(const GridFormat::CellType& ct, const std::vector<std::size_t>& corners)
    {
        // special treatment for polylines not supported by the Dune adapter supplied by gridformat
        if (ct == GridFormat::CellType::polyline && corners.size() > 2)
            for (std::size_t i = 1; i < corners.size(); ++i)
                ParentType::insert_cell(GridFormat::CellType::segment, {corners[i-1], corners[i]});

        else
            ParentType::insert_cell(ct, corners);
    }
};

// reader adapter that fulfills the GridFormat::GridReader concept
// and can read polylines as separate cells
class VTPReader : public GridFormat::GridReader {
 private:
    void _open(const std::string& filename, typename GridReader::FieldNames& fields) override {
        auto helper = GridFormat::VTK::XMLReaderHelper::make_from(filename, "PolyData");

        _num_points = GridFormat::from_string<std::size_t>(helper.get("PolyData/Piece").get_attribute("NumberOfPoints"));
        _num_verts = helper.get("PolyData/Piece").get_attribute_or(std::size_t{0}, "NumberOfVerts");
        _num_polylines = helper.get("PolyData/Piece").get_attribute_or(std::size_t{0}, "NumberOfLines");
        _num_segments = 0; // initialize to zero and compute later
        _num_strips = helper.get("PolyData/Piece").get_attribute_or(std::size_t{0}, "NumberOfStrips");
        _num_polys = helper.get("PolyData/Piece").get_attribute_or(std::size_t{0}, "NumberOfPolys");

        if (_num_strips > 0)
            throw GridFormat::NotImplemented("Triangle strips are not (yet) supported");

        GridFormat::VTK::XMLDetail::copy_field_names_from(helper.get("PolyData"), fields);
        _helper.emplace(std::move(helper));

        if (_num_polylines > 0)
            _visit_cells("Lines", GridFormat::CellType::polyline, _num_polylines, // count segments in polylines
                [&](const GridFormat::CellType& ct, const std::vector<std::size_t>& corners) {
                    _num_segments += corners.size() - 1;
                }
            );
    }

    void _close() override {
        _helper.reset();
        _num_points = 0;
        _num_verts = 0;
        _num_polylines = 0;
        _num_segments = 0;
        _num_strips = 0;
        _num_polys = 0;
    }

    std::string _name() const override {
        return "VTPReader";
    }

    std::size_t _number_of_cells() const override {
        return _num_verts + _num_segments + _num_strips + _num_polys;
    }

    std::size_t _number_of_points() const override {
        return _num_points;
    }

    std::size_t _number_of_pieces() const override {
        return 1;
    }

    bool _is_sequence() const override {
        return false;
    }

    GridFormat::FieldPtr _points() const override {
        return _helper.value().make_points_field("PolyData/Piece/Points", _number_of_points());
    }

    void _visit_cells(const typename GridFormat::GridReader::CellVisitor& visitor) const override {
        if (_num_verts > 0)
            _visit_cells("Verts", GridFormat::CellType::vertex, _num_verts, visitor);
        if (_num_polylines > 0)
            _visit_cells("Lines", GridFormat::CellType::polyline, _num_polylines, visitor);
        if (_num_polys > 0)
            _visit_cells("Polys", GridFormat::CellType::polygon, _num_polys, visitor);
    }

    GridFormat::FieldPtr _cell_field(std::string_view name) const override {
        const auto num_vtk_cells = _num_verts + _num_polylines + _num_strips + _num_polys;
        const auto field_ptr = _helper.value().make_data_array_field(name, "PolyData/Piece/CellData", num_vtk_cells);
        if (_num_polylines == 0)
            return field_ptr;

        // make a field that wraps the usual field and repeats entries for polylines
        std::vector<std::size_t> extents(field_ptr->layout().dimension(), _number_of_cells());
        std::ranges::copy(field_ptr->layout() | std::views::drop(1), extents.begin() + 1);
        const auto precision = field_ptr->precision();

        return make_field_ptr(GridFormat::LazyField{
            std::move(field_ptr),
            GridFormat::MDLayout{extents},
            precision,
            [this, _nv=_number_of_cells()] (const auto& fp) {
                GridFormat::Serialization result{_nv*fp->precision().size_in_bytes()};
                std::size_t in = 0, out = 0;
                fp->precision().visit([&] <typename T> (const GridFormat::Precision<T>& p) {
                    const auto source = fp->template export_to<std::vector<T>>();
                    auto target = result.as_span_of(p);
                    _visit_cells([&] (const auto& ct, const auto& corners) {
                        if (ct == GridFormat::CellType::polyline) { // repeat source value for each segment
                            for (const auto endIdx = out + corners.size()-1; out < endIdx; ++out)
                                target[out] = source[in];
                            ++in;
                        } else {
                            target[out++] = source[in++];
                        }
                    });
                });

                return result;
            }
        });
    }

    GridFormat::FieldPtr _point_field(std::string_view name) const override {
        return _helper.value().make_data_array_field(name, "PolyData/Piece/PointData", _number_of_points());
    }

    GridFormat::FieldPtr _meta_data_field(std::string_view name) const override {
        return _helper.value().make_data_array_field(name, "PolyData/FieldData");
    }

    void _visit_cells(std::string_view type_name,
                      const GridFormat::CellType& cell_type,
                      const std::size_t expected_size,
                      const typename GridFormat::GridReader::CellVisitor& visitor) const {
        const std::string path = "PolyData/Piece/" + std::string{type_name};
        const auto offsets = _helper.value().make_data_array_field(
            "offsets", path, expected_size
        )->template export_to<std::vector<std::size_t>>();
        const auto connectivity = _helper.value().make_data_array_field(
            "connectivity", path
        )->template export_to<std::vector<std::size_t>>();

        std::vector<std::size_t> corners;
        for (std::size_t i = 0; i < expected_size; ++i) {
            corners.clear();
            const std::size_t offset_begin = (i == 0 ? 0 : offsets.at(i-1));
            const std::size_t offset_end = offsets.at(i);
            if (connectivity.size() < offset_end)
                throw GridFormat::SizeError("Connectivity array read from the file is too small");
            if (offset_end < offset_begin)
                throw GridFormat::ValueError("Invalid offset array");
            std::copy_n(connectivity.begin() + offset_begin, offset_end - offset_begin, std::back_inserter(corners));

            // look for cell types that allow for a more specific classification
            const GridFormat::CellType ct = cell_type == GridFormat::CellType::polygon ? (
                    corners.size() == 3 ? GridFormat::CellType::triangle
                                        : (corners.size() == 4 ? GridFormat::CellType::quadrilateral : cell_type)
                ) : (cell_type == GridFormat::CellType::polyline ? (
                    corners.size() == 2 ? GridFormat::CellType::segment
                                        : cell_type
                ) : cell_type);

            visitor(ct, corners);
        }
    }

    std::optional<GridFormat::VTK::XMLReaderHelper> _helper;
    std::size_t _num_points;
    std::size_t _num_verts;
    std::size_t _num_polylines;
    std::size_t _num_segments;
    std::size_t _num_strips;
    std::size_t _num_polys;
};

template<GridFormat::Concepts::Communicator C>
auto makeGridformatReaderFactory(const C& c) {
    return [fac = GridFormat::AnyReaderFactory<C>{c}] (const std::string& f)
    -> std::unique_ptr<GridFormat::GridReader>
    {
        // custom reader for vtp files
        if (f.ends_with(".vtp"))
            return std::make_unique<::Dumux::Detail::VTKReader::VTPReader>();
        return fac.make_for(f);
    };
}

} // end namespace Dumux::Detail::VTKReader

namespace Dumux {

/*!
 * \ingroup InputOutput
 * \brief A vtk file reader using gridformat
 */
class VTKReader
{
public:
    enum class DataType { cellData, pointData, fieldData };

    //! the cell / point data type for point data read from a grid file
    using Data = std::unordered_map<std::string, std::vector<double>>;

    explicit VTKReader(const std::string& fileName)
    : reader_(GridFormat::Reader::from(fileName,
        Detail::VTKReader::makeGridformatReaderFactory(
            Dune::MPIHelper::instance().getCommunicator()
        )
    ))
    {}

    /*!
     * \brief Reviews data from the vtk file to check if there is a data array with a specified name
     * \param name the name attribute of the data array to read
     * \param type the data array type
     */
    bool hasData(const std::string& name, const DataType& type) const
    {
        if (type == DataType::cellData)
            return std::ranges::any_of(
                cell_field_names(reader_),
                [&] (const auto& n) { return name == n; }
            );
        else if (type == DataType::pointData)
            return std::ranges::any_of(
                point_field_names(reader_),
                [&] (const auto& n) { return name == n; }
            );
        else if (type == DataType::fieldData)
            return std::ranges::any_of(
                meta_data_field_names(reader_),
                [&] (const auto& n) { return name == n; }
            );
        else
            DUNE_THROW(Dune::IOError, "Unknown data type for field '" << name << "'");
    }

    /*!
     * \brief read data from the vtk file to a container, e.g. std::vector<double>
     * \tparam Container a container type that has begin(), end(), push_back(), e.g. std::vector<>
     * \param name the name attribute of the data array to read
     * \param type the data array type
     */
    template<class Container>
    Container readData(const std::string& name, const DataType& type) const
    {
        if (type == DataType::cellData)
        {
            Container values(reader_.number_of_cells());
            reader_.cell_field(name)->export_to(values);
            return values;
        }
        else if (type == DataType::pointData)
        {
            Container values(reader_.number_of_points());
            reader_.point_field(name)->export_to(values);
            return values;
        }
        else if (type == DataType::fieldData)
        {
            DUNE_THROW(Dune::NotImplemented, "Reading meta data not yet implemented");
        }
        else
            DUNE_THROW(Dune::IOError, "Unknown data type for field '" << name << "'");
    }

    /*!
     * \brief Read a grid from a vtk/vtu/vtp file, ignoring cell and point data
     * \param verbose if the output should be verbose
     */
    template<class Grid>
    std::unique_ptr<Grid> readGrid(bool verbose = false) const
    {
        Dune::GridFactory<Grid> factory;
        return readGrid(factory, verbose);
    }

    /*!
     * \brief Read a grid from a vtk/vtu/vtp file, ignoring cell and point data
     * \note use this signature if the factory might be needed outside to interpret the data via the factory's insertion indices
     * \param verbose if the output should be verbose
     * \param factory the (empty) grid factory
     */
    template<class Grid>
    std::unique_ptr<Grid> readGrid(Dune::GridFactory<Grid>& factory, bool verbose = false) const
    {
        if (Dune::MPIHelper::instance().rank() == 0)
        {
            if (verbose)
                std::cout << "Reading " << Grid::dimension << "d grid from vtk file " << reader_.filename() << "." << std::endl;

            {
                Dumux::Detail::VTKReader::GridFactoryAdapter adapter{ factory };
                reader_.export_grid(adapter);
            }
        }

        return std::unique_ptr<Grid>(factory.createGrid());
    }

    /*!
     * \brief Read a grid from a vtk/vtu/vtp file, reading all cell and point data
     * \note the factory will be needed outside to interpret the data via the factory's insertion indices
     * \param factory the (empty) grid factory
     * \param cellData the cell data arrays to be filled
     * \param pointData the point data arrays to be filled
     * \param verbose if the output should be verbose
     */
    template<class Grid>
    std::unique_ptr<Grid> readGrid(Dune::GridFactory<Grid>& factory, Data& cellData, Data& pointData, bool verbose = false) const
    {
        auto grid = readGrid(factory, verbose);
        if (Dune::MPIHelper::instance().rank() == 0)
        {
            for (const auto& [name, field_ptr] : cell_fields(reader_))
                field_ptr->export_to(cellData[name]);

            for (const auto& [name, field_ptr] : point_fields(reader_))
                field_ptr->export_to(pointData[name]);
        }
        return std::unique_ptr<Grid>(std::move(grid));
    }

private:
    GridFormat::Reader reader_;
};

} // namespace Dumux

#else // DUMUX_HAVE_GRIDFORMAT

#include <dumux/io/xml/tinyxml2.h>
// fallback to simple vtk reader using tinyxml2
// this reader only support ascii files and limited VTK file formats

namespace Dumux::Detail::VTKReader {

/*!
 * \brief Get the piece node an xml document
 * \note Returns nullptr if the piece node wasn't found
 * \param doc an xml document
 * \param fileName a file name the doc was created from
 */
const tinyxml2::XMLElement* getPieceNode(const tinyxml2::XMLDocument& doc, const std::string& fileName)
{
    using namespace tinyxml2;

    const XMLElement* pieceNode = doc.FirstChildElement("VTKFile");
    if (pieceNode == nullptr)
        DUNE_THROW(Dune::IOError, "Couldn't get 'VTKFile' node in " << fileName << ".");

    pieceNode = pieceNode->FirstChildElement("UnstructuredGrid");
    if (pieceNode == nullptr)
        pieceNode = doc.FirstChildElement("VTKFile")->FirstChildElement("PolyData");
    if (pieceNode == nullptr)
        pieceNode = doc.FirstChildElement("VTKFile")->FirstChildElement("PUnstructuredGrid");
    if (pieceNode == nullptr)
        pieceNode = doc.FirstChildElement("VTKFile")->FirstChildElement("PPolyData");
    if (pieceNode == nullptr)
        DUNE_THROW(Dune::IOError, "Couldn't get 'UnstructuredGrid', 'PUnstructuredGrid', 'PolyData', or 'PPolyData' node in " << fileName << ".");

    return pieceNode->FirstChildElement("Piece");
}

/*!
 * \brief get the vtk filename for the current processor
 */
std::string getProcessPieceFileName(const std::string& pvtkFileName)
{
    using namespace tinyxml2;

    XMLDocument pDoc;
    const auto eResult = pDoc.LoadFile(pvtkFileName.c_str());
    if (eResult != XML_SUCCESS)
        DUNE_THROW(Dune::IOError, "Couldn't open XML file " << pvtkFileName << ".");

    // get the first piece node
    const XMLElement* pieceNode = getPieceNode(pDoc, pvtkFileName);
    const auto myrank = Dune::MPIHelper::instance().rank();
    for (int rank = 0; rank < myrank; ++rank)
    {
        pieceNode = pieceNode->NextSiblingElement("Piece");
        if (pieceNode == nullptr)
            DUNE_THROW(Dune::IOError, "Couldn't find 'Piece' node for rank "
                                    << rank << " in " << pvtkFileName << ".");
    }

    const char *vtkFileName = pieceNode->Attribute("Source");
    if (vtkFileName == nullptr)
        DUNE_THROW(Dune::IOError, "Couldn't get 'Source' attribute of 'Piece' node no. " << myrank << " in " << pvtkFileName);

    return vtkFileName;
}

} // end namespace Dumux::Detail::VTKReader

namespace Dumux {

/*!
 * \ingroup InputOutput
 * \brief A vtk file reader using tinyxml2 as xml backend
 */
class VTKReader
{
public:
    /*!
     * \brief The data array types
     */
    enum class DataType { cellData, pointData };

    //! the cell / point data type for point data read from a grid file
    using Data = std::unordered_map<std::string, std::vector<double>>;

    /*!
     * \brief The constructor creates a tinyxml2::XMLDocument from file
     */
    explicit VTKReader(const std::string& fileName)
    {
        using namespace tinyxml2;
        const auto ext = std::filesystem::path(fileName).extension().string();
        // If in parallel and the file to read is a parallel piece collection (pvtu/pvtp)
        // read only the piece belonging to the own process. For this to work, the files
        // have to have exactly the same amount of pieces than processes.
        fileName_ = Dune::MPIHelper::instance().size() > 1 && ext[1] == 'p' ?
            Detail::VTKReader::getProcessPieceFileName(fileName) : fileName;

        const auto eResult = doc_.LoadFile(fileName_.c_str());
        if (eResult != tinyxml2::XML_SUCCESS)
            DUNE_THROW(Dune::IOError, "Couldn't open XML file " << fileName_ << ".");

        const XMLElement* pieceNode = getPieceNode_();
        if (pieceNode == nullptr)
            DUNE_THROW(Dune::IOError, "Couldn't get 'Piece' node in " << fileName_ << ".");
    }

    /*!
     * \brief Reviews data from the vtk file to check if there is a data array with a specified name
     * \param name the name attribute of the data array to read
     * \param type the data array type
     */
    bool hasData(const std::string& name, const DataType& type) const
    {
        using namespace tinyxml2;

        const XMLElement* pieceNode = getPieceNode_();
        const XMLElement* dataNode = getDataNode_(pieceNode, type);
        if (dataNode == nullptr)
            return false;

        const XMLElement* dataArray = findDataArray_(dataNode, name);
        if (dataArray == nullptr)
            return false;

        return true;
    }

    /*!
     * \brief read data from the vtk file to a container, e.g. std::vector<double>
     * \tparam Container a container type that has begin(), end(), push_back(), e.g. std::vector<>
     * \param name the name attribute of the data array to read
     * \param type the data array type
     */
    template<class Container>
    Container readData(const std::string& name, const DataType& type) const
    {
        using namespace tinyxml2;

        const XMLElement* pieceNode = getPieceNode_();
        const XMLElement* dataNode = getDataNode_(pieceNode, type);
        if (dataNode == nullptr)
            DUNE_THROW(Dune::IOError, "Couldn't get 'PointData' or 'CellData' node in " << fileName_ << ".");

        const XMLElement* dataArray = findDataArray_(dataNode, name);
        if (dataArray == nullptr)
            DUNE_THROW(Dune::IOError, "Couldn't find the data array " << name << ".");

        return parseDataArray_<Container>(dataArray);
    }

    /*!
     * \brief Read a grid from a vtk/vtu/vtp file, ignoring cell and point data
     * \param verbose if the output should be verbose
     */
    template<class Grid>
    std::unique_ptr<Grid> readGrid(bool verbose = false) const
    {
        static_assert(!Dune::Capabilities::isCartesian<Grid>::v, "Grid reader only supports unstructured grid implementations");

        // make a grid factory
        Dune::GridFactory<Grid> factory;

        // only read on rank 0
        if (Dune::MPIHelper::instance().rank() == 0)
        {
            if (verbose)
                std::cout << "Reading " << Grid::dimension << "d grid from vtk file " << fileName_ << "." << std::endl;
            readGrid_(factory, verbose);
        }

        return std::unique_ptr<Grid>(factory.createGrid());
    }

    /*!
     * \brief Read a grid from a vtk/vtu/vtp file, ignoring cell and point data
     * \note use this signature if the factory might be needed outside to interpret the data via the factory's insertion indices
     * \param verbose if the output should be verbose
     * \param factory the (empty) grid factory
     */
    template<class Grid>
    std::unique_ptr<Grid> readGrid(Dune::GridFactory<Grid>& factory, bool verbose = false) const
    {
        static_assert(!Dune::Capabilities::isCartesian<Grid>::v, "Grid reader only supports unstructured grid implementations");

        if (Dune::MPIHelper::instance().rank() == 0)
        {
            if (verbose)
                std::cout << "Reading " << Grid::dimension << "d grid from vtk file " << fileName_ << "." << std::endl;
            readGrid_(factory, verbose);
        }

        return std::unique_ptr<Grid>(factory.createGrid());
    }

    /*!
     * \brief Read a grid from a vtk/vtu/vtp file, reading all cell and point data
     * \note the factory will be needed outside to interpret the data via the factory's insertion indices
     * \param factory the (empty) grid factory
     * \param cellData the cell data arrays to be filled
     * \param pointData the point data arrays to be filled
     * \param verbose if the output should be verbose
     */
    template<class Grid>
    std::unique_ptr<Grid> readGrid(Dune::GridFactory<Grid>& factory, Data& cellData, Data& pointData, bool verbose = false) const
    {
        static_assert(!Dune::Capabilities::isCartesian<Grid>::v, "Grid reader only supports unstructured grid implementations");

        if (Dune::MPIHelper::instance().rank() == 0)
        {
            if (verbose)
                std::cout << "Reading " << Grid::dimension << "d grid from vtk file " << fileName_ << "." << std::endl;
            readGrid_(factory, verbose);
            readGridData_(cellData, pointData, verbose);
        }

        return std::unique_ptr<Grid>(factory.createGrid());
    }

private:
    /*!
     * \brief Read a grid from a vtk/vtu/vtp file
     * \param factory the (empty) grid factory
     * \param verbose if the output should be verbose
     */
    template<class Grid>
    void readGrid_(Dune::GridFactory<Grid>& factory, bool verbose = false) const
    {
        using namespace tinyxml2;

        const XMLElement* pieceNode = getPieceNode_();
        const XMLElement* pointsNode = pieceNode->FirstChildElement("Points")->FirstChildElement("DataArray");
        if (pointsNode == nullptr)
            DUNE_THROW(Dune::IOError, "Couldn't get data array of points in " << fileName_ << ".");

        using Point3D = Dune::FieldVector<double, 3>;
        std::vector<Point3D> points3D;
        std::stringstream dataStream(pointsNode->GetText());
        std::istream_iterator<Point3D> it(dataStream);
        std::copy(it, std::istream_iterator<Point3D>(), std::back_inserter(points3D));

        // adapt point dimensions if grid dimension is smaller than 3
        auto points = adaptPointDimension_<Grid::dimensionworld>(std::move(points3D));

        if (verbose) std::cout << "Found " << points.size() << " vertices." << std::endl;

        // insert vertices to the grid factory
        for (auto&& point : points)
            factory.insertVertex(std::move(point));

        const XMLElement* cellsNode = pieceNode->FirstChildElement("Cells");
        const XMLElement* linesNode = pieceNode->FirstChildElement("Lines");
        if (cellsNode)
        {
            const XMLElement* connectivityNode = findDataArray_(cellsNode, "connectivity");
            const XMLElement* offsetsNode = findDataArray_(cellsNode, "offsets");
            const XMLElement* typesNode = findDataArray_(cellsNode, "types");

            const auto connectivity = parseDataArray_<std::vector<unsigned int>>(connectivityNode);
            const auto offsets = parseDataArray_<std::vector<unsigned int>>(offsetsNode);
            const auto types = parseDataArray_<std::vector<unsigned int>>(typesNode);

            if (verbose) std::cout << "Found " << offsets.size() << " elements." << std::endl;

            unsigned int lastOffset = 0;
            for (unsigned int i = 0; i < offsets.size(); ++i)
            {
                const auto geomType = vtkToDuneGeomType_(types[i]);
                unsigned int offset = offsets[i];
                std::vector<unsigned int> corners; corners.resize(offset-lastOffset);
                for (unsigned int j = 0; j < offset-lastOffset; ++j)
                    corners[Dune::VTK::renumber(geomType, j)] = connectivity[lastOffset+j];
                factory.insertElement(geomType, std::move(corners));
                lastOffset = offset;
            }
        }
        // for poly data
        else if (linesNode)
        {
            // sanity check
            if (Grid::dimension != 1)
                DUNE_THROW(Dune::IOError, "Grid expects dimension " << Grid::dimension
                                           << " but " << fileName_ << " contains a 1D grid.");

            const XMLElement* connectivityNode = findDataArray_(linesNode, "connectivity");
            const XMLElement* offsetsNode = findDataArray_(linesNode, "offsets");

            const auto connectivity = parseDataArray_<std::vector<unsigned int>>(connectivityNode);
            const auto offsets = parseDataArray_<std::vector<unsigned int>>(offsetsNode);

            if (verbose) std::cout << "Found " << offsets.size() << " polylines." << std::endl;

            unsigned int lastOffset = 0;
            for (unsigned int i = 0; i < offsets.size(); ++i)
            {
                // a polyline can have many points in the VTK format
                // split the line in segments with two points
                unsigned int offset = offsets[i];
                for (unsigned int j = 0; j < offset-lastOffset-1; ++j)
                    factory.insertElement(Dune::GeometryTypes::line,
                        std::vector<unsigned int>({ connectivity[lastOffset+j], connectivity[lastOffset+j+1] }));
                lastOffset = offset;
            }
        }
        else
            DUNE_THROW(Dune::IOError, "No Cells or Lines element found in " << fileName_);
    }

    /*!
     * \brief Read a grid data from a vtk/vtu/vtp file
     * \param cellData the cell data arrays to be filled
     * \param pointData the point data arrays to be filled
     * \param verbose if the output should be verbose
     */
    void readGridData_(Data& cellData, Data& pointData, bool verbose = false) const
    {
        using namespace tinyxml2;

        const XMLElement* pieceNode = getPieceNode_();
        const XMLElement* cellDataNode = getDataNode_(pieceNode, DataType::cellData);
        if (cellDataNode != nullptr)
        {
            const XMLElement* cellsNode = pieceNode->FirstChildElement("Cells");
            const XMLElement* linesNode = pieceNode->FirstChildElement("Lines");

            if (cellsNode)
            {
                const XMLElement* dataArray = cellDataNode->FirstChildElement("DataArray");
                for (; dataArray != nullptr; dataArray = dataArray->NextSiblingElement("DataArray"))
                {
                    const char *attributeText = dataArray->Attribute("Name");

                    if (attributeText == nullptr)
                        DUNE_THROW(Dune::IOError, "Couldn't get Name attribute of a cell data array.");

                    cellData[std::string(attributeText)] = parseDataArray_<std::vector<double>>(dataArray);

                    if (verbose)
                        std::cout << "Read cell data field " << attributeText << std::endl;
                }
            }
            // for poly data
            else if (linesNode)
            {
                // first parse all the cell data (each cell in this sense can be a polyline)
                const XMLElement* dataArray = cellDataNode->FirstChildElement("DataArray");
                if (dataArray)
                {
                    Data polyLineCellData;
                    for (; dataArray != nullptr; dataArray = dataArray->NextSiblingElement("DataArray"))
                    {
                        const char *attributeText = dataArray->Attribute("Name");

                        if (attributeText == nullptr)
                            DUNE_THROW(Dune::IOError, "Couldn't get Name attribute of a cell data array.");

                        polyLineCellData[std::string(attributeText)] = parseDataArray_<std::vector<double>>(dataArray);

                        if (verbose)
                            std::cout << "Read cell data field " << attributeText << std::endl;
                    }

                    // a polyline can have many points in the VTK format
                    // we split the line in segments with two points
                    // so we also need to duplicate the cell data to fit the increased line number
                    const XMLElement* offsetsNode = findDataArray_(linesNode, "offsets");
                    const auto offsets = parseDataArray_<std::vector<unsigned int>>(offsetsNode);

                    if (offsets.size() != polyLineCellData.begin()->second.size())
                        DUNE_THROW(Dune::IOError, "Expected the same number of cell data entries (is "
                                                   << polyLineCellData.begin()->second.size()
                                                   << ") as polylines (" << offsets.size() << ")!");

                    // count the number of Dune cells to be able to resize the data vectors
                    unsigned int lastOffset = 0;
                    std::size_t numCells = 0;
                    for (unsigned int i = 0; i < offsets.size(); ++i)
                    {
                        unsigned int offset = offsets[i];
                        for (unsigned int j = 0; j < offset-lastOffset-1; ++j)
                            ++numCells;
                        lastOffset = offset;
                    }

                    // create the data arrays
                    for (const auto& dArray : polyLineCellData)
                    {
                        cellData[dArray.first] = std::vector<double>(numCells);
                        auto& cd = cellData[dArray.first];
                        const auto& pd = dArray.second;

                        lastOffset = 0;
                        std::size_t cellIdx = 0;
                        for (unsigned int i = 0; i < offsets.size(); ++i)
                        {
                            unsigned int offset = offsets[i];
                            for (unsigned int j = 0; j < offset-lastOffset-1; ++j)
                                cd[cellIdx++] = pd[i];
                            lastOffset = offset;
                        }
                    }
                }
            }
            else
                DUNE_THROW(Dune::IOError, "No Cells or Lines element found in " << fileName_);
        }

        const XMLElement* pointDataNode = getDataNode_(pieceNode, DataType::pointData);
        if (pointDataNode != nullptr)
        {
            const XMLElement* dataArray = pointDataNode->FirstChildElement("DataArray");
            for (; dataArray != nullptr; dataArray = dataArray->NextSiblingElement("DataArray"))
            {
                const char *attributeText = dataArray->Attribute("Name");

                if (attributeText == nullptr)
                    DUNE_THROW(Dune::IOError, "Couldn't get Name attribute of a point data array.");

                pointData[std::string(attributeText)] = parseDataArray_<std::vector<double>>(dataArray);

                if (verbose)
                    std::cout << "Read point data field " << attributeText << std::endl;
            }
        }
    }

    /*!
     * \brief Get the piece node of the XMLDocument
     * \note Returns nullptr if the piece node wasn't found
     */
    const tinyxml2::XMLElement* getPieceNode_() const
    { return Detail::VTKReader::getPieceNode(doc_, fileName_); }

    /*!
     * \brief Get the piece node of the XMLDocument
     * \param pieceNode the pieceNode of the vtk file
     * \param type the vtk data type (cell data or point data)
     * \note Returns nullptr if the data node wasn't found
     */
    const tinyxml2::XMLElement* getDataNode_(const tinyxml2::XMLElement* pieceNode, const DataType& type) const
    {
        using namespace tinyxml2;

        const XMLElement* dataNode = nullptr;
        if (type == DataType::pointData)
            dataNode = pieceNode->FirstChildElement("PointData");
        else if (type == DataType::cellData)
            dataNode = pieceNode->FirstChildElement("CellData");
        else
            DUNE_THROW(Dune::IOError, "Only cell and point data are supported.");

        return dataNode;
    }

    /*!
     * \brief Find a data array with a specific name
     * \param dataNode a cell or point data node
     * \param name the name of the data array to be found
     * \note Returns nullptr if the data array wasn't found
     */
    const tinyxml2::XMLElement* findDataArray_(const tinyxml2::XMLElement* dataNode, const std::string& name) const
    {
        using namespace tinyxml2;

        // loop over XML node siblings to find the correct data array
        const XMLElement* dataArray = dataNode->FirstChildElement("DataArray");
        for (; dataArray != nullptr; dataArray = dataArray->NextSiblingElement("DataArray"))
        {
            const char *attributeText = dataArray->Attribute("Name");

            if (attributeText == nullptr)
                DUNE_THROW(Dune::IOError, "Couldn't get Name attribute of a data array.");

            if (attributeText == name)
                break;
        }

        return dataArray;
    }

    /*!
     * \brief Parses the text of a data array into a container
     * \tparam Container a container type that has begin(), end(), push_back(), e.g. std::vector<double>
     * \param dataArray the data array node to be parsed
     */
    template<class Container>
    Container parseDataArray_(const tinyxml2::XMLElement* dataArray) const
    {
        std::stringstream dataStream(dataArray->GetText());
        return readStreamToContainer<Container>(dataStream);
    }

    /*!
     * \brief Return the Dune::GeometryType for a given VTK geometry type
     * \param vtkCellType the vtk cell type
     */
    Dune::GeometryType vtkToDuneGeomType_(unsigned int vtkCellType) const
    {
        switch (vtkCellType)
        {
            case Dune::VTK::GeometryType::vertex: return Dune::GeometryTypes::vertex;
            case Dune::VTK::GeometryType::line: return Dune::GeometryTypes::line;
            case Dune::VTK::GeometryType::triangle: return Dune::GeometryTypes::triangle;
            case Dune::VTK::GeometryType::quadrilateral: return Dune::GeometryTypes::quadrilateral;
            case Dune::VTK::GeometryType::tetrahedron: return Dune::GeometryTypes::tetrahedron;
            case Dune::VTK::GeometryType::hexahedron: return Dune::GeometryTypes::hexahedron;
            case Dune::VTK::GeometryType::prism: return Dune::GeometryTypes::prism;
            case Dune::VTK::GeometryType::pyramid: return Dune::GeometryTypes::pyramid;
            default: DUNE_THROW(Dune::NotImplemented, "VTK cell type " << vtkCellType);
        }
    }

    template<int dim, std::enable_if_t<(dim < 3), int> = 0>
    std::vector<Dune::FieldVector<double, dim>>
    adaptPointDimension_(std::vector<Dune::FieldVector<double, 3>>&& points3D) const
    {
        std::vector<Dune::FieldVector<double, dim>> points(points3D.size());
        for (std::size_t i = 0; i < points.size(); ++i)
            for (int j = 0; j < dim; ++j)
                points[i][j] = points3D[i][j];
        return points;
    }

    template<int dim, std::enable_if_t<(dim == 3), int> = 0>
    std::vector<Dune::FieldVector<double, dim>>
    adaptPointDimension_(std::vector<Dune::FieldVector<double, 3>>&& points3D) const
    { return std::move(points3D); }

    std::string fileName_; //!< the vtk file name
    tinyxml2::XMLDocument doc_; //!< the xml document created from file with name fileName_
};

} // end namespace Dumux

#endif // DUMUX_HAVE_GRIDFORMAT

#endif
