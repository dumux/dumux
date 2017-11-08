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
 * \brief Provides a grid creator for all supported grid managers with
 *        input file interfaces.
 *
 * \todo add Petrel grids with dune-cornerpoint
 */
#ifndef DUMUX_GRID_CREATOR_HH
#define DUMUX_GRID_CREATOR_HH

#include <array>
#include <bitset>
#include <memory>
#include <sstream>

#include <dune/common/exceptions.hh>
#include <dune/common/classname.hh>
#include <dune/common/parallel/collectivecommunication.hh>
#include <dune/common/parallel/mpihelper.hh>
#include <dune/grid/io/file/dgfparser/dgfparser.hh>
#include <dune/grid/io/file/gmshreader.hh>
#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/utility/structuredgridfactory.hh>

// YaspGrid specific includes
#include <dune/grid/yaspgrid.hh>
#include <dune/grid/io/file/dgfparser/dgfyasp.hh>

 // OneDGrid specific includes
#include <dune/grid/onedgrid.hh>
#include <dune/grid/io/file/dgfparser/dgfoned.hh>

// UGGrid specific includes
#if HAVE_UG
#include <dune/grid/uggrid.hh>
#include <dune/grid/io/file/dgfparser/dgfug.hh>
#endif

// ALUGrid specific includes
#if HAVE_DUNE_ALUGRID
#include <dune/alugrid/grid.hh>
#include <dune/alugrid/dgf.hh>
#endif

// FoamGrid specific includes
#if HAVE_DUNE_FOAMGRID
#include <dune/foamgrid/foamgrid.hh>
#include <dune/common/version.hh>
#if DUNE_VERSION_NEWER(DUNE_COMMON,2,6)
#include <dune/foamgrid/dgffoam.hh>
#else
#include <dune/foamgrid/dgffoam.cc>
#endif
#endif

#include <dumux/common/propertysystem.hh>
#include <dumux/common/parameters.hh>

namespace Dumux
{
namespace Properties
{
// poperty forward declarations
NEW_PROP_TAG(Grid);
NEW_PROP_TAG(GridParameterGroup);
NEW_PROP_TAG(AdaptiveGrid);
NEW_PROP_TAG(Scalar);
NEW_PROP_TAG(ImplicitIsBox);
}


/*!
 * \brief Provides the grid creator base interface (public) and methods common
 *        to most grid creator specializations (protected).
 */
template <class TypeTag, class Grid>
class GridCreatorBase
{
public:

    /*!
     * \brief Returns a reference to the grid.
     */
    static Grid &grid()
    {
        if(enableDgfGridPointer_)
            return *dgfGridPtr();
        else
            return *gridPtr();
    }

    /*!
     * \brief Call the parameters function of the DGF grid pointer if available
     */
    template <class Entity>
    static const std::vector<double>& parameters(const Entity& entity)
    {
        if(enableDgfGridPointer_)
            return dgfGridPtr().parameters(entity);
        else
            DUNE_THROW(Dune::InvalidStateException, "The parameters method is only available if the grid was constructed with a DGF file!");
    }

    /*!
     * \brief Call the parameters function of the DGF grid pointer if available
     */
    template <class GridImp, class IntersectionImp>
    static const Dune::DGFBoundaryParameter::type& parameters(const Dune::Intersection<GridImp, IntersectionImp>& intersection)
    {
        if(enableDgfGridPointer_)
            return dgfGridPtr().parameters(intersection);
        else
            DUNE_THROW(Dune::InvalidStateException, "The parameters method is only available if the grid was constructed with a DGF file!");
    }

    /*!
     * \brief Return the boundary domain marker (Gmsh physical entity number) of an intersection
              Only available when using Gmsh with GridParameterGroup.DomainMarkers = 1.
     * \param boundarySegmentIndex The boundary segment index of the intersection (intersection.boundarySegmentIndex()
     */
    static int getBoundaryDomainMarker(int boundarySegmentIndex)
    {
        if(enableGmshDomainMarkers_)
            return boundaryMarkers_[boundarySegmentIndex];
        else
            DUNE_THROW(Dune::InvalidStateException, "The getBoundaryDomainMarker method is only available if DomainMarkers for Gmsh were enabled!"
                                                     << " If your Gmsh file contains domain markers / physical entities,"
                                                     << " enable them by setting " << GET_PROP_VALUE(TypeTag, GridParameterGroup)
                                                     << ".DomainMarkers = 1 in the input file.");
    }

    /*!
     * \brief Return the element domain marker (Gmsh physical entity number) of an element.
              Only available when using Gmsh with GridParameterGroup.DomainMarkers = 1.
     * \param elementIdx The element index
     */
    static int getElementDomainMarker(int elementIdx)
    {
        if(enableGmshDomainMarkers_)
        {
            if(elementIdx >= grid().levelGridView(0).size(0))
                DUNE_THROW(Dune::RangeError, "Requested element index is bigger than the number of level 0 elements!");
            return elementMarkers_[elementIdx];
        }
        else
            DUNE_THROW(Dune::InvalidStateException, "The getElementDomainMarker method is only available if DomainMarkers for Gmsh were enabled!"
                                                     << " If your Gmsh file contains domain markers / physical entities,"
                                                     << " enable them by setting " << GET_PROP_VALUE(TypeTag, GridParameterGroup)
                                                     << ".DomainMarkers = 1 in the input file.");
    }

    /*!
     * \brief Call loadBalance() function of the grid.
     */
    static void loadBalance()
    {
        if(enableDgfGridPointer_)
            dgfGridPtr().loadBalance();
        else
            gridPtr()->loadBalance();
    }

protected:

    /*!
     * \brief Returns a reference to the grid pointer (std::shared_ptr<Grid>)
     */
    static std::shared_ptr<Grid> &gridPtr()
    {
        if(!enableDgfGridPointer_)
        {
            static std::shared_ptr<Grid> gridPtr_;
            return gridPtr_;
        }
        else
            DUNE_THROW(Dune::InvalidStateException, "You are using DGF. To get the grid pointer use method dgfGridPtr()!");
    }

    /*!
     * \brief Returns a reference to the DGF grid pointer (Dune::GridPtr<Grid>).
     */
    static Dune::GridPtr<Grid> &dgfGridPtr()
    {
        if(enableDgfGridPointer_)
        {
            static Dune::GridPtr<Grid> dgfGridPtr_;
            return dgfGridPtr_;
        }
        else
            DUNE_THROW(Dune::InvalidStateException, "The DGF grid pointer is only available if the grid was constructed with a DGF file!");
    }

    /*!
     * \brief Returns the filename extension of a given filename
     */
    static std::string getFileExtension(const std::string& fileName)
    {
        std::size_t i = fileName.rfind('.', fileName.length());
        if (i != std::string::npos)
        {
            return(fileName.substr(i+1, fileName.length() - i));
        }
        else
        {
            DUNE_THROW(Dune::IOError, "Please provide and extension for your grid file ('"<< fileName << "')!");
        }
        return "";
    }

    /*!
     * \brief Makes a grid from a file. We currently support *.dgf (Dune Grid Format) and *.msh (Gmsh mesh format).
     */
    static void makeGridFromFile(const std::string& fileName,
                                 const std::string& modelParamGroup = "")
    {
        // We found a file in the input file...does it have a supported extension?
        const std::string extension = getFileExtension(fileName);
        if(extension != "dgf" && extension != "msh")
            DUNE_THROW(Dune::IOError, "Grid type " << Dune::className<Grid>() << " only supports DGF (*.dgf) and Gmsh (*.msh) grid files but the specified filename has extension: *."<< extension);

        // make the grid
        if(extension == "dgf")
        {
            enableDgfGridPointer_ = true;
            dgfGridPtr() = Dune::GridPtr<Grid>(fileName.c_str(), Dune::MPIHelper::getCommunicator());
        }
        if(extension == "msh")
        {
            // get some optional parameters
            const bool verbose = getParamFromGroup<bool>(modelParamGroup, "Grid.Verbosity", false);
            const bool boundarySegments = getParamFromGroup<bool>(modelParamGroup, "Grid.BoundarySegments", false);
            const bool domainMarkers = getParamFromGroup<bool>(modelParamGroup, "Grid.DomainMarkers", false);

            if(domainMarkers)
            {
                enableGmshDomainMarkers_ = true;
                if (Dune::MPIHelper::getCollectiveCommunication().rank() == 0)
                {
                    gridPtr() = std::shared_ptr<Grid>(Dune::GmshReader<Grid>::read(fileName, boundaryMarkers_, elementMarkers_, verbose, boundarySegments));
                }
                else
                {
                    Dune::GridFactory<Grid> factory;
                    gridPtr() = std::shared_ptr<Grid>(factory.createGrid());
                }
            }
            else
            {
                if (Dune::MPIHelper::getCollectiveCommunication().rank() == 0)
                {
                    gridPtr() = std::shared_ptr<Grid>(Dune::GmshReader<Grid>::read(fileName, verbose, boundarySegments));
                }
                else
                {
                    Dune::GridFactory<Grid> factory;
                    gridPtr() = std::shared_ptr<Grid>(factory.createGrid());
                }
            }
        }
    }

    /*!
     * \brief Makes a grid from a DGF file. This is used by grid managers that only support DGF.
     */
    static void makeGridFromDgfFile(const std::string& fileName)
    {
        // We found a file in the input file...does it have a supported extension?
        const std::string extension = getFileExtension(fileName);
        if(extension != "dgf")
            DUNE_THROW(Dune::IOError, "Grid type " << Dune::className<Grid>() << " only supports DGF (*.dgf) but the specified filename has extension: *."<< extension);

        enableDgfGridPointer_ = true;
        dgfGridPtr() = Dune::GridPtr<Grid>(fileName.c_str(), Dune::MPIHelper::getCommunicator());
    }

    /*!
     * \brief The cell types for structured grids
     */
    enum CellType {Simplex, Cube};

    /*!
     * \brief Makes a structured cube grid using the structured grid factory
     */
    template <int dim, int dimworld>
    static void makeStructuredGrid(CellType cellType,
                                   const std::string& modelParamGroup = "")
    {
        using GlobalPosition = Dune::FieldVector<typename Grid::ctype, dimworld>;
        const auto upperRight = getParamFromGroup<GlobalPosition>(modelParamGroup, "Grid.UpperRight");
        const auto lowerLeft = getParamFromGroup<GlobalPosition>(modelParamGroup, "Grid.LowerLeft", GlobalPosition(0.0));

        using CellArray = std::array<unsigned int, dim>;
        CellArray cells; cells.fill(1);
        cells = getParamFromGroup<CellArray>(modelParamGroup, "Grid.Cells", cells);

        // make the grid
        if (cellType == CellType::Cube)
        {
            gridPtr() = Dune::StructuredGridFactory<Grid>::createCubeGrid(lowerLeft, upperRight, cells);
        }
        else if (cellType == CellType::Simplex)
        {
            gridPtr() = Dune::StructuredGridFactory<Grid>::createSimplexGrid(lowerLeft, upperRight, cells);
        }
        else
        {
            DUNE_THROW(Dune::GridError, "Unknown cell type for making structured grid! Choose Cube or Simplex.");
        }
    }

    /*!
     * \brief Refines a grid after construction if GridParameterGroup.Refinement is set in the input file
     */
    static void maybeRefineGrid(const std::string& modelParamGroup)
    {
        if (haveParamInGroup(modelParamGroup, "Grid.Refinement"))
            grid().globalRefine(getParamFromGroup<int>(modelParamGroup, "Grid.Refinement"));
    }

    /*!
    * \brief A state variable if the DGF Dune::GridPtr has been enabled.
    *        It is always enabled if a DGF grid file was used to create the grid.
    */
    static bool enableDgfGridPointer_;

    /*!
    * \brief A state variable if domain markers have been read from a Gmsh file.
    */
    static bool enableGmshDomainMarkers_;

    /*!
    * \brief Element and boundary domain markers obtained from Gmsh physical entities
    *        They map from element indices / boundary ids to the physical entity number
    */
    static std::vector<int> elementMarkers_;
    static std::vector<int> boundaryMarkers_;
};

template <class TypeTag, class Grid>
bool GridCreatorBase<TypeTag, Grid>::enableDgfGridPointer_ = false;

template <class TypeTag, class Grid>
bool GridCreatorBase<TypeTag, Grid>::enableGmshDomainMarkers_ = false;

template <class TypeTag, class Grid>
std::vector<int> GridCreatorBase<TypeTag, Grid>::elementMarkers_;

template <class TypeTag, class Grid>
std::vector<int> GridCreatorBase<TypeTag, Grid>::boundaryMarkers_;

/*!
 * \brief Provides the grid creator implementation for all supported grid managers that constructs a grid
 *        from information in the input file. This class is specialised below for all
 *        supported grid managers. It inherits the functionality of the base class.
 */
template <class TypeTag, class Grid>
class GridCreatorImpl : public GridCreatorBase<TypeTag, Grid>
{
public:
    /*!
     * \brief Make the grid. This is implemented by specializations of this class.
     */
    static void makeGrid(const std::string& modelParamGroup = "")
    {
        DUNE_THROW(Dune::NotImplemented,
            "The GridCreator for grid type " << Dune::className<Grid>() << " is not implemented! Consider providing your own GridCreator.");
    }
};

/*!
 * \brief Provides the grid creator (this is the class called by the user) for all supported grid managers that constructs a grid
 *        from information in the input file. This class is specialised below for all
 *        supported grid managers. It inherits the functionality of the base class.
 */
template <class TypeTag>
class GridCreator : public GridCreatorImpl<TypeTag, typename GET_PROP_TYPE(TypeTag, Grid)>
{};

/////////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Specializations //////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////////////

/*!
 * \brief Helper class for determining the default overlap in case of parallel yasp grids
 */
template <class TypeTag>
class YaspOverlapHelper
{
public:
    // trick to set overlap different for implicit models...
    template<class T = TypeTag>
    static typename std::enable_if<Properties::propertyDefined<T, T, PTAG_(ImplicitIsBox)>::value, int>::type
    getOverlap(const std::string& modelParamGroup)
    {
        // the default is dependent on the discretization:
        // our box models only work with overlap 0
        // our cc models only work with overlap > 0
        static const int isBox = GET_PROP_VALUE(TypeTag, ImplicitIsBox);
        int overlap = getParamFromGroup<int>(modelParamGroup, "Grid.Overlap", isBox ? 0 : 1);

        if (isBox && overlap != 0)
            DUNE_THROW(Dune::NotImplemented, "Parallel overlapping grids for box models.");
        if (!isBox && overlap < 1)
            DUNE_THROW(Dune::NotImplemented, "Parallel non-overlapping grids for cc models.");

        return overlap;
    }

    //... then for other models (e.g. sequential)
    template<class T = TypeTag>
    static typename std::enable_if<!Properties::propertyDefined<T, T, PTAG_(ImplicitIsBox)>::value, int>::type
    getOverlap(const std::string& modelParamGroup)
    {
        return getParamFromGroup<int>(modelParamGroup, "Grid.Overlap", 1);
    }

};

/*!
 * \brief Provides a grid creator for YaspGrids
 *        from information in the input file
 *
 * All keys are expected to be in group GridParameterGroup.
 * The following keys are recognized:
 * - File : a DGF file to load the coarse grid from
 * - UpperRight : extension of the domain
 * - Cells : the number of cells in each direction
 * - Periodic : true or false for each direction
 * - Overlap : overlap size in cells
 * - Partitioning : a non-standard load-balancing, number of processors per direction
 * - KeepPyhsicalOverlap : whether to keep the physical overlap
 *     in physical size or in number of cells upon refinement
 * - Refinement : the number of global refines to apply initially.
 *
 */
template<class TypeTag, class ct, int dim>
class GridCreatorImpl<TypeTag, Dune::YaspGrid<dim, Dune::EquidistantCoordinates<ct, dim> > >
          : public GridCreatorBase<TypeTag, Dune::YaspGrid<dim, Dune::EquidistantCoordinates<ct, dim> > >
{
public:
    typedef typename Dune::YaspGrid<dim, Dune::EquidistantCoordinates<ct, dim> > Grid;
    typedef GridCreatorBase<TypeTag, Grid> ParentType;

    /*!
     * \brief Make the grid. This is implemented by specializations of this method.
     */
    static void makeGrid(const std::string& modelParamGroup = "")
    {
        // First try to create it from a DGF file in GridParameterGroup.File
        if (haveParamInGroup(modelParamGroup, "Grid.File"))
        {
            ParentType::makeGridFromDgfFile(getParamFromGroup<std::string>(modelParamGroup, "Grid.File"));
            postProcessing_(modelParamGroup);
            return;
        }

        // Then look for the necessary keys to construct from the input file
        if (!haveParamInGroup(modelParamGroup, "Grid.UpperRight"))
            DUNE_THROW(ParameterException, "Please supply the mandatory parameter "
                                          << modelParamGroup <<  ".UpperRight or a grid file in "
                                          << modelParamGroup << ".File.");

        // get the upper right corner coordinates
        const auto upperRight = getParamFromGroup<Dune::FieldVector<ct, dim>>(modelParamGroup, "Grid.UpperRight");

        // number of cells in each direction
        std::array<int, dim> cells; cells.fill(1);
        cells = getParamFromGroup<std::array<int, dim>>(modelParamGroup, "Grid.Cells", cells);

        // periodic boundaries
        const auto periodic = getParamFromGroup<std::bitset<dim>>(modelParamGroup, "Grid.Periodic", std::bitset<dim>());

        // get the overlap dependent on some template parameters
        const int overlap = YaspOverlapHelper<TypeTag>::getOverlap(modelParamGroup);

        bool default_lb = !haveParamInGroup(modelParamGroup, "Grid.Partitioning");

        // make the grid
        if (default_lb)
        {
            // construct using default load balancing
            ParentType::gridPtr() = std::make_shared<Grid>(upperRight, cells, periodic, overlap);
        }
        else
        {
            // construct using user defined partitioning
            const auto partitioning = getParamFromGroup<std::array<int, dim>>(modelParamGroup, "Grid.Partitioning");
            Dune::YaspFixedSizePartitioner<dim> lb(partitioning);
            ParentType::gridPtr() = std::make_shared<Grid>(upperRight, cells, periodic, overlap, typename Grid::CollectiveCommunicationType(), &lb);
        }
        postProcessing_(modelParamGroup);
    }

private:
    /*!
     * \brief Postprocessing for YaspGrid
     */
    static void postProcessing_(const std::string& modelParamGroup)
    {
        // Check if should refine the grid
        bool keepPhysicalOverlap = getParamFromGroup<bool>(modelParamGroup, "Grid.KeepPhysicalOverlap", true);
        ParentType::grid().refineOptions(keepPhysicalOverlap);
        ParentType::maybeRefineGrid(modelParamGroup);
    }
};

/*!
 * \brief Provides a grid creator for YaspGrids with non-zero offset
 *        from information in the input file
 *
 * All keys are expected to be in group GridParameterGroup.
 * The following keys are recognized:
 * - LowerLeft : lower left corner coordinates
 * - UpperRight : upper right corner coordinates
 * - Cells : the number of cells in each direction
 * - Periodic : true or false for each direction
 * - Overlap : overlap size in cells
 * - Partitioning : a non-standard load-balancing, number of processors per direction
 * - KeepPyhsicalOverlap : whether to keep the physical overlap
 *     in physical size or in number of cells upon refinement
 * - Refinement : the number of global refines to apply initially.
 *
 */
template<class TypeTag, class ct, int dim>
class GridCreatorImpl<TypeTag, Dune::YaspGrid<dim, Dune::EquidistantOffsetCoordinates<ct, dim> > >
          : public GridCreatorBase<TypeTag, Dune::YaspGrid<dim, Dune::EquidistantOffsetCoordinates<ct, dim> > >
{
public:
    typedef typename Dune::YaspGrid<dim, Dune::EquidistantOffsetCoordinates<ct, dim> > Grid;
    typedef GridCreatorBase<TypeTag, Grid> ParentType;

    /*!
     * \brief Make the grid. This is implemented by specializations of this method.
     */
    static void makeGrid(const std::string& modelParamGroup = "")
    {
        // First try to create it from a DGF file in GridParameterGroup.File
        if (haveParamInGroup(modelParamGroup, "Grid.File"))
        {
            ParentType::makeGridFromDgfFile(getParamFromGroup<std::string>(modelParamGroup, "Grid.File"));
            postProcessing_(modelParamGroup);
            return;
        }

        // Then look for the necessary keys to construct from the input file
        if (!haveParamInGroup(modelParamGroup, "Grid.UpperRight"))
            DUNE_THROW(ParameterException, "Please supply the mandatory parameter "
                                          << modelParamGroup <<  ".UpperRight or a grid file in "
                                          << modelParamGroup << ".File.");

        using GlobalPosition = Dune::FieldVector<ct, dim>;

        // get the upper right corner coordinates
        const auto upperRight = getParamFromGroup<GlobalPosition>(modelParamGroup, "Grid.UpperRight");

        // get the lowerLeft corner coordinates (have a default value)
        const auto lowerLeft = getParamFromGroup<GlobalPosition>(modelParamGroup, "Grid.LowerLeft", GlobalPosition(0.0));

        // number of cells in each direction
        std::array<int, dim> cells; cells.fill(1);
        cells = getParamFromGroup<std::array<int, dim>>(modelParamGroup, "Grid.Cells", cells);

        // periodic boundaries
        const auto periodic = getParamFromGroup<std::bitset<dim>>(modelParamGroup, "Grid.Periodic", std::bitset<dim>());

        // get the overlap dependent on some template parameters
        const int overlap = YaspOverlapHelper<TypeTag>::getOverlap(modelParamGroup);

        bool default_lb = !haveParamInGroup(modelParamGroup, "Grid.Partitioning");

        // make the grid
        if (default_lb)
        {
            // construct using default load balancing
            ParentType::gridPtr() = std::make_shared<Grid>(lowerLeft, upperRight, cells, periodic, overlap);
        }
        else
        {
            // construct using user defined partitioning
            const auto partitioning = getParamFromGroup<std::array<int, dim>>(modelParamGroup, "Grid.Partitioning");
            Dune::YaspFixedSizePartitioner<dim> lb(partitioning);
            ParentType::gridPtr() = std::make_shared<Grid>(lowerLeft, upperRight, cells, periodic, overlap, typename Grid::CollectiveCommunicationType(), &lb);
        }
        postProcessing_(modelParamGroup);
    }

private:
    /*!
     * \brief Postprocessing for YaspGrid
     */
    static void postProcessing_(const std::string& modelParamGroup)
    {
        // Check if should refine the grid
        const bool keepPhysicalOverlap = getParamFromGroup<bool>(modelParamGroup, "Grid.KeepPhysicalOverlap", true);
        ParentType::grid().refineOptions(keepPhysicalOverlap);
        ParentType::maybeRefineGrid(modelParamGroup);
    }
};

/*!
 * \brief Provides a grid creator for YaspGrids with different zones and grading
 *
 * All keys are expected to be in group GridParameterGroup.
 * The following keys are recognized:
 * - Positions0 : position array for x-coordinate
 * - Positions1 : position array for y-coordinate
 * - Positions2 : position array for z-coordinate
 * - Cells0 : number of cells array for x-coordinate
 * - Cells1 : number of cells array for y-coordinate
 * - Cells2 : number of cells array for z-coordinate
 * - Grading0 : grading factor array for x-coordinate
 * - Grading1 : grading factor array for y-coordinate
 * - Grading2 : grading factor array for z-coordinate
 * - Verbosity : whether the grid construction should output to standard out
 * - Periodic : true or false for each direction
 * - Overlap : overlap size in cells
 * - Partitioning : a non-standard load-balancing, number of processors per direction
 * - KeepPyhsicalOverlap : whether to keep the physical overlap
 *     in physical size or in number of cells upon refinement
 * - Refinement : the number of global refines to apply initially.
 *
 * The grading factor \f$ g \f$ specifies the ratio between the next and the current cell size:
 * \f$ g = \frac{h_{i+1}}{h_i} \f$.
 * Negative grading factors are converted to
 * \f$ g = -\frac{1}{g_\textrm{negative}} \f$
 * to avoid issues with imprecise fraction numbers.
 */
template<class TypeTag, class ct, int dim>
class GridCreatorImpl<TypeTag, Dune::YaspGrid<dim, Dune::TensorProductCoordinates<ct, dim> > >
          : public GridCreatorBase<TypeTag, Dune::YaspGrid<dim, Dune::TensorProductCoordinates<ct, dim> > >
{
public:
    typedef typename Dune::YaspGrid<dim, Dune::TensorProductCoordinates<ct, dim> > Grid;
    typedef GridCreatorBase<TypeTag, Grid> ParentType;

    /*!
     * \brief Make the grid. This is implemented by specializations of this method.
     */
    static void makeGrid()
    {
        // Only construction from the input file is possible
        // Look for the necessary keys to construct from the input file
        try {
            typedef ct Scalar;
            typedef std::vector<int> IntVector;
            typedef std::vector<Scalar> ScalarVector;

            // The positions
            std::array<ScalarVector, dim> positions;
            for (int i = 0; i < dim; ++i)
            {
              std::string paramName = GET_PROP_VALUE(TypeTag, GridParameterGroup);
              paramName += ".Positions";
              paramName += std::to_string(i);
              positions[i] = GET_RUNTIME_PARAM_CSTRING(TypeTag, ScalarVector, paramName.c_str());
            }

            // the number of cells (has a default)
            std::array<IntVector, dim> cells;
            for (int i = 0; i < dim; ++i)
            {
              std::string paramName = GET_PROP_VALUE(TypeTag, GridParameterGroup);
              paramName += ".Cells";
              paramName += std::to_string(i);
              try { cells[i] = GET_RUNTIME_PARAM_CSTRING(TypeTag, IntVector, paramName.c_str()); }
              catch (ParameterException &e) { cells[i].resize(positions[i].size()-1, 1.0); }
            }

            // grading factor (has a default)
            std::array<ScalarVector, dim> grading;
            for (int i = 0; i < dim; ++i)
            {
              std::string paramName = GET_PROP_VALUE(TypeTag, GridParameterGroup);
              paramName += ".Grading";
              paramName += std::to_string(i);
              try { grading[i] = GET_RUNTIME_PARAM_CSTRING(TypeTag, ScalarVector, paramName.c_str()); }
              catch (ParameterException &e) { grading[i].resize(positions[i].size()-1, 1.0); }
            }

            // The optional parameters (they have a default)
            typedef std::bitset<dim> BitSet;
            BitSet periodic;
            try { periodic = GET_RUNTIME_PARAM_FROM_GROUP_CSTRING(TypeTag, BitSet, GET_PROP_VALUE(TypeTag, GridParameterGroup).c_str(), Periodic);}
            catch (ParameterException &e) { }

            // the default is dependent on the discretization
            int overlap = YaspOverlapHelper<TypeTag>::getOverlap();

            bool default_lb = false;
            typedef std::array<int, dim> CellArray;
            CellArray partitioning;
            try { partitioning = GET_RUNTIME_PARAM_FROM_GROUP_CSTRING(TypeTag, CellArray, GET_PROP_VALUE(TypeTag, GridParameterGroup).c_str(), Partitioning);}
            catch (ParameterException &e) { default_lb = true; }

            //make the grid
            // sanity check of the input parameters
            bool verbose = false;
            try { verbose = GET_RUNTIME_PARAM_FROM_GROUP_CSTRING(TypeTag, bool, GET_PROP_VALUE(TypeTag, GridParameterGroup).c_str(), Verbosity);}
            catch (ParameterException &e) {  }

            for (unsigned int dimIdx = 0; dimIdx < dim; ++dimIdx)
            {
                if (cells[dimIdx].size() + 1 != positions[dimIdx].size())
                {
                    DUNE_THROW(Dune::RangeError, "Make sure to specify correct \"Cells\" and \"Positions\" arrays");
                }
                if (grading[dimIdx].size() + 1 != positions[dimIdx].size())
                {
                    DUNE_THROW(Dune::RangeError, "Make sure to specify correct \"Grading\" and \"Positions\" arrays");
                }
                Scalar temp = std::numeric_limits<Scalar>::lowest();
                for (unsigned int posIdx = 0; posIdx < positions[dimIdx].size(); ++posIdx)
                {
                    if (temp > positions[dimIdx][posIdx])
                    {
                        DUNE_THROW(Dune::RangeError, "Make sure to specify a monotone increasing \"Positions\" array");
                    }
                    temp = positions[dimIdx][posIdx];
                }
            }

            std::array<ScalarVector, dim> globalPositions;
            using std::pow;
            for (int dimIdx = 0; dimIdx < dim; dimIdx++)
            {
                for (int zoneIdx = 0; zoneIdx < cells[dimIdx].size(); ++zoneIdx)
                {
                    Scalar lower = positions[dimIdx][zoneIdx];
                    Scalar upper = positions[dimIdx][zoneIdx+1];
                    int numCells = cells[dimIdx][zoneIdx];
                    Scalar gradingFactor = grading[dimIdx][zoneIdx];
                    Scalar length = upper - lower;
                    Scalar height = 1.0;
                    bool increasingCellSize = false;

                    if (verbose)
                    {
                        std::cout << "dim " << dimIdx
                                  << " lower "  << lower
                                  << " upper "  << upper
                                  << " numCells "  << numCells
                                  << " grading "  << gradingFactor;
                    }

                    if (gradingFactor > 1.0)
                    {
                        increasingCellSize = true;
                    }

                    // take absolute values and reverse cell size increment to achieve
                    // reverse behavior for negative values
                    if (gradingFactor < 0.0)
                    {
                        using std::abs;
                        gradingFactor = abs(gradingFactor);
                        if (gradingFactor < 1.0)
                        {
                            increasingCellSize = true;
                        }
                    }

                    // if the grading factor is exactly 1.0 do equal spacing
                    if (gradingFactor > 1.0 - 1e-7 && gradingFactor < 1.0 + 1e-7)
                    {
                        height = 1.0 / numCells;
                        if (verbose)
                        {
                            std::cout << " -> h "  << height * length << std::endl;
                        }
                    }
                    // if grading factor is not 1.0, do power law spacing
                    else
                    {
                        height = (1.0 - gradingFactor) / (1.0 - pow(gradingFactor, numCells));

                        if (verbose)
                        {
                            std::cout << " -> grading_eff "  << gradingFactor
                                      << " h_min "  << height * pow(gradingFactor, 0) * length
                                      << " h_max "  << height * pow(gradingFactor, numCells-1) * length
                                      << std::endl;
                        }
                    }

                    std::vector<Scalar> localPositions;
                    localPositions.push_back(0);
                    for (int i = 0; i < numCells-1; i++)
                    {
                        Scalar hI = height;
                        if (!(gradingFactor < 1.0 + 1e-7 && gradingFactor > 1.0 - 1e-7))
                        {
                            if (increasingCellSize)
                            {
                                hI *= pow(gradingFactor, i);
                            }
                            else
                            {
                                hI *= pow(gradingFactor, numCells-i-1);
                            }
                        }
                        localPositions.push_back(localPositions[i] + hI);
                    }

                    for (int i = 0; i < localPositions.size(); i++)
                    {
                        localPositions[i] *= length;
                        localPositions[i] += lower;
                    }


                    for (unsigned int i = 0; i < localPositions.size(); ++i)
                    {
                        globalPositions[dimIdx].push_back(localPositions[i]);
                    }
                }
                globalPositions[dimIdx].push_back(positions[dimIdx].back());
            }

            if (default_lb)
                ParentType::gridPtr() = std::make_shared<Grid>(globalPositions, periodic, overlap);
            else
            {
                typename Dune::YaspFixedSizePartitioner<dim> lb(partitioning);
                ParentType::gridPtr() = std::make_shared<Grid>(globalPositions, periodic, overlap, typename Grid::CollectiveCommunicationType(), &lb);
            }
            postProcessing_();
        }
        catch (ParameterException &e) {
                DUNE_THROW(ParameterException, "Please supply the mandatory parameters:" << std::endl
                                              << GET_PROP_VALUE(TypeTag, GridParameterGroup) << ".Positions0, ..." << std::endl);
        }
        catch (...) { throw; }
    }

private:
    /*!
     * \brief Postprocessing for YaspGrid
     */
    static void postProcessing_()
    {
        // Check if should refine the grid
        bool keepPhysicalOverlap = true;
        try { keepPhysicalOverlap = GET_RUNTIME_PARAM_FROM_GROUP_CSTRING(TypeTag, bool, GET_PROP_VALUE(TypeTag, GridParameterGroup).c_str(), KeepPhysicalOverlap);}
        catch (ParameterException &e) { }
        ParentType::grid().refineOptions(keepPhysicalOverlap);
        ParentType::maybeRefineGrid();
    }
};

/*!
 * \brief Provides a grid creator for OneDGrids
 *        from information in the input file
 *
 * All keys are expected to be in group GridParameterGroup.
 * The following keys are recognized:
 * - LeftBoundary : start coordinate
 * - RightBoundary : end coordinate
 * - Cells : the number of cell
 * - RefinementType : local or copy
 * - Refinement : the number of global refines to apply initially.
 *
 */
template<class TypeTag>
class GridCreatorImpl<TypeTag, Dune::OneDGrid>
          : public GridCreatorBase<TypeTag, Dune::OneDGrid>
{
public:
    typedef typename Dune::OneDGrid Grid;
    typedef GridCreatorBase<TypeTag, Grid> ParentType;

    /*!
     * \brief Make the grid. This is implemented by specializations of this method.
     */
    static void makeGrid()
    {
        // First try to create it from a DGF file in GridParameterGroup.File
        try {
            const std::string fileName = GET_RUNTIME_PARAM_FROM_GROUP_CSTRING(TypeTag, std::string, GET_PROP_VALUE(TypeTag, GridParameterGroup).c_str(), File);
            ParentType::makeGridFromDgfFile(fileName);
            postProcessing_();
            return;
        }
        catch (ParameterException &e) {}
        catch (...) { throw; }

        // Look for the necessary keys to construct from the input file
        try {
            // The required parameters
            typedef typename Grid::ctype Coordinate;
            const Coordinate leftBoundary = GET_RUNTIME_PARAM_FROM_GROUP_CSTRING(TypeTag, Coordinate, GET_PROP_VALUE(TypeTag, GridParameterGroup).c_str(), LeftBoundary);
            const Coordinate rightBoundary = GET_RUNTIME_PARAM_FROM_GROUP_CSTRING(TypeTag, Coordinate, GET_PROP_VALUE(TypeTag, GridParameterGroup).c_str(), RightBoundary);

            // The optional parameters
            int cells = 1;
            try { cells = GET_RUNTIME_PARAM_FROM_GROUP_CSTRING(TypeTag, int, GET_PROP_VALUE(TypeTag, GridParameterGroup).c_str(), Cells);}
            catch (ParameterException &e) { }

            ParentType::gridPtr() = std::make_shared<Grid>(cells, leftBoundary, rightBoundary);
            postProcessing_();
        }
        catch (ParameterException &e) {}
        catch (...) { throw; }

        // Look for the necessary keys to construct from the input file with just a coordinates vector
        try {
            // The required parameters
            typedef std::vector<typename Grid::ctype> Coordinates;
            const Coordinates coordinates = GET_RUNTIME_PARAM_FROM_GROUP_CSTRING(TypeTag, Coordinates, GET_PROP_VALUE(TypeTag, GridParameterGroup).c_str(), Coordinates);
            // make the grid
            ParentType::gridPtr() = std::make_shared<Grid>(coordinates);
            postProcessing_();
        }
        catch (ParameterException &e) {
            DUNE_THROW(ParameterException, "Please supply the mandatory parameters "
                                         << GET_PROP_VALUE(TypeTag, GridParameterGroup) << ".LeftBoundary, "
                                         << GET_PROP_VALUE(TypeTag, GridParameterGroup) << ".RightBoundary or "
                                         << GET_PROP_VALUE(TypeTag, GridParameterGroup) << ".Coordinates or a grid file in "
                                         << GET_PROP_VALUE(TypeTag, GridParameterGroup) << ".File.");
        }
        catch (...) { throw; }
    }

private:
    /*!
     * \brief Do some operatrion after making the grid, like global refinement
     */
    static void postProcessing_()
    {
        // Check for refinement type
        std::string refType = "Local";
        try { refType = GET_RUNTIME_PARAM_FROM_GROUP_CSTRING(TypeTag, std::string, GET_PROP_VALUE(TypeTag, GridParameterGroup).c_str(), RefinementType);
              if (refType != "Local" && refType != "Copy")
                   DUNE_THROW(Dune::IOError, "OneGrid only supports 'Local' or 'Copy' as refinment type. Not '"<< refType<<"'!");
        }
        catch (ParameterException &e) {}
        catch (...) { throw; }

        if (refType == "Local")
            ParentType::grid().setRefinementType(Dune::OneDGrid::RefinementType::LOCAL);
        if (refType == "Copy")
            ParentType::grid().setRefinementType(Dune::OneDGrid::RefinementType::COPY);

        // Check if should refine the grid
        ParentType::maybeRefineGrid();
    }
};

#if HAVE_UG

/*!
 * \brief Provides a grid creator for UGGrids
 *        from information in the input file
 *
 * All keys are expected to be in group GridParameterGroup.

 * The following keys are recognized:
 * - File : A DGF or gmsh file to load from, type detection by file extension
 * - LowerLeft : lowerLeft corner of a structured grid
 * - UpperRight : upperright corner of a structured grid
 * - Cells : number of elements in a structured grid
 * - CellType : "Cube" or "Simplex" to be used for structured grids
 * - Refinement : the number of global refines to perform
 * - Verbosity : whether the grid construction should output to standard out
 * - HeapSize: The heapsize used to allocate memory
 * - BoundarySegments : whether to insert boundary segments into the grid
 *
 */
template<class TypeTag, int dim>
class GridCreatorImpl<TypeTag, Dune::UGGrid<dim> >
          : public GridCreatorBase<TypeTag, Dune::UGGrid<dim> >
{
public:
    typedef typename Dune::UGGrid<dim> Grid;
    typedef GridCreatorBase<TypeTag, Grid> ParentType;

    /*!
     * \brief Make the UGGrid.
     */
    static void makeGrid()
    {
        // First try to create it from a DGF or msh file in GridParameterGroup.File
        try {
            const std::string fileName = GET_RUNTIME_PARAM_FROM_GROUP_CSTRING(TypeTag, std::string, GET_PROP_VALUE(TypeTag, GridParameterGroup).c_str(), File);
            preProcessing_();
            ParentType::makeGridFromFile(fileName);
            postProcessing_();
            return;
        }
        catch (ParameterException &e) {}
        catch (...) { throw; }

        // Then look for the necessary keys to construct from the input file
        try {
            preProcessing_();
            // Check for cell type
            std::string cellType = "Cube";
            try { cellType = GET_RUNTIME_PARAM_FROM_GROUP_CSTRING(TypeTag, std::string, GET_PROP_VALUE(TypeTag, GridParameterGroup).c_str(), CellType);
                if (cellType != "Cube" && cellType != "Simplex")
                    DUNE_THROW(Dune::IOError, "UGGrid only supports 'Cube' or 'Simplex' as cell type. Not '"<< cellType<<"'!");
            }
            catch (ParameterException &e) {}
            catch (...) { throw; }

            // make the grid
            if (cellType == "Cube")
                ParentType::template makeStructuredGrid<dim, dim>(ParentType::CellType::Cube);
            if (cellType == "Simplex")
                ParentType::template makeStructuredGrid<dim, dim>(ParentType::CellType::Simplex);
            postProcessing_();
        }
        catch (ParameterException &e) {
                DUNE_THROW(ParameterException, "Please supply the mandatory parameters "
                                              << GET_PROP_VALUE(TypeTag, GridParameterGroup) << ".UpperRight or a grid file in "
                                              << GET_PROP_VALUE(TypeTag, GridParameterGroup) << ".File.");
        }
        catch (...) { throw; }
    }

private:
    /*!
     * \brief Do some operatrion before making the grid
     */
    static void preProcessing_()
    {
        bool setDefaultHeapSize = true;
        unsigned defaultHeapSize;
        try { defaultHeapSize = GET_RUNTIME_PARAM_FROM_GROUP_CSTRING(TypeTag, unsigned, GET_PROP_VALUE(TypeTag, GridParameterGroup).c_str(), HeapSize);}
        catch (ParameterException &e) { setDefaultHeapSize = false; }

        if(setDefaultHeapSize)
            Grid::setDefaultHeapSize(defaultHeapSize);
    }

    /*!
     * \brief Do some operatrion after making the grid, like global refinement
     */
    static void postProcessing_()
    {
        // Check for refinement type
        std::string refType = "Local";
        try { refType = GET_RUNTIME_PARAM_FROM_GROUP_CSTRING(TypeTag, std::string, GET_PROP_VALUE(TypeTag, GridParameterGroup).c_str(), RefinementType);
              if (refType != "Local" && refType != "Copy")
                   DUNE_THROW(Dune::IOError, "UGGrid only supports 'Local' or 'Copy' as refinment type. Not '"<< refType<<"'!");
        }
        catch (ParameterException &e) {}
        catch (...) { throw; }

        if (refType == "Local")
            ParentType::grid().setRefinementType(Dune::UGGrid<dim>::RefinementType::LOCAL);
        if (refType == "Copy")
            ParentType::grid().setRefinementType(Dune::UGGrid<dim>::RefinementType::COPY);

        // Check for closure type
        std::string closureType = "Green";
        try { closureType = GET_RUNTIME_PARAM_FROM_GROUP_CSTRING(TypeTag, std::string, GET_PROP_VALUE(TypeTag, GridParameterGroup).c_str(), ClosureType);
              if (closureType != "None" && closureType != "Green")
                   DUNE_THROW(Dune::IOError, "UGGrid only supports 'Green' or 'None' as closure type. Not '"<< closureType<<"'!");
        }
        catch (ParameterException &e) {}
        catch (...) { throw; }

        if (closureType == "Green")
            ParentType::grid().setClosureType(Dune::UGGrid<dim>::ClosureType::GREEN);
        if (closureType == "None")
            ParentType::grid().setClosureType(Dune::UGGrid<dim>::ClosureType::NONE);

        // Check if should refine the grid
        ParentType::maybeRefineGrid();
    }
};

#endif // HAVE_UG

#if HAVE_DUNE_ALUGRID

/*!
 * \brief Provides a grid creator for Dune ALUGrids
 *        from information in the input file
 *
 * All keys are expected to be in group GridParameterGroup.

 * The following keys are recognized:
 * - File : A DGF or gmsh file to load from, type detection by file extension
 * - LowerLeft : lowerLeft corner of a structured grid
 * - UpperRight : upperright corner of a structured grid
 * - Cells : number of elements in a structured grid
 * - Refinement : the number of global refines to perform
 * - Verbosity : whether the grid construction should output to standard out
 * - BoundarySegments : whether to insert boundary segments into the grid
 *
 */
template<class TypeTag, int dim, int dimworld, Dune::ALUGridElementType elType, Dune::ALUGridRefinementType refinementType>
class GridCreatorImpl<TypeTag, Dune::ALUGrid<dim, dimworld, elType, refinementType> >
          : public GridCreatorBase<TypeTag, Dune::ALUGrid<dim, dimworld, elType, refinementType> >
{
public:
    typedef typename Dune::ALUGrid<dim, dimworld, elType, refinementType> Grid;
    typedef GridCreatorBase<TypeTag, Grid> ParentType;

    /*!
     * \brief Make the grid. This is implemented by specializations of this method.
     */
    static void makeGrid()
    {
#if HAVE_DUNE_ALUGRID
        // First check if a restart for an adaptive grid is required
        try {
            if (GET_PROP_VALUE(TypeTag, AdaptiveGrid))
            {
                typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
                Scalar restartTime = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, TimeManager, Restart);
                // we came until here so the restart key was found. Restore the grid.
                int rank = 0;
#if HAVE_MPI
                MPI_Comm_rank(Dune::MPIHelper::getCommunicator(), &rank);
#endif
                try {
                    std::string name = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Problem, Name);
                    std::ostringstream oss;
                    oss << name << "_time=" << restartTime << "_rank=" << rank << ".grs";
                    std::cout << "Restoring an ALUGrid from " << oss.str() << std::endl;
                    ParentType::gridPtr() = std::shared_ptr<Grid>(Dune::BackupRestoreFacility<Grid>::restore(oss.str()));
                    return;
                }
                catch (ParameterException &e)
                {
                    std::cerr << e.what() << std::endl;
                    std::cerr << "Restart functionality for an adaptive grid requested, but failed." << std::endl;
                    std::cerr << "Did you forget to provide Problem.Name in your .input file?" << std::endl;
                    throw;
                }
            }
        }
        catch (ParameterException &e) {}
        catch (...) { throw; }
#endif

        // Then try to create it from a DGF or msh file in GridParameterGroup.File
        try {
            const std::string fileName = GET_RUNTIME_PARAM_FROM_GROUP_CSTRING(TypeTag, std::string, GET_PROP_VALUE(TypeTag, GridParameterGroup).c_str(), File);
            ParentType::makeGridFromFile(fileName);
            ParentType::maybeRefineGrid();
            return;
        }
        catch (ParameterException &e) {}
        catch (...) { throw; }

        // Then look for the necessary keys to construct from the input file
        try {
            // make a structured grid
            if (elType == Dune::cube)
                ParentType::template makeStructuredGrid<dim, dimworld>(ParentType::CellType::Cube);
            if (elType == Dune::simplex)
                ParentType::template makeStructuredGrid<dim, dimworld>(ParentType::CellType::Simplex);
            ParentType::maybeRefineGrid();
        }
        catch (ParameterException &e) {
                DUNE_THROW(ParameterException, "Please supply the mandatory parameters "
                                              << GET_PROP_VALUE(TypeTag, GridParameterGroup) << ".UpperRight or a grid file in "
                                              << GET_PROP_VALUE(TypeTag, GridParameterGroup) << ".File.");
        }
        catch (...) { throw; }
    }
};

#endif // HAVE_DUNE_ALUGRID

#if HAVE_DUNE_FOAMGRID

/*!
 * \brief Provides a grid creator for FoamGrids
 *        from information in the input file
 *
 * All keys are expected to be in group GridParameterGroup.

 * The following keys are recognized:
 * - File : A DGF or gmsh file to load from, type detection by file extension
 * - Verbosity : whether the grid construction should output to standard out
 * - LowerLeft : lowerLeft corner of a structured grid
 * - UpperRight : upperright corner of a structured grid
 * - Cells : number of elements in a structured grid
 *
 */
template<class TypeTag, int dim, int dimworld>
class GridCreatorImpl<TypeTag, Dune::FoamGrid<dim, dimworld> >
          : public GridCreatorBase<TypeTag, Dune::FoamGrid<dim, dimworld> >
{
public:
    typedef typename Dune::FoamGrid<dim, dimworld> Grid;
    typedef GridCreatorBase<TypeTag, Grid> ParentType;

    /*!
     * \brief Make the grid. This is implemented by specializations of this method.
     */
    static void makeGrid(const std::string& modelParamGroup = "")
    {
        // try to create it from file
        if (haveParamInGroup(modelParamGroup, "Grid.File"))
        {
            ParentType::makeGridFromFile(getParamFromGroup<std::string>(modelParamGroup, "Grid.File"), modelParamGroup);
            ParentType::maybeRefineGrid(modelParamGroup);
            return;
        }

        // Then look for the necessary keys to construct a structured grid from the input file
        try {
            ParentType::template makeStructuredGrid<dim, dimworld>(ParentType::CellType::Simplex, modelParamGroup);
            ParentType::maybeRefineGrid(modelParamGroup);
        }
        catch (ParameterException &e) {
                DUNE_THROW(ParameterException, "Please supply the mandatory parameters "
                                              << "Grid.UpperRight or a grid file in Grid.File.");
        }
        catch (...) { throw; }
    }
};

/*!
 * \brief Provides a grid creator for FoamGrids of dim 1
 *        from information in the input file
 *
 * All keys are expected to be in group GridParameterGroup.

 * The following keys are recognized:
 * - File : A DGF or gmsh file to load from, type detection by file extension
 * - Verbosity : whether the grid construction should output to standard out
 * - LowerLeft : lowerLeft corner of a structured grid
 * - UpperRight : upperright corner of a structured grid
 * - Cells : number of elements in a structured grid
 *
 */
template<class TypeTag, int dimworld>
class GridCreatorImpl<TypeTag, Dune::FoamGrid<1, dimworld> >
          : public GridCreatorBase<TypeTag, Dune::FoamGrid<1, dimworld> >
{
public:
    typedef typename Dune::FoamGrid<1, dimworld> Grid;
    typedef GridCreatorBase<TypeTag, Grid> ParentType;

    /*!
     * \brief Make the grid. This is implemented by specializations of this method.
     */
    static void makeGrid(const std::string& modelParamGroup = "")
    {
        // try to create it from file
        if (haveParamInGroup(modelParamGroup, "Grid.File"))
        {
            ParentType::makeGridFromFile(getParamFromGroup<std::string>(modelParamGroup, "Grid.File"), modelParamGroup);
            ParentType::maybeRefineGrid(modelParamGroup);
            return;
        }

        // The required parameters
        using GlobalPosition = Dune::FieldVector<typename Grid::ctype, dimworld>;
        const auto upperRight = getParamFromGroup<GlobalPosition>(modelParamGroup, "Grid.UpperRight");
        const auto lowerLeft = getParamFromGroup<GlobalPosition>(modelParamGroup, "Grid.UpperRight", GlobalPosition(0.0));
        using CellArray = std::array<unsigned int, 1>;
        const auto cells = getParamFromGroup<CellArray>(modelParamGroup, "Grid.Cells", std::array<unsigned int, 1>{{1}});

        // make the grid (structured interval grid in dimworld space)
        Dune::GridFactory<Grid> factory;
        constexpr auto geomType = Dune::GeometryTypes::simplex(1);

        // create a step vector
        GlobalPosition step = upperRight;
        step -= lowerLeft, step /= cells[0];

        // create the vertices
        GlobalPosition globalPos = lowerLeft;
        for (unsigned int vIdx = 0; vIdx <= cells[0]; vIdx++, globalPos += step)
            factory.insertVertex(globalPos);

        // create the cells
        for(unsigned int eIdx = 0; eIdx < cells[0]; eIdx++)
            factory.insertElement(geomType, {eIdx, eIdx+1});

        ParentType::gridPtr() = std::shared_ptr<Grid>(factory.createGrid());
        ParentType::maybeRefineGrid(modelParamGroup);
    }
};

#endif // HAVE_DUNE_FOAMGRID

} // namespace Dumux

#endif
