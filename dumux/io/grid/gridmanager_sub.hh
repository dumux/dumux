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
 * \brief Grid manager specialization for SubGrid
 */
#ifndef DUMUX_IO_GRID_MANAGER_SUB_HH
#define DUMUX_IO_GRID_MANAGER_SUB_HH

#include <memory>

#include <dune/common/deprecated.hh>
#include <dune/common/shared_ptr.hh>
#include <dune/common/concept.hh>

// SubGrid specific includes
#if HAVE_DUNE_SUBGRID
#include <dune/subgrid/subgrid.hh>
#include <dumux/io/rasterimagereader.hh>
// TODO: deprecated: remove this after 3.1 is released
#include <dune/grid/io/file/vtk.hh>
// TODO: deprecated: remove this after 3.1 is released
#include <dune/grid/io/file/dgfparser/dgfwriter.hh>
#endif

#ifndef DUMUX_IO_GRID_MANAGER_HH
#include <dumux/io/grid/gridmanager.hh>
#endif

#include <dumux/common/parameters.hh>
#include <dumux/common/boundaryflag.hh>

#if HAVE_DUNE_SUBGRID
namespace Dumux {
namespace Concept {

/*!
 * \ingroup InputOutput
 * \brief The element selector concept
 */
template<class Element>
struct ElementSelector
{
    template<class F>
    auto require(F&& f) -> decltype(
        bool(f(std::declval<const Element&>()))
    );
};
} // end namespace Concept

/*!
 * \ingroup InputOutput
 * \brief The base class for grid managers for dune-subgrid.
 */
template <class HostGrid, class HostGridManager = GridManager<HostGrid>>
class SubGridManagerBase
: public GridManagerBase<Dune::SubGrid<HostGrid::dimension, HostGrid>>
{
    static constexpr int dim = HostGrid::dimension;
    using HostElement = typename HostGrid::template Codim<0>::Entity;

public:
    using Grid = Dune::SubGrid<dim, HostGrid>;

    /*!
     * \brief Make the grid using an externally created host grid.
     */
    template<class ES,
             typename std::enable_if_t<Dune::models<Concept::ElementSelector<HostElement>, ES>(), int> = 0>
    void init(HostGrid& hostGrid,
              const ES& selector,
              const std::string& paramGroup = "")
    {
        this->gridPtr() = createGrid_(hostGrid, selector, paramGroup);
        loadBalance();
    }

    /*!
     * \brief Make the grid and create the host grid internally.
     */
    template<class ES,
             typename std::enable_if_t<Dune::models<Concept::ElementSelector<HostElement>, ES>(), int> = 0>
    void init(const ES& selector,
              const std::string& paramGroup = "")
    {
        initHostGrid_(paramGroup);
        this->gridPtr() = createGrid_(hostGridManager_->grid(), selector, paramGroup);
        loadBalance();
    }

    /*!
     * \brief Make the subgrid.
     */
    template<class ElementSelector>
    DUNE_DEPRECATED_MSG("Create an instance of this class and use subgridManager.init(hostGrid, selector, paramGroup)")
    static std::unique_ptr<Grid> makeGrid(HostGrid& hostGrid,
                                          const ElementSelector& selector,
                                          const std::string& paramGroup = "")
    {
        std::cerr << "Deprecation warning: SubGridManager::makeGrid is deprecated."
                  << "Create an instance of this class and use subgridManager.init(hostGrid, selector, paramGroup)." << std::endl;

        auto subgridPtr = createGrid_(hostGrid, selector, paramGroup);

        // TODO: remove this after 3.1 is released
        // If desired, write out the final subgrid as a dgf file.
        if (getParamFromGroup<bool>(paramGroup, "Grid.WriteSubGridToDGF", false))
        {
            std::cerr << "Deprecation warning: SubGridManager: Grid.WriteSubGridToDGF is deprecated."
                      << "Use Dune::VTKWriter to write out your grid manually." << std::endl;

            const auto postfix = getParamFromGroup<std::string>(paramGroup, "Problem.Name", "");
            const std::string name = postfix == "" ? "subgrid" : "subgrid_" + postfix;
            Dune::DGFWriter<typename Grid::LeafGridView> writer(subgridPtr->leafGridView());
            writer.write(name + ".dgf");
        }

        // TODO: remove this after 3.1 is released
        // If desired, write out the hostgrid as vtk file.
        if (getParamFromGroup<bool>(paramGroup, "Grid.WriteSubGridToVtk", false))
        {
            std::cerr << "Deprecation warning: SubGridManager: Grid.WriteSubGridToVtk is deprecated."
                      << "Use Dune::VTKWriter to write out your grid manually." << std::endl;

            const auto postfix = getParamFromGroup<std::string>(paramGroup, "Problem.Name", "");
            const std::string name = postfix == "" ? "subgrid" : "subgrid_" + postfix;
            Dune::VTKWriter<typename Grid::LeafGridView> vtkWriter(subgridPtr->leafGridView());
            vtkWriter.write(name);
        }

        // Return a unique pointer to the subgrid.
        return subgridPtr;
    }

    /*!
     * \brief Call loadBalance() function of the grid.
     */
    void loadBalance()
    {
        if (Dune::MPIHelper::getCollectiveCommunication().size() > 1)
            this->grid().loadBalance();
    }

protected:

    /*!
     * \brief Make the subgrid.
     */
    template<class ES,
             typename std::enable_if_t<Dune::models<Concept::ElementSelector<HostElement>, ES>(), int> = 0>
    static std::unique_ptr<Grid> createGrid_(HostGrid& hostGrid,
                                             const ES& selector,
                                             const std::string& paramGroup = "")
    {
        // A unique pointer to the subgrid.
        auto subgridPtr = std::make_unique<Grid>(hostGrid);

        // A container to store the host grid elements' ids.
        std::set<typename HostGrid::Traits::GlobalIdSet::IdType> elementsForSubgrid;
        const auto& globalIDset = subgridPtr->getHostGrid().globalIdSet();

        // Construct the subgrid.
        subgridPtr->createBegin();

        // Loop over all elements of the host grid and use the selector to
        // choose which elements to add to the subgrid.
        auto hostGridView = subgridPtr->getHostGrid().leafGridView();
        for (const auto& e : elements(hostGridView))
            if (selector(e))
                elementsForSubgrid.insert(globalIDset.template id<0>(e));

        if (elementsForSubgrid.empty())
            DUNE_THROW(Dune::GridError, "No elements in subgrid");

        subgridPtr->insertSetPartial(elementsForSubgrid);
        subgridPtr->createEnd();

        // Return a unique pointer to the subgrid.
        return subgridPtr;
    }

    void initHostGrid_(const std::string& paramGroup)
    {
        hostGridManager_ = std::make_unique<HostGridManager>();
        hostGridManager_->init(paramGroup);
    }

    /*!
     * \brief Returns a reference to the host grid.
     */
    HostGrid& hostGrid_()
    {
        return hostGridManager_->grid();
    }

    std::unique_ptr<HostGridManager> hostGridManager_;
};

/*!
 * \ingroup InputOutput
 * \brief Provides a grid manager for SubGrids
 *        from information in the input file
 *
 * The following keys are recognized:
 * - All parameters that the host grid knows
 */
template<int dim, class HostGrid>
class GridManager<Dune::SubGrid<dim, HostGrid>>
: public SubGridManagerBase<HostGrid, GridManager<HostGrid>>
{};

/*!
 * \ingroup InputOutput
 * \brief Provides a grid manager for SubGrids
 *        from information in the input file
 *
 * The following keys are recognized:
 * - All parameters that the host grid knows
 * - Image: the image file if the sub grid is constructed from a raster image
 *          (in that case the host grid has to be any 2D YaspGrid)
 */
template<int dim, class Coordinates>
class GridManager<Dune::SubGrid<dim, Dune::YaspGrid<dim, Coordinates>>>
: public SubGridManagerBase<Dune::YaspGrid<dim, Coordinates>, GridManager<Dune::YaspGrid<dim, Coordinates>>>
{
    using ParentType = SubGridManagerBase<Dune::YaspGrid<dim, Coordinates>,
                                          GridManager<Dune::YaspGrid<dim, Coordinates>>>;
public:
    using typename ParentType::Grid;
    using ParentType::init;

    /*!
     * \brief Make the subgrid without host grid and element selector
     * This means we try to construct the element selector from the input file
     * \param paramGroup the parameter file group to check
     */
    void init(const std::string& paramGroup = "")
    {
        // check if there is an image file we can construct the element selector from
        if (hasParamInGroup(paramGroup, "Grid.Image"))
        {
            const auto imgFileName = getParamFromGroup<std::string>(paramGroup, "Grid.Image");
            const auto ext = this->getFileExtension(imgFileName);
            if (ext == "pbm")
            {
                if (dim != 2)
                    DUNE_THROW(Dune::GridError, "Portable Bitmap Format only supports dim == 2");

                // read image
                const auto img = NetPBMReader::readPBM(imgFileName);
                createGridFromImage_(img, paramGroup);
            }
            else
                DUNE_THROW(Dune::IOError, "The SubGridManager doesn't support image files with extension: *." << ext);

        }
        else
            DUNE_THROW(Dune::IOError, "SubGridManager couldn't construct element selector. Specify Grid.Image in the input file!");
    }

private:
    template<class Img>
    void createGridFromImage_(const Img& img, const std::string& paramGroup)
    {
        // check if the number of cells matches
        std::array<int, dim> cells; cells.fill(1);
        cells = getParamFromGroup<std::array<int, dim>>(paramGroup, "Grid.Cells", cells);

        if (img.header().nCols != cells[0])
            DUNE_THROW(Dune::GridError, "Host grid has wrong number of cells in x-direction. Expected "
                         << img.header().nCols << ", got " << cells[0]);
        if (img.header().nRows != cells[1])
            DUNE_THROW(Dune::GridError, "Host grid has wrong number of cells in y-direction. Expected "
                         << img.header().nRows << ", got " << cells[1]);


        // construct the host grid
        this->initHostGrid_(paramGroup);

        // check if the marker is customized, per default
        // we mark all cells that are encoded as 0
        const bool marked = getParamFromGroup<bool>(paramGroup, "Grid.Marker", 0);

        // Create the selector
        auto elementSelector = [this, &img, &marked](const auto& element)
        {
            const auto eIdx = this->hostGrid_().leafGridView().indexSet().index(element);
            return img[eIdx] == marked;
        };

        // create the grid
        this->gridPtr() = this->createGrid_(this->hostGrid_(), elementSelector, paramGroup);
        this->loadBalance();
    }
};

//! dune-subgrid doesn't have this implemented
template<int dim, class HostGrid>
class BoundaryFlag<Dune::SubGrid<dim, HostGrid>>
{
public:
    BoundaryFlag() : flag_(-1) {}

    template<class Intersection>
    BoundaryFlag(const Intersection& i) : flag_(-1) {}

    using value_type = int;

    value_type get() const
    { DUNE_THROW(Dune::NotImplemented, "Sub-grid doesn't implement boundary segment indices!"); }

private:
    int flag_;
};

} // end namespace Dumux
#endif // HAVE_DUNE_SUBGRID
#endif
