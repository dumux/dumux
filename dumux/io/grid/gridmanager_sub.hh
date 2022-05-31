// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup InputOutput
 * \brief Grid manager specialization for SubGrid
 */
#ifndef DUMUX_IO_GRID_MANAGER_SUB_HH
#define DUMUX_IO_GRID_MANAGER_SUB_HH

#include <memory>
#include <utility>
#include <algorithm>

#include <dune/common/shared_ptr.hh>
#include <dune/common/concept.hh>
#include <dune/grid/yaspgrid.hh>

// SubGrid specific includes
#if HAVE_DUNE_SUBGRID
#include <dune/subgrid/subgrid.hh>
#include <dumux/io/rasterimagereader.hh>
#include <dumux/io/rasterimagewriter.hh>
#endif

#include <dumux/io/grid/gridmanager_base.hh>

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
    using GlobalPosition = typename HostElement::Geometry::GlobalCoordinate;

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
        this->gridPtr() = std::make_unique<Grid>(hostGrid);
        updateSubGrid_(selector);
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
        this->gridPtr() = std::make_unique<Grid>(hostGridManager_->grid());
        updateSubGrid_(selector);
        loadBalance();
    }

    /*!
     * \brief Call loadBalance() function of the grid.
     */
    void loadBalance()
    {
        if (Dune::MPIHelper::getCommunication().size() > 1)
            this->grid().loadBalance();
    }

    /*!
     * \brief Update the existing subgrid.
     */
    template<class ES,
             typename std::enable_if_t<Dune::models<Concept::ElementSelector<HostElement>, ES>(), int> = 0>
    void update(const ES& selector)
    {
        updateSubGrid_(selector);
        loadBalance();
    }

protected:

    /*!
     * \brief Update the subgrid.
     */
    template<class ES,
             typename std::enable_if_t<Dune::models<Concept::ElementSelector<HostElement>, ES>(), int> = 0>
    void updateSubGrid_(const ES& selector)
    {
        auto& subgridPtr = this->gridPtr();

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
    }

    void initHostGrid_(const std::string& paramGroup)
    {
        hostGridManager_ = std::make_unique<HostGridManager>();
        hostGridManager_->init(paramGroup);
    }

    void initHostGrid_(const GlobalPosition& lowerLeft,
                       const GlobalPosition& upperRight,
                       const std::array<int, dim>& cells,
                       const std::string& paramGroup,
                       const int overlap = 1)
    {
        hostGridManager_ = std::make_unique<HostGridManager>();
        hostGridManager_->init(lowerLeft, upperRight, cells, paramGroup, overlap);
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
        else if (hasParamInGroup(paramGroup, "Grid.BinaryMask"))
        {
            const auto maskFileName = getParamFromGroup<std::string>(paramGroup, "Grid.BinaryMask");
            std::ifstream mask(maskFileName, std::ios_base::binary);
            std::vector<char> buffer(std::istreambuf_iterator<char>(mask), std::istreambuf_iterator<char>{});
            const auto cells = getParamFromGroup<std::array<int, dim>>(paramGroup, "Grid.Cells");
            if (std::accumulate(cells.begin(), cells.end(), 1, std::multiplies<int>{}) != buffer.size())
                DUNE_THROW(Dune::IOError, "Grid dimensions doesn't match number of cells specified");

            using GlobalPosition = typename ParentType::Grid::template Codim<0>::Geometry::GlobalCoordinate;
            const auto lowerLeft = GlobalPosition(0.0);
            const auto upperRight = [&]{
                if (hasParamInGroup(paramGroup, "Grid.PixelDimensions"))
                {
                    auto upperRight = getParamFromGroup<GlobalPosition>(paramGroup, "Grid.PixelDimensions");
                    for (int i = 0; i < upperRight.size(); ++i)
                        upperRight[i] *= cells[i];
                    upperRight += lowerLeft;
                    return upperRight;
                }
                else
                    return getParamFromGroup<GlobalPosition>(paramGroup, "Grid.UpperRight");
            }();

            // make sure there is no grid refinement specified
            if (getParamFromGroup<int>(paramGroup, "Grid.Refinement", 0) > 0)
                DUNE_THROW(Dune::IOError,
                    "Binary mask doesn't work together with Grid.Refinement."
                    << " Use grid.globalRefine() after grid construction."
                );

            // construct the host grid
            this->initHostGrid_(lowerLeft, upperRight, cells, paramGroup);

            // check if the marker is customized, per default
            // we use all cells that are encoded as 0
            const char marker = getParamFromGroup<char>(paramGroup, "Grid.Marker", 0);
            const auto elementSelector = [&](const auto& element)
            {
                const auto eIdx = this->hostGrid_().leafGridView().indexSet().index(element);
                return buffer[eIdx] != marker;
            };

            // create the grid
            this->gridPtr() = std::make_unique<Grid>(this->hostGrid_());
            this->updateSubGrid_(elementSelector);
            this->loadBalance();
        }
        else
            DUNE_THROW(Dune::IOError, "SubGridManager couldn't construct element selector."
                << " Specify Grid.Image or Grid.BinaryMask in the input file!");
    }

private:
    template<class Img>
    void createGridFromImage_(const Img& img, const std::string& paramGroup)
    {
        using GlobalPosition = typename ParentType::Grid::template Codim<0>::Geometry::GlobalCoordinate;
        const bool repeated = hasParamInGroup(paramGroup,"Grid.Repeat");

        // get the number of cells
        const std::array<int, dim> imageDimensions{static_cast<int>(img.header().nCols), static_cast<int>(img.header().nRows)};
        std::array<int, dim> cells{imageDimensions[0], imageDimensions[1]};

        std::array<int, dim> repeatsDefault; repeatsDefault.fill(1);
        const auto numRepeats = getParamFromGroup<std::array<int, dim>>(paramGroup, "Grid.Repeat", repeatsDefault);
        for (int i = 0; i < dim; i++)
            cells[i] = cells[i] * numRepeats[i];

        // get the corner coordinates
        const auto [lowerLeft, upperRight] = [&]()
        {
            const auto lowerLeft = getParamFromGroup<GlobalPosition>(paramGroup, "Grid.LowerLeft", GlobalPosition(0.0));
            if (hasParamInGroup(paramGroup, "Grid.PixelDimensions"))
            {
                auto upperRight = getParamFromGroup<GlobalPosition>(paramGroup, "Grid.PixelDimensions");
                for (int i = 0; i < upperRight.size(); ++i)
                    upperRight[i] *= cells[i];
                upperRight += lowerLeft;
                return std::make_pair(lowerLeft, upperRight);
            }
            else
                return std::make_pair(lowerLeft, getParamFromGroup<GlobalPosition>(paramGroup, "Grid.UpperRight"));
        }();

        // construct the host grid
        this->initHostGrid_(lowerLeft, upperRight, cells, paramGroup);

        // check if the marker is customized, per default
        // we mark all cells that are encoded as 0
        const bool marked = getParamFromGroup<bool>(paramGroup, "Grid.Marker", false);

        // Create the element selector for a single image
        const auto elementSelector = [&](const auto& element)
        {
            auto eIdx = this->hostGrid_().leafGridView().indexSet().index(element);

            // if the hostgrid was refined, get the index of the original, un-refined
            // host grid element which fits with the image's indices
            if (element.hasFather())
            {
                auto level0Element = element.father();
                while (level0Element.hasFather())
                    level0Element = level0Element.father();

                assert(level0Element.level() == 0);
                eIdx = this->hostGrid_().levelGridView(0).indexSet().index(level0Element);
            }
            return img[eIdx] == marked;
        };

        // Create the element selector for a repeated image
        const auto repeatedElementSelector = [&](const auto& element)
        {
            auto eIdx = this->hostGrid_().leafGridView().indexSet().index(element);

            // if the hostgrid was refined, get the index of the original, un-refined
            // host grid element which fits with the image's indices
            if (element.hasFather())
            {
                auto level0Element = element.father();
                while (level0Element.hasFather())
                    level0Element = level0Element.father();

                assert(level0Element.level() == 0);
                eIdx = this->hostGrid_().levelGridView(0).indexSet().index(level0Element);
            }

            // figure out the size of the repeat, and the size of the target repeated grid
            const int numCols = imageDimensions[0];
            const int numRows = imageDimensions[1];

            // map the eIdx to the original img index
            const int repeatUnitIndex = eIdx % (numCols * numRepeats[0] * numRows);
            const int imgI = repeatUnitIndex % numCols;
            const int imgJ = repeatUnitIndex / (numCols * numRepeats[0]);
            return img[ (imgJ * numCols + imgI) ] == marked;
        };

        // create the grid
        if (repeated)
        {
            this->gridPtr() = std::make_unique<Grid>(this->hostGrid_());
            this->updateSubGrid_(repeatedElementSelector);
        }
        else
        {
            this->gridPtr() = std::make_unique<Grid>(this->hostGrid_());
            this->updateSubGrid_(elementSelector);
        }

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
