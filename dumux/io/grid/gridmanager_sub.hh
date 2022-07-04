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
#include <utility>
#include <algorithm>
#include <stack>
#include <vector>

#include <dune/common/shared_ptr.hh>
#include <dune/common/concept.hh>

// SubGrid specific includes
#if HAVE_DUNE_SUBGRID
#include <dune/subgrid/subgrid.hh>
#include <dumux/io/rasterimagereader.hh>
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
     * \brief Call loadBalance() function of the grid.
     */
    void loadBalance()
    {
        if (Dune::MPIHelper::getCommunication().size() > 1)
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

            postProcessBinaryMask_(buffer, cells, paramGroup);

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
            this->gridPtr() = this->createGrid_(this->hostGrid_(), elementSelector, paramGroup);
            this->loadBalance();
        }
        else
        {
            std::cerr << "No element selector provided and Grid.Image or Grid.BinaryMask not found" << std::endl;
            std::cerr << "Constructing sub-grid with all elements of the host grid" << std::endl;

            // create host grid
            using GlobalPosition = typename ParentType::Grid::template Codim<0>::Geometry::GlobalCoordinate;
            const auto lowerLeft = getParamFromGroup<GlobalPosition>(paramGroup, "Grid.LowerLeft", GlobalPosition(0.0));
            const auto upperRight = getParamFromGroup<GlobalPosition>(paramGroup, "Grid.UpperRight");
            const auto cells = getParamFromGroup<std::array<int, dim>>(paramGroup, "Grid.Cells");
            this->initHostGrid_(lowerLeft, upperRight, cells, paramGroup);
            const auto elementSelector = [](const auto& element){ return true; };
            // create the grid
            this->gridPtr() = this->createGrid_(this->hostGrid_(), elementSelector, paramGroup);
            this->loadBalance();
        }
    }

private:
    template<class Img>
    void createGridFromImage_(const Img& img, const std::string& paramGroup)
    {
        // get the number of cells
        const std::array<int, dim> cells{static_cast<int>(img.header().nCols), static_cast<int>(img.header().nRows)};

        // get the corner coordinates
        const auto [lowerLeft, upperRight] = [&]()
        {
            using GlobalPosition = typename ParentType::Grid::template Codim<0>::Geometry::GlobalCoordinate;
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

        // Create the selector
        auto elementSelector = [this, &img, &marked](const auto& element)
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

        // create the grid
        this->gridPtr() = this->createGrid_(this->hostGrid_(), elementSelector, paramGroup);

        this->loadBalance();
    }

private:
    void postProcessBinaryMask_(std::vector<char>& mask, const std::array<int, dim>& cells, const std::string& paramGroup) const
    {
        if (hasParamInGroup(paramGroup, "Grid.KeepOnlyConnected"))
        {
            std::vector<int> marker(mask.size(), 0);

            const std::string direction = getParamFromGroup<std::string>(paramGroup, "Grid.KeepOnlyConnected");
            if (direction == "Y")
            {
                std::cout << "Keeping only cells connected to both boundaries in y-direction" << std::endl;
                for (unsigned int i = 0; i < cells[0]; ++i)
                    for (unsigned int k = 0; k < cells[2]; ++k)
                        fill_(mask, marker, cells, i, 0, k, 1);

                for (unsigned int i = 0; i < cells[0]; ++i)
                    for (unsigned int k = 0; k < cells[2]; ++k)
                        fill_(mask, marker, cells, i, cells[1]-1, k, 2);
            }
            else if (direction == "Z")
            {
                std::cout << "Keeping only cells connected to both boundaries in z-direction" << std::endl;
                for (unsigned int i = 0; i < cells[0]; ++i)
                    for (unsigned int j = 0; j < cells[1]; ++j)
                        fill_(mask, marker, cells, i, j, 0, 1);

                for (unsigned int i = 0; i < cells[0]; ++i)
                    for (unsigned int j = 0; j < cells[1]; ++j)
                        fill_(mask, marker, cells, i, j, cells[2]-1, 2);
            }
            else
            {
                std::cout << "Keeping only cells connected to both boundaries in x-direction" << std::endl;
                for (unsigned int j = 0; j < cells[1]; ++j)
                    for (unsigned int k = 0; k < cells[2]; ++k)
                        fill_(mask, marker, cells, 0, j, k, 1);

                for (unsigned int j = 0; j < cells[1]; ++j)
                    for (unsigned int k = 0; k < cells[2]; ++k)
                        fill_(mask, marker, cells, cells[0]-1, j, k, 2);
            }

            for (unsigned int i = 0; i < cells[0]; ++i)
            {
                for (unsigned int j = 0; j < cells[1]; ++j)
                {
                    for (unsigned int k = 0; k < cells[2]; ++k)
                    {
                        const auto ijk = getIJK_(i, j, k, cells);
                        if (marker[ijk] < 2)
                            mask[ijk] = false;
                    }
                }
            }
        }
    }

    void fill_(std::vector<char>& mask, std::vector<int>& marker, const std::array<int, dim>& cells,
               unsigned int i, unsigned int j, unsigned int k, int m) const
    {
        std::stack<std::tuple<unsigned int, unsigned int, unsigned int>> st;
        st.push(std::make_tuple(i, j, k));
        while(!st.empty())
        {
            std::tie(i, j, k) = st.top();
            st.pop();

            if (i >= cells[0] || j >= cells[1] || k >= cells[2])
                continue;

            const auto ijk = getIJK_(i, j, k, cells);
            if (!mask[ijk] || marker[ijk] >= m)
               continue;

            marker[ijk] += 1;

            // we rely on -1 wrapping over for unsigned int = 0
            st.push(std::make_tuple(i+1, j, k));
            st.push(std::make_tuple(i-1, j, k));
            st.push(std::make_tuple(i, j+1, k));
            st.push(std::make_tuple(i, j-1, k));
            st.push(std::make_tuple(i, j, k+1));
            st.push(std::make_tuple(i, j, k-1));
        }
    }

    int getIJK_(int i, int j, int k, const std::array<int, dim>& cells) const
    {
        return i + j*cells[0] + k*cells[0]*cells[1];
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
