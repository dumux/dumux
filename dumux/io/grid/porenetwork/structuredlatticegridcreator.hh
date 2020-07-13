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
 * \brief Creates a network grid from a strucutred lattice. Connections can be randomly deleted.
 */
#ifndef DUMUX_IO_STRUCTURED_LATTICE_GRID_CREATOR_HH
#define DUMUX_IO_STRUCTURED_LATTICE_GRID_CREATOR_HH

#include <vector>
#include <memory>
#include <type_traits>
#include <random>

#include <dune/common/concept.hh>
#include <dune/common/exceptions.hh>
#include <dune/common/version.hh>
#include <dune/grid/common/gridfactory.hh>
#include <dune/grid/yaspgrid.hh>
#include <dune/geometry/referenceelements.hh>

// FoamGrid specific includes
#if HAVE_DUNE_FOAMGRID
#include <dune/foamgrid/foamgrid.hh>
#include <dune/foamgrid/dgffoam.hh>
#endif

#include <dumux/io/grid/gridmanager_yasp.hh>
#include <dumux/common/parameters.hh>

namespace Dumux {

namespace Concept {

/*!
 * \ingroup InputOutput
 * \brief The element selector concept
 */
template<class GlobalPosition>
struct LowDimElementSelector
{
    template<class F>
    auto require(F&& f) -> decltype(
        bool(f(std::declval<const GlobalPosition&>(), std::declval<const GlobalPosition&>()))
    );
};
} // end namespace Concept

template<int dimWorld>
class StructuredLatticeGridCreator
{
    using GridType = Dune::FoamGrid<1, dimWorld>;
    using Element = typename GridType::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using CoordScalar = typename GridType::ctype;
    using HostGrid = Dune::YaspGrid<dimWorld, Dune::TensorProductCoordinates<CoordScalar, dimWorld>>;
    using HostGridView = typename HostGrid::LeafGridView;
    using ReferenceElements = typename Dune::ReferenceElements<CoordScalar, dimWorld>;
    // grid factory expects unsigned int, will yield compiler warning for std::size_t
    using VertexIndexPair = std::pair<unsigned int, unsigned int>;

    struct Diagonal
    {
        VertexIndexPair localVertexIndices;
        unsigned int directionNumber;
    };
public:

    using Grid = GridType;

    void init(const std::string& paramGroup = "")
    {
        auto useAllElements = [](const GlobalPosition& a, const GlobalPosition& b){ return true; };
        init(useAllElements, paramGroup);
    }

    template<class LowDimElementSelector,
             typename std::enable_if_t<Dune::models<Concept::LowDimElementSelector<GlobalPosition>, LowDimElementSelector>(), int> = 0>
    void init(const LowDimElementSelector& lowDimElementSelector,
              const std::string& paramGroup = "")
    {
        paramGroup_ = paramGroup;
        removeThroatsOnBoundary_ = getParamFromGroup<std::vector<std::size_t>>(paramGroup, "Grid.RemoveThroatsOnBoundary",
                                                     std::vector<std::size_t>());

        setElementSelector_(lowDimElementSelector);
        initRandomNumberGenerator_();

        // create the host grid
        const auto hostGridInputData = getHostGridInputData_();
        using HostGrid = Dune::YaspGrid<dimWorld, Dune::TensorProductCoordinates<CoordScalar, dimWorld>>;
        using HastGridManager = GridManager<HostGrid>;
        HastGridManager hostGridManager;
        hostGridManager.init(hostGridInputData.positions, hostGridInputData.cells, hostGridInputData.grading, paramGroup_);
        hostGridView_ = std::make_unique<HostGridView>(hostGridManager.grid().leafGridView());

        hostGridLowerLeft_ = hostGridInputData.lowerLeft;
        hostGridUpperRight_ = hostGridInputData.upperRight;

        convertHostGridToLowDimGrid_();
    }

    /*!
     * \brief Call loadBalance() function of the grid.
     */
    void loadBalance()
    {
        if (Dune::MPIHelper::getCollectiveCommunication().size() > 1)
            gridPtr_->loadBalance();
    }

    /*!
     * \brief Returns a reference to the grid.
     */
    Grid& grid()
    { return *gridPtr(); }

    /*!
     * \brief Returns a reference to the grid pointer (std::shared_ptr<Grid>)
     */
    std::shared_ptr<Grid>& gridPtr()
    { return gridPtr_; }

private:

    template<class LowDimElementSelector>
    void setElementSelector_(const LowDimElementSelector& lowDimElementSelector)
    {
        if (!removeThroatsOnBoundary_.empty())
        {
            auto finalSelector = [&](const GlobalPosition& a, const GlobalPosition& b)
            {
                static const auto lambdaToRemoveThroatsOnBoundary = getLambdaToRemoveThroatsOnBoundary_();
                using std::min;
                return min(lambdaToRemoveThroatsOnBoundary(a, b), lowDimElementSelector(a, b));
            };
            considerLowDimElement_ = finalSelector;
        }
        else
            considerLowDimElement_ = lowDimElementSelector;
    }

    void initRandomNumberGenerator_()
    {
        directionProbability_ = getDirectionProbability_();

        if (hasParamInGroup(paramGroup_, "Grid.DeletionRandomNumberSeed"))
        {
            const auto seed = getParamFromGroup<std::size_t>(paramGroup_, "Grid.DeletionRandomNumberSeed");
            generator_.seed(seed);
        }
        else
        {
            std::random_device rd;
            generator_.seed(rd());
        }
    }

    void convertHostGridToLowDimGrid_()
    {
        resizeContainers_();

        for (const auto& element : elements(*hostGridView_))
        {
            const auto& geometry = element.geometry();
            const auto& refElement = ReferenceElements::general(geometry.type());
            treatEdges_(element, refElement, geometry);
            treatDiagonals_(element, refElement, geometry);
        }
        if (elementSet_.empty())
            DUNE_THROW(Dune::GridError, "Trying to create pore network with zero elements!");

        // make the actual network grid
        Dune::GridFactory<Grid> factory;
        for (auto&& vertex : vertexSet_)
            factory.insertVertex(vertex);
        for (auto&& element : elementSet_)
#if DUNE_VERSION_NEWER(DUNE_COMMON,2,7)
            factory.insertElement(Dune::GeometryTypes::cube(1), {element.first, element.second});
#else
            factory.insertElement(Dune::GeometryType(1), {element.first, element.second});
#endif

        gridPtr_ = std::shared_ptr<Grid>(factory.createGrid());
    }

    void resizeContainers_()
    {
        vertexInserted_.resize(hostGridView_->size(dimWorld), false);
        hostVertexIdxToVertexIdx_.resize(hostGridView_->size(dimWorld));
        edgeInserted_.resize(hostGridView_->size(dimWorld-1), false);
        faceInserted_.resize(hostGridView_->size(dimWorld-2), false);
    }

    template<class HostGridElement>
    bool maybeInsertElementAndVertices_(const HostGridElement& hostGridElement,
                                        const typename HostGridElement::Geometry& hostGridElementGeometry,
                                        std::size_t localVertexIdx0,
                                        std::size_t localVertexIdx1,
                                        double directionProbability)
    {
        if (neglectElement_(hostGridElementGeometry, localVertexIdx0, localVertexIdx1, directionProbability))
            return false;
        else
        {
            insertElementAndVertices_(hostGridElement, hostGridElementGeometry, localVertexIdx0, localVertexIdx1, directionProbability);
            return true;
        }
    }

    template<class HostGridElement>
    void insertElementAndVertices_(const HostGridElement& hostGridElement,
                                   const typename HostGridElement::Geometry& hostGridElementGeometry,
                                   std::size_t localVertexIdx0,
                                   std::size_t localVertexIdx1,
                                   double directionProbability)
    {
        // get global vertex indices w.r.t. host grid
        const auto vIdxGlobal0 = hostGridView_->indexSet().subIndex(hostGridElement, localVertexIdx0, dimWorld);
        const auto vIdxGlobal1 = hostGridView_->indexSet().subIndex(hostGridElement, localVertexIdx1, dimWorld);

        auto getGobalVertexIdx = [&](auto globalIdx, auto localIdx)
        {
            if (!vertexInserted_[globalIdx])
            {
                vertexInserted_[globalIdx] = true;
                hostVertexIdxToVertexIdx_[globalIdx] = vertexSet_.size();
                vertexSet_.push_back(hostGridElementGeometry.corner(localIdx));
            }

            return hostVertexIdxToVertexIdx_[globalIdx];
        };

        const auto newVertexIdx0 = getGobalVertexIdx(vIdxGlobal0, localVertexIdx0);
        const auto newVertexIdx1 = getGobalVertexIdx(vIdxGlobal1, localVertexIdx1);
        elementSet_.emplace_back(newVertexIdx0, newVertexIdx1);
    }

    auto getDirectionProbability_() const
    {
        std::array<double, numDirections_()> directionProbability;
        directionProbability.fill(0.0); // do not delete any throats

        // convenience option to delete all diagonal throats
        if (getParamFromGroup<bool>(paramGroup_, "Grid.RegularLattice", false))
        {
            directionProbability.fill(1.0); // delete all throats ...
            for (int i = 0; i < dimWorld; ++i)
                directionProbability[i] = 0.0; // ... but not the ones parallel to the main axes

            return directionProbability;
        }

        // get user input or print out help message explaining correct usage
        if (hasParamInGroup(paramGroup_, "Grid.DeletionProbability"))
        {
            try
            {
                directionProbability = getParamFromGroup<decltype(directionProbability)>(paramGroup_, "Grid.DeletionProbability");
            }
            catch(Dumux::ParameterException &e)
            {
                throwDirectionError_();
            }
        }

        return directionProbability;
    }

    void throwDirectionError_() const
    {
        // define directions (to be used for user specified anisotropy)
        Dune::FieldVector<std::string, numDirections_()> directions;
        if (dimWorld < 3) // 2D
        {
            // x, y
            directions[0] = "1: (1, 0)\n";
            directions[1] = "2: (0, 1)\n";
            // diagonals through cell midpoint
            directions[2] = "3: (1, 1)\n";
            directions[3] = "4: (1, -1)\n";
        }
        else // 3D
        {
            // x, y, z
            directions[0] = " 1: (1, 0, 0)\n";
            directions[1] = "2: (0, 1, 0)\n";
            directions[2] = "3: (0, 0, 1)\n";
            //face diagonals
            directions[3] = "4: (1, 1, 0)\n";
            directions[4] = "5: (1, -1, 0)\n";
            directions[5] = "6: (1, 0, 1)\n";
            directions[6] = "7: (1, 0, -1)\n";
            directions[7] = "8: (0, 1, 1)\n";
            directions[8] = "9: (0, 1, -1)\n";
            // diagonals through cell midpoint
            directions[9] = "10: (1, 1, 1)\n";
            directions[10] = "11: (1, 1, -1)\n";
            directions[11] = "12: (-1, 1, 1)\n";
            directions[12] = "13: (-1, -1, 1)\n";
        }
        DUNE_THROW(Dumux::ParameterException, "You must specifiy probabilities for all directions (" << numDirections_() << ") \n" << directions << "\nExample (3D):\n\n"
        << "DeletionProbability = 0.5 0.5 0 0 0 0 0 0 0 0 0 0 0 \n\n"
        << "deletes approximately 50% of all throats in x and y direction, while no deletion in any other direction takes place.\n" );
    }

    static constexpr std::size_t numDirections_()
    { return (dimWorld < 3) ? 4 : 13; }

    template<class Geometry>
    bool neglectElement_(Geometry& hostGridElementGeometry,
                         std::size_t localVertexIdx0,
                         std::size_t localVertexIdx1,
                         double directionProbability)
    {
        if (randomNumer_() < directionProbability)
            return true;

        // TODO Change order of execution: check considerLowDimElement_ before checking the random number
        // TODO This change will alter the reference solution because randomNumer_() gets called less, therefore the random numbers are different
        if (!considerLowDimElement_(hostGridElementGeometry.corner(localVertexIdx0), hostGridElementGeometry.corner(localVertexIdx1)))
            return true;

        return false;
    }

    auto randomNumer_()
    { return uniformdist_(generator_); }

    auto getHostGridInputData_() const
    {
        struct InputData
        {
            std::array<std::vector<int>, dimWorld> cells;
            std::array<std::vector<CoordScalar>, dimWorld> positions;
            std::array<std::vector<CoordScalar>, dimWorld> grading;
            GlobalPosition lowerLeft;
            GlobalPosition upperRight;
        } inputData;

        // convenience references
        auto& cells = inputData.cells;
        auto& positions = inputData.positions;
        auto& grading = inputData.grading;
        auto& lowerLeft = inputData.lowerLeft;
        auto& upperRight = inputData.upperRight;

        // try to get the pore positions explicitly ...
        for (int i = 0; i < dimWorld; ++i)
            positions[i] = getParamFromGroup<std::vector<CoordScalar>>(paramGroup_, "Grid.Positions" + std::to_string(i), std::vector<CoordScalar>{});
        if (std::none_of(positions.begin(), positions.end(), [&](auto& v){ return v.empty(); }))
        {
            for (int i = 0; i < dimWorld; ++i)
            {
                cells[i].resize(positions[i].size()-1, 1.0);
                grading[i].resize(positions[i].size()-1, 1.0);
            }
        }
        else // .. or calculate positions from input data
        {
            const auto lowerLeft = getParamFromGroup<GlobalPosition>(paramGroup_, "Grid.LowerLeft", GlobalPosition(0.0));
            const auto upperRight = getParamFromGroup<GlobalPosition>(paramGroup_, "Grid.UpperRight");
            const auto numPores = getParamFromGroup<std::array<unsigned int, dimWorld>>(paramGroup_, "Grid.NumPores");
            for (int i = 0; i < dimWorld; ++i)
            {
                positions[i].push_back(lowerLeft[i]);
                positions[i].push_back(upperRight[i]);
                cells[i].push_back(numPores[i] - 1);
                grading[i].resize(positions[i].size()-1, 1.0);
                grading[i] = getParamFromGroup<std::vector<CoordScalar>>(paramGroup_, "Grid.Grading" + std::to_string(i), grading[i]);
            }
        }

        // get the lower left position
        lowerLeft = [&]()
        {
            GlobalPosition result;
            for (int i = 0; i < dimWorld; ++i)
                result[i] = positions[i][0];
            return result;
        }();

        // get the upper right position
        upperRight = [&]()
        {
            GlobalPosition result;
            for (int i = 0; i < dimWorld; ++i)
                result[i] = positions[i].back();
            return result;
        }();

        return inputData;
    }

    template<class HostGridElement, class ReferenceElement, class Geometry>
    void treatEdges_(const HostGridElement& element, const ReferenceElement& refElement, const Geometry& geometry)
    {
        // go over all edges and add them as elements if they passed all the tests
        for (unsigned int edgeIdx = 0; edgeIdx < element.subEntities(dimWorld-1); ++edgeIdx)
        {
            const auto vIdxLocal0 = refElement.subEntity(edgeIdx, dimWorld-1, 0, dimWorld);
            const auto vIdxLocal1 = refElement.subEntity(edgeIdx, dimWorld-1, 1, dimWorld);
            const auto edgeIdxGlobal = hostGridView_->indexSet().subIndex(element, edgeIdx, dimWorld-1);

            if (edgeInserted_[edgeIdxGlobal])
                continue;
            else
                edgeInserted_[edgeIdxGlobal] = true;

            std::size_t directionNumber = 0;

            if(dimWorld == 2 ) // 2D
            {
                if(edgeIdx < 2) // y-direction
                    directionNumber = 0;
                else // x-direction
                    directionNumber = 1;
            }
            else // 3D
            {
                if(edgeIdx < 4) // z-direction
                    directionNumber = 2;
                else if(edgeIdx == 4 || edgeIdx == 5 || edgeIdx == 8 || edgeIdx == 9) // y-direction
                    directionNumber = 1;
                else if(edgeIdx == 6 || edgeIdx == 7 || edgeIdx == 10 || edgeIdx == 11) // x-direction
                    directionNumber = 0;
            }
            maybeInsertElementAndVertices_(element, geometry, vIdxLocal0, vIdxLocal1, directionProbability_[directionNumber]);
        }
    }

    template<class HostGridElement, class ReferenceElement, class Geometry>
    void treatDiagonals_(const HostGridElement& element, const ReferenceElement& refElement, const Geometry& geometry)
    {
        if constexpr (dimWorld < 3)
            treatDiagonalConnections2D_(element, geometry);
        else
            treatDiagonalConnections3D_(element, refElement, geometry);
    }

    template<class HostGridElement, class Geometry>
    void treatDiagonalConnections2D_(const HostGridElement& element,
                                     const Geometry& geometry)
    {
         static constexpr std::array<Diagonal, 2> diagonals{ Diagonal{std::make_pair(0, 3), 2},
                                                             Diagonal{std::make_pair(1, 2), 3} };

        treatIntersectingDiagonals_(element, geometry, diagonals);
    }

    template<class HostGridElement, class ReferenceElement, class Geometry>
    void treatDiagonalConnections3D_(const HostGridElement& element,
                                     const ReferenceElement& refElement,
                                     const Geometry& geometry)
    {
        // set diagonals on host grid element faces
        for (auto faceIdx = 0; faceIdx < element.subEntities(dimWorld-2); ++faceIdx)
        {
            const auto faceIdxGlobal = hostGridView_->indexSet().subIndex(element, faceIdx, dimWorld-2);
            // face already checked?
            if (faceInserted_[faceIdxGlobal])
                continue;
            else
                faceInserted_[faceIdxGlobal] = true;

            // get local vertex indices
            std::array<unsigned int, 4> vIdxLocal;
            for (int i = 0; i < 4; i++)
                vIdxLocal[i]  = refElement.subEntity(faceIdx, dimWorld-2, i, dimWorld);

            const auto directionNumbers = [&]()
            {
                if (faceIdx < 2)
                    return std::make_pair<unsigned int, unsigned int>(8,7);
                else if (faceIdx < 4)
                    return std::make_pair<unsigned int, unsigned int>(6,5);
                else
                    return std::make_pair<unsigned int, unsigned int>(4,3);
            }();

            const std::array<Diagonal, 2> diagonals{ Diagonal{std::make_pair(vIdxLocal[1], vIdxLocal[2]), directionNumbers.first},
                                                     Diagonal{std::make_pair(vIdxLocal[0], vIdxLocal[3]), directionNumbers.second} };

            treatIntersectingDiagonals_(element, geometry, diagonals);
        }

        static constexpr std::array<Diagonal, 4> diagonals{ Diagonal{std::make_pair(0, 7), 9},
                                                            Diagonal{std::make_pair(3, 4), 10},
                                                            Diagonal{std::make_pair(1, 6), 11},
                                                            Diagonal{std::make_pair(2, 5), 12}, };

        treatIntersectingDiagonals_(element, geometry, diagonals);
    }

    template<class HostGridElement, class Geometry, class Diagonals>
    void treatIntersectingDiagonals_(const HostGridElement& element,
                                     const Geometry& geometry,
                                     const Diagonals& diagonals)
    {
        static const bool allowIntersectingDiagonals = getParamFromGroup<bool>(paramGroup_, "Grid.AllowIntersectingDiagonals", true);

        auto treat = [&](const Diagonal& diagonal)
        {
            return maybeInsertElementAndVertices_(element, geometry,
                                                  diagonal.localVertexIndices.first, diagonal.localVertexIndices.second,
                                                  directionProbability_[diagonal.directionNumber]);
        };

        if (allowIntersectingDiagonals)
        {
            // insert all diagonals
            for (const auto& diagonal : diagonals)
                treat(diagonal);
        }
        else
        {
            auto order = createOrderedList_(diagonals.size());
            std::shuffle(order.begin(), order.end(), generator_);
            for (auto i : order)
            {
                const auto& diagonal = diagonals[i];
                const bool inserted = treat(diagonal);
                if (inserted)
                    return;
            }
        }
    }

    std::vector<std::size_t> createOrderedList_(const std::size_t size) const
    {
        std::vector<std::size_t> result(size);
        std::iota(result.begin(), result.end(), 0);
        return result;
    }

    auto getLambdaToRemoveThroatsOnBoundary_() const
    {
        auto negletElementsOnBoundary = [&](const GlobalPosition& a, const GlobalPosition& b)
        {
            const auto center = 0.5 * (a + b);
            const auto eps = (a-b).two_norm() * 1e-5;

            bool neglectElement = false;
            for (auto i : removeThroatsOnBoundary_)
            {
                switch(i)
                {
                    case 0: neglectElement = center[0] < hostGridLowerLeft_[0] + eps; break;
                    case 1: neglectElement = center[0] > hostGridUpperRight_[0] - eps; break;
                    case 2: neglectElement = center[1] < hostGridLowerLeft_[1] + eps; break;
                    case 3: neglectElement = center[1] > hostGridUpperRight_[1] - eps; break;
                    case 4: neglectElement = center[2] < hostGridLowerLeft_[2] + eps; break;
                    case 5: neglectElement = center[2] > hostGridUpperRight_[2] - eps; break;
                }

                if (neglectElement)
                    return false;
            }

            return true;
        };

        return negletElementsOnBoundary;
    }

    std::string paramGroup_;
    GlobalPosition hostGridLowerLeft_;
    GlobalPosition hostGridUpperRight_;
    std::vector<std::size_t> removeThroatsOnBoundary_;
    std::unique_ptr<const HostGridView> hostGridView_;
    std::function<bool(const GlobalPosition&, const GlobalPosition&)> considerLowDimElement_;

    std::vector<GlobalPosition> vertexSet_;
    std::vector<VertexIndexPair> elementSet_;

    std::vector<bool> vertexInserted_;
    std::vector<std::size_t> hostVertexIdxToVertexIdx_;
    std::vector<std::size_t> edgeInserted_;
    std::vector<std::size_t> faceInserted_;

    mutable std::mt19937 generator_;
    std::uniform_real_distribution<> uniformdist_{0, 1};
    std::array<double, numDirections_()> directionProbability_;

    std::shared_ptr<Grid> gridPtr_;
};

} // namespace Dumux

#endif // DUMUX_IO_STRUCTURED_LATTICE_GRID_CREATOR_HH
