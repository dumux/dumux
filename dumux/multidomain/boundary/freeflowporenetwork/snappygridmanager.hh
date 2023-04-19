// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup FreeFlowPoreNetworkCoupling
 * \brief A grid creator that matches a free-flow grid to a PNM grid.
 */
#ifndef DUMUX_MULTIDOMAIN_BOUNDARY_FREEFLOW_PORENETWORK_SNAPPY_GRID_MANAGER_HH
#define DUMUX_MULTIDOMAIN_BOUNDARY_FREEFLOW_PORENETWORK_SNAPPY_GRID_MANAGER_HH

#include <bitset>
#include <optional>
#include <dune/common/float_cmp.hh>
#include <dune/geometry/axisalignedcubegeometry.hh>
#include <dumux/geometry/intersectspointgeometry.hh>
#include <dumux/io/grid/gridmanager.hh>

namespace Dumux::PoreNetwork {

namespace Detail {
/*!
 * \ingroup FreeFlowPoreNetworkCoupling
 * \brief A helper for the grid creator that matches a free-flow grid to a PNM grid.
 */
template<class GridView>
class SnappyGridManagerHelper
{
    using Scalar = typename GridView::ctype;
    static constexpr auto dim = GridView::dimension;
    static constexpr auto dimWorld = GridView::dimensionworld;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
    using Line = Dune::AxisAlignedCubeGeometry< Scalar, 1, dimWorld>;
    using Plane = Dune::AxisAlignedCubeGeometry< Scalar, dimWorld-1, dimWorld>;

public:

    /*!
     * \brief Returns the direction index of a unit vector (0 = x, 1 = y, 2 = z)
     */
    static unsigned int directionIndex(const GlobalPosition& vector)
    {
        const auto eps = 1e-8*vector.two_norm();
        return std::find_if(vector.begin(), vector.end(), [eps](const auto& x) { return std::abs(x) > eps; } ) - vector.begin();
    }

    static auto couplingPlaneBoundingBox(const GlobalPosition& gridLowerLeft,
                                         const GlobalPosition& gridUpperRight,
                                         const GlobalPosition& couplingPlaneNormal,
                                         const std::string& modelParamGroup)
    {
        const auto couplingPlaneNormalDirectionIndex = directionIndex(couplingPlaneNormal);

        // get the spatial extent of the coupling plane
        const auto couplingPlaneLowerLeft = getParamFromGroup<GlobalPosition>(modelParamGroup, "Grid.CouplingPlaneLowerLeft", gridLowerLeft);
        const auto couplingPlaneUpperRight = getParamFromGroup<GlobalPosition>(modelParamGroup, "Grid.CouplingPlaneUpperRight",
                                                                               [&]()
                                                                               {
                                                                                   GlobalPosition tmp(gridUpperRight);
                                                                                   tmp[couplingPlaneNormalDirectionIndex] = gridLowerLeft[couplingPlaneNormalDirectionIndex];
                                                                                   return tmp;
                                                                               }());
        return std::make_pair(couplingPlaneLowerLeft, couplingPlaneUpperRight);
    }

    static Plane makeCouplingPlane(const GlobalPosition& gridLowerLeft,
                                   const GlobalPosition& gridUpperRight,
                                   const GlobalPosition& couplingPlaneNormal,
                                   const std::string& modelParamGroup)
    {
        auto [couplingPlaneLowerLeft, couplingPlaneUpperRight] = couplingPlaneBoundingBox(gridLowerLeft, gridUpperRight, couplingPlaneNormal, modelParamGroup);
        const auto couplingPlaneNormalDirectionIndex = directionIndex(couplingPlaneNormal);
        auto inPlaneAxes = std::move(std::bitset<dimWorld>{}.set());
        inPlaneAxes.set(couplingPlaneNormalDirectionIndex, false);
        return Plane(std::move(couplingPlaneLowerLeft), std::move(couplingPlaneUpperRight), inPlaneAxes);
    }

    /*!
     * \brief Creates lines parallel to the bounding box axis of the coupling plane.
     *        Creates one line for dim == 2 and two lines for dim == 3. The array size is dim, therefore the entry
     *        for the normal to the coupling plane is empty.
     *
     * \param gridLowerLeft The lower left position of the bulk grid
     * \param gridUpperRight The upper right position of the bulk grid
     * \param couplingPlaneNormal The normal vector of the coupling plane
     * \param modelParamGroup The bulk model parameter group
     */
    static std::array<std::optional<Line>, dim> makeAxisParallelLinesFromCouplingPlane(const GlobalPosition& gridLowerLeft,
                                                                                       const GlobalPosition& gridUpperRight,
                                                                                       const GlobalPosition& couplingPlaneNormal,
                                                                                       const std::string& modelParamGroup)
    {
        const auto couplingPlaneNormalDirectionIndex = directionIndex(couplingPlaneNormal);

        // get the spatial extent of the coupling plane
        auto [a, b] = couplingPlaneBoundingBox(gridLowerLeft, gridUpperRight, couplingPlaneNormal, modelParamGroup);
        // structured bindings cannot be captured by a lambda, therefore we introduce the helper variables a and b here
        const auto& couplingPlaneLowerLeft = a;
        const auto& couplingPlaneUpperRight = b;

        // Create an array of lines. The entry for dim == couplingPlaneNormalDirectionIndex remains empty
        std::array<std::optional<Line>, dim> lines;
        for (int i = 0; i < dim; ++i)
        {
            if (i == couplingPlaneNormalDirectionIndex)
                continue;

            lines[i] = [&]()
            {
                GlobalPosition tmp(couplingPlaneLowerLeft);
                tmp[i] = couplingPlaneUpperRight[i];
                auto axis = std::bitset<dimWorld>();
                axis.set(i, true);
                return Line(couplingPlaneLowerLeft, tmp, axis);
            }();
        }

        return lines;
    }

    /*!
     * \brief Returns the lowDim positions intersecting with a given line.
     *
     * \param bulkGridLowerLeft The lower left position of the bulk grid
     * \param bulkGridUpperRight The upper right position of the bulk grid
     * \param couplingPlaneNormal The normal vector of the coupling plane
     * \param lowDimGridView The lowDim gridView
     * \param lowDimGridData The lowDim grid data
     * \param modelParamGroup The bulk model parameter group
     */
    template<class LowDimGridView, class LowDimGridData>
    static auto getPointsOnLine(const Dune::FieldVector<Scalar, 3>& bulkGridLowerLeft,
                                const Dune::FieldVector<Scalar, 3>& bulkGridUpperRight,
                                const Dune::FieldVector<Scalar, 3>& couplingPlaneNormal,
                                const LowDimGridView& lowDimGridView,
                                const LowDimGridData& lowDimGridData,
                                const std::string& modelParamGroup)
    {
        const auto couplingPlane = makeCouplingPlane(bulkGridLowerLeft, bulkGridUpperRight, couplingPlaneNormal, modelParamGroup);
        const auto couplingPlaneNormalDirectionIndex = directionIndex(couplingPlaneNormal);
        using ScalarVector = std::vector<Scalar>;

        std::array<ScalarVector, 3> coordinates;
        for (auto& c : coordinates)
            c.reserve(lowDimGridView.size(1));

        ScalarVector poreRadius;
        poreRadius.reserve(lowDimGridView.size(1));

        auto makeUnique = [](ScalarVector& v)
        {
            std::sort(v.begin(), v.end());
            v.erase(std::unique(v.begin(), v.end()), v.end());
        };

        for (const auto& vertex : vertices(lowDimGridView))
        {
            const auto& pos = vertex.geometry().center();
            if (intersectsPointGeometry(pos, couplingPlane))
            {
                poreRadius.push_back(lowDimGridData.parameters(vertex)[lowDimGridData.parameterIndex("PoreInscribedRadius")]);
                for (int i = 0; i < 3; ++i)
                    coordinates[i].push_back(pos[i]);
            }
        }

        if (std::any_of(poreRadius.begin(), poreRadius.end(), [&](auto& r) { return Dune::FloatCmp::eq(r, poreRadius[0]); }))
            DUNE_THROW(Dune::GridError, "SnappyGridCreator in 3D currently only supports equal pore radii at the interface");

        for (auto& c : coordinates)
            makeUnique(c);

        struct Data { Scalar pos; Scalar radius; };
        std::array<std::optional<std::vector<Data>>, 3> result;

        for (int i = 0; i < 3; ++i)
        {
            if (i != couplingPlaneNormalDirectionIndex)
            {
                std::vector<Data> tmp(coordinates[i].size());
                for (int poreIdx = 0; poreIdx < tmp.size(); ++poreIdx)
                    tmp[poreIdx] = Data{coordinates[i][poreIdx], poreRadius[0]};

                result[i] = std::move(tmp);
            }
        }

        return result;
    }

    /*!
     * \brief Returns the lowDim positions intersecting with a given line.
     *
     * \param bulkGridLowerLeft The lower left position of the bulk grid
     * \param bulkGridUpperRight The upper right position of the bulk grid
     * \param couplingPlaneNormal The normal vector of the coupling plane
     * \param lowDimGridView The lowDim gridView
     * \param lowDimGridData The lowDim grid data
     * \param modelParamGroup The bulk model parameter group
     */
    template<class LowDimGridView, class LowDimGridData>
    static auto getPointsOnLine(const Dune::FieldVector<Scalar, 2>& bulkGridLowerLeft,
                                const Dune::FieldVector<Scalar, 2>& bulkGridUpperRight,
                                const Dune::FieldVector<Scalar, 2>& couplingPlaneNormal,
                                const LowDimGridView& lowDimGridView,
                                const LowDimGridData& lowDimGridData,
                                const std::string& modelParamGroup)
    {
        std::vector<bool> vertexVisited(lowDimGridView.size(LowDimGridView::dimension), false);

        const auto axisParallelLines = makeAxisParallelLinesFromCouplingPlane(bulkGridLowerLeft, bulkGridUpperRight, couplingPlaneNormal, modelParamGroup);

        //The line to check for intersections
        const auto line = [&]()
        {
            for (const auto& l : axisParallelLines)
                if (l.has_value())
                    return *l;
            DUNE_THROW(Dune::InvalidStateException, "No line found");
        }();

        const auto orientation = line.corner(1) - line.corner(0);
        const auto dirIdx  = directionIndex(orientation);

        struct Data { Scalar pos; Scalar radius;};
        std::vector<Data> data;

        for (const auto& element : elements(lowDimGridView))
        {
            for (std::size_t localVertexIdx = 0; localVertexIdx < 2; ++localVertexIdx)
            {
                if (const auto& pos = element.geometry().corner(localVertexIdx); intersectsPointGeometry(pos, line))
                {
                    const auto& vertex = element.template subEntity<LowDimGridView::dimension>(localVertexIdx);
                    const auto vIdx = lowDimGridView.indexSet().index(vertex);
                    const auto poreRadius = lowDimGridData.parameters(vertex)[lowDimGridData.parameterIndex("PoreInscribedRadius")];
                    assert(poreRadius > 0.0);

                    if (vertexVisited[vIdx])
                        continue;
                    else
                        vertexVisited[vIdx] = true;

                    data.push_back({pos[dirIdx], poreRadius});
                }
            }
        }

        std::sort(data.begin(), data.end(),
                  [](const auto& pore0, const auto& pore1)
                  {
                      return pore0.pos < pore1.pos;
                  });

        // check if any intersections were found
        if (data.empty())
        {
            std::cout << "line boundaries: " << std::endl;
            for (int k = 0; k < line.corners(); ++k)
                std::cout << "(" << line.corner(k) << ")  ";

            std::cout << std::endl;
            DUNE_THROW(Dune::GridError, "No coupled pores found");
        }

        std::array<std::optional<std::vector<Data>>, 2> result;
        result[dirIdx] = data;
        return result;
    }
};

} // end namespace Detail

/*!
 * \ingroup FreeFlowPoreNetworkCoupling
 * \brief  A grid creator that matches a free-flow grid to a PNM grid.
 */
template<int dim, class OtherGridCreator, class DiscMethod = DiscretizationMethods::None>
class SnappyGridManager : public Dumux::GridManager<Dune::YaspGrid<dim, Dune::TensorProductCoordinates<typename OtherGridCreator::Grid::ctype, OtherGridCreator::Grid::LeafGridView::dimensionworld>>>
{
    using Scalar = typename OtherGridCreator::Grid::ctype;
    using OtherGrid = typename OtherGridCreator::Grid;
    using IntVector = std::vector<int>;
    using ScalarVector = std::vector<Scalar>;

    static constexpr auto dimWorld = OtherGrid::LeafGridView::dimensionworld;
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

    using ParentType = Dumux::GridManager<Dune::YaspGrid<dim, Dune::TensorProductCoordinates<typename OtherGridCreator::Grid::ctype, OtherGridCreator::Grid::LeafGridView::dimensionworld>>>;
    using LowDimGridData = typename OtherGridCreator::GridData;

    struct GridConstructionData
    {
            std::array<ScalarVector, dim> positions;
            std::array<IntVector, dim> cells;
            std::array<ScalarVector, dim> grading;
            std::array<std::optional<ScalarVector>, dim> interfacePositions;
    };

public:
    using Grid = Dune::YaspGrid<dim, Dune::TensorProductCoordinates<Scalar, dim>>;
    using ParentType::ParentType;

    //! Make the grid
    void init(const OtherGrid& otherGrid, const LowDimGridData& otherData, const std::string& modelParamGroup = "")
    {
        modelParamGroup_ = modelParamGroup;
        std::array<ScalarVector, dim> positions;
        std::array<IntVector, dim> cells;
        std::array<ScalarVector, dim> grading;
        std::array<std::optional<ScalarVector>, dim> interFacePositions;
        using SnappyGridManagerHelper = Detail::SnappyGridManagerHelper<typename Grid::LeafGridView>;

        const auto lowerLeft = getParamFromGroup<GlobalPosition>(modelParamGroup, "Grid.LowerLeft");
        const auto upperRight = getParamFromGroup<GlobalPosition>(modelParamGroup, "Grid.UpperRight");

        static const auto couplingPlaneNormal = getParamFromGroup<GlobalPosition>(modelParamGroup,
                                                                                  "Grid.CouplinglineNormal",
                                                                                  [](){ GlobalPosition tmp(0.0); tmp[tmp.size()-1] = 1.0; return tmp; }());

        const auto couplingPlaneNormalDirectionIndex = SnappyGridManagerHelper::directionIndex(couplingPlaneNormal);
        const auto couplingPlane = SnappyGridManagerHelper::makeCouplingPlane(lowerLeft, upperRight, couplingPlaneNormal, modelParamGroup);

        std::cout << "plane: \n";
        for (int i = 0; i < couplingPlane.corners(); ++i)
            std::cout << couplingPlane.corner(i) << std::endl;

        const auto pointsOnAxisParallelLines = SnappyGridManagerHelper::getPointsOnLine(lowerLeft, upperRight, couplingPlaneNormal, otherGrid.leafGridView(), otherData, modelParamGroup_);

        for (int i = 0; i < dim; ++i)
        {
            // set the lower left positions
            positions[i].push_back(lowerLeft[i]);

            if (i != couplingPlaneNormalDirectionIndex) // treat cells parallel to the coupling plane
            {
                if (!pointsOnAxisParallelLines[i].has_value())
                    DUNE_THROW(Dune::GridError, "Something went wrong with the coupling plane normal");

                const auto& pointsOnLine = (*pointsOnAxisParallelLines[i]);

                std::cout << "point " << i << std::endl;
                for (auto x : pointsOnLine)
                    std::cout << x.pos << std::endl;

                // check for user-defined additional points in the upstream area
                const ScalarVector upstreamPositions = getParamFromGroup<ScalarVector>(modelParamGroup, "Grid.UpstreamPositions" + std::to_string(i), ScalarVector{});

                addUpstreamPositions_(i, upstreamPositions, positions, pointsOnLine);
                addUpstreamCells_(i, upstreamPositions, cells);

                addCouplingPositions_(i, positions, pointsOnLine, interFacePositions, lowerLeft, upperRight);
                addCouplingCells_(i, cells, pointsOnLine);

                // check for user-defined additional points in the downstream area
                const ScalarVector downstreamPositions = getParamFromGroup<ScalarVector>(modelParamGroup, "Grid.DownstreamPositions" + std::to_string(i), ScalarVector{});

                addDownstreamPositions_(i, downstreamPositions, positions, upperRight);
                addDownstreamCells_(i, downstreamPositions, cells);

                // set the upper right positions
                positions[i].push_back(upperRight[i]);

                // handle grading of the cells
                // use 1 by default
                grading[i].resize(positions[i].size()-1, 1.0);

                addUpstreamGrading_(i, upstreamPositions, grading);
                addDownstreamGrading_(i, downstreamPositions, positions, grading);
            }
            else // i == couplingPlaneNormalDirectionIndex
            {
                // treat the cells normal to the coupling plane
                if (i != couplingPlaneNormalDirectionIndex)
                    DUNE_THROW(Dune::GridError, "Something went wrong with the coupling plane normal");

                const ScalarVector normalPositions = getParamFromGroup<ScalarVector>(modelParamGroup, "Grid.Positions" + std::to_string(i), ScalarVector{});

                // add positions, cells and grading
                addNormalPositions_(i, normalPositions, positions, upperRight);
                positions[i].push_back(upperRight[i]);
                addNormalCells_(i, normalPositions, cells);
                grading[i].resize(positions[i].size()-1, 1.0);
                addNormalGrading_(i, normalPositions, grading);
            }

            if (positions[i].size() != cells[i].size() + 1)
                DUNE_THROW(Dune::GridError, "Wrong number of cells in " << i);
        }

        // forward to the actual grid creator
        ParentType::init(positions, cells, grading, modelParamGroup);

        gridConstructionData_ = GridConstructionData{std::move(positions), std::move(cells), std::move(grading), std::move(interFacePositions)};
    }

    //! Return data used to create the snappy grid (needed for Dune::TensorProductCoordinates) and the locations of pores intersecting with the interface
    const GridConstructionData& getGridConstructionData() const
    { return gridConstructionData_; }

private:

    /////////////////////////////////////////////////////
    //// Inlet //////////////////////////////////////////
    /////////////////////////////////////////////////////

    template<class PointsOnLine>
    void addUpstreamPositions_(const int directionIndex,
                               const ScalarVector& upstreamPositions,
                               std::array<ScalarVector, dim>& positions,
                               const PointsOnLine& points)
    {
        if (!upstreamPositions.empty())
        {
            const Scalar gridLowerLeft = positions[directionIndex][0];
            for (const Scalar pos : upstreamPositions)
            {
                if ((pos < points[0].pos - points[0].radius) && (pos > gridLowerLeft))
                    positions[directionIndex].push_back(pos);
                else
                    DUNE_THROW(Dune::RangeError, "Make sure to set positions only in the inlet");
            }
        }
    }

    void addUpstreamCells_(const int directionIndex,
                           const ScalarVector& upstreamPositions,
                           std::array<IntVector, dim>& cells)
    {
        const IntVector cellsUpstream = getParamFromGroup<IntVector>(modelParamGroup_,
                                                                     "Grid.UpstreamCells"  + std::to_string(directionIndex),
                                                                     IntVector{});

        if (cellsUpstream.empty())
            return;

        if (cellsUpstream.size() != upstreamPositions.size() + 1)
            DUNE_THROW(Dumux::ParameterException, "UpstreamCells" << std::to_string(directionIndex) << " must equal UpstreamPositions" << std::to_string(directionIndex) <<  " + 1");

        for (int c : cellsUpstream)
        {
            if (c < 1)
                DUNE_THROW(Dumux::ParameterException, "UpstreamCells" << std::to_string(directionIndex) << " must have values greater than zero.");
            cells[directionIndex].push_back(c);
        }
    }

    void addUpstreamGrading_(const int directionIndex,
                             const ScalarVector& upstreamPositions,
                             std::array<ScalarVector, dim>& grading)
    {
        if (upstreamPositions.empty())
            return;

        const ScalarVector upstreamGrading = getParamFromGroup<ScalarVector>(modelParamGroup_, "Grid.UpstreamGrading"  + std::to_string(directionIndex), ScalarVector{});

        if (!upstreamGrading.empty())
        {
            if (upstreamGrading.size() != upstreamPositions.size() + 1)
                DUNE_THROW(Dune::RangeError, "UpstreamGrading"  << std::to_string(directionIndex) << " must equal UpstreamPositions" << std::to_string(directionIndex) <<  " + 1");

            for (int i = 0; i < upstreamPositions.size() + 1; ++i)
                grading[directionIndex][i] = upstreamGrading[i];
        }
    }

    /////////////////////////////////////////////////////
    //// Outlet /////////////////////////////////////////
    /////////////////////////////////////////////////////

    void addDownstreamPositions_(const int directionIndex,
                                 const ScalarVector& downstreamPositions,
                                 std::array<ScalarVector, dim>& gridPositions,
                                 const GlobalPosition& gridUpperRight)
    {
        if (!downstreamPositions.empty())
        {
            for (const Scalar pos : downstreamPositions)
            {
                if ((pos > gridPositions[directionIndex].back()) && (pos < gridUpperRight[directionIndex]))
                    gridPositions[directionIndex].push_back(pos);
                else
                    DUNE_THROW(Dune::RangeError, "Make sure to set ascending positions only in the outlet");
            }
        }
    }

    void addDownstreamCells_(const int directionIndex,
                             const ScalarVector& downstreamPositions,
                             std::array<IntVector, dim>& cells)
    {
        const IntVector downstreamcells = getParamFromGroup<IntVector>(modelParamGroup_,
                                                                       "Grid.DownstreamCells"  + std::to_string(directionIndex),
                                                                       IntVector{});

        if (downstreamcells.empty())
            return;

        if (downstreamcells.size() != downstreamPositions.size() + 1)
            DUNE_THROW(Dumux::ParameterException, "DownstreamCells" << std::to_string(directionIndex) << " must equal DownstreamPositions" << std::to_string(directionIndex) <<  " + 1");

        for (int c : downstreamcells)
        {
            if (c < 1)
                DUNE_THROW(Dumux::ParameterException, "DownstreamCells" << std::to_string(directionIndex) << " must have values greater than zero.");
            cells[directionIndex].push_back(c);
        }
    }

    void addDownstreamGrading_(const int directionIndex,
                               const ScalarVector& downstreamPositions,
                               std::array<ScalarVector, dim>& gridPositions,
                               std::array<ScalarVector, dim>& grading)
    {
        if (downstreamPositions.empty())
            return;

        const ScalarVector downstreamGrading = getParamFromGroup<ScalarVector>(modelParamGroup_, "Grid.DownstreamGrading"  + std::to_string(directionIndex), ScalarVector{});

        if (!downstreamGrading.empty())
        {
            if (downstreamGrading.size() != downstreamPositions.size() + 1)
                DUNE_THROW(Dune::RangeError, "DownstreamGrading"  << std::to_string(directionIndex) << " must equal DownstreamPositions" << std::to_string(directionIndex) <<  " + 1");

            const int offSet = gridPositions[directionIndex].size() - downstreamPositions.size() - 2;
            for (int i = 0; i < downstreamPositions.size() + 1; ++i)
                grading[directionIndex][offSet + i] = downstreamGrading[i];
        }
    }

    /////////////////////////////////////////////////////
    //// Coupling interface//////////////////////////////
    /////////////////////////////////////////////////////

    template<class PointsOnLine>
    void addCouplingPositions_(const int directionIndex,
                               std::array<ScalarVector, dim>& positions,
                               const PointsOnLine& points,
                               std::array<std::optional<ScalarVector>, dim>& interFacePositions,
                               const GlobalPosition& gridLowerLeft,
                               const GlobalPosition& gridUpperRight)
    {
        // create an empty vector for the given directionIndex
        interFacePositions[directionIndex] = ScalarVector{};

        // set the points for the pore body positions
        for (const auto& point : points)
        {
            const auto left = point.pos - point.radius;
            const auto right = point.pos + point.radius;

            if (left < positions[directionIndex].back())
                DUNE_THROW(Dune::RangeError, "Throat radii are too large, they intersect!");

            if (left > gridLowerLeft[directionIndex])
                positions[directionIndex].push_back(left);

            if (right < gridUpperRight[directionIndex])
                positions[directionIndex].push_back(right);

            interFacePositions[directionIndex]->push_back(left);
            interFacePositions[directionIndex]->push_back(right);
            assert(interFacePositions[directionIndex].has_value());
        }
    }

    template<class PointsOnLine>
    void addCouplingCells_(const int directionIndex,
                           std::array<IntVector, dim>& cells,
                           const PointsOnLine& points)
    {
        // set the number of cells above the porous medium
        const int cellsPerThroat = getParamFromGroup<int>(modelParamGroup_, "Grid.CellsPerThroat");
        const int fixedCellsBetweenThroats = getParamFromGroup<int>(modelParamGroup_, "Grid.FixedCellsBetweenThroats", -1);

        // set the cells between the throats
        for (int i = 0; i < points.size(); ++i)
        {
            cells[directionIndex].push_back(cellsPerThroat);
            if (i < points.size() -1)
            {
                if (fixedCellsBetweenThroats > 0)
                    cells[directionIndex].push_back(fixedCellsBetweenThroats);
                else
                {
                    const auto spacingLeft = points[i].radius*2.0 / cellsPerThroat;
                    const auto spacingRight = points[i+1].radius*2.0 / cellsPerThroat;
                    const auto avgSpacing = (spacingLeft + spacingRight) / 2;
                    const auto lengthBetween = (points[i+1].pos - (points[i+1].radius))
                                             - (points[i].pos + (points[i].radius));
                    const int cellsBetween = std::ceil(lengthBetween / avgSpacing);
                    cells[directionIndex].push_back(cellsBetween);
                }
            }
        }
    }

    /////////////////////////////////////////////////////
    //// Normal direction ///////////////////////////////
    /////////////////////////////////////////////////////

    void addNormalPositions_(const int directionIndex,
                             const ScalarVector& normalPositions,
                             std::array<ScalarVector, dim>& gridPositions,
                             const GlobalPosition& gridUpperRight)
    {
        if (!normalPositions.empty())
        {
            for (const Scalar pos : normalPositions)
            {
                if ((pos > gridPositions[directionIndex].back()) && (pos < gridUpperRight[directionIndex]))
                    gridPositions[directionIndex].push_back(pos);
                else
                    DUNE_THROW(Dune::RangeError, "Make sure to set ascending normal positions");
            }
        }
    }

    void addNormalCells_(const int directionIndex,
                         const ScalarVector& normalPositions,
                         std::array<IntVector, dim>& cells)
    {
        const IntVector cellsNormal = getParamFromGroup<IntVector>(modelParamGroup_, "Grid.Cells"  + std::to_string(directionIndex));

        if (cellsNormal.size() != normalPositions.size() + 1)
            DUNE_THROW(Dumux::ParameterException, "Cells" << std::to_string(directionIndex) << " must be of size " << normalPositions.size() + 1);

        for (int c : cellsNormal)
            cells[directionIndex].push_back(c);
    }

    void addNormalGrading_(const int directionIndex,
                           const ScalarVector& normalPositions,
                           std::array<ScalarVector, dim>& grading)
    {
        const ScalarVector normalGrading = getParamFromGroup<ScalarVector>(modelParamGroup_, "Grid.Grading"  + std::to_string(directionIndex), ScalarVector{});

        if (!normalGrading.empty())
        {
            if (normalGrading.size() != normalPositions.size() + 1)
                DUNE_THROW(Dune::RangeError, "Grading"  << std::to_string(directionIndex) << " must be of size " << normalPositions.size() + 1);

            for (int i = 0; i < normalPositions.size() + 1; ++i)
                grading[directionIndex][i] = normalGrading[i];
        }
    }


    std::string modelParamGroup_ = "";
    GridConstructionData gridConstructionData_;
};

} // end namespace Dumux::PoreNetwork

#endif
