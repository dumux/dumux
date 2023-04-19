// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup Discretization
 * \copydoc Dumux::WallDistance
 */
#ifndef DUMUX_DISCRETIZATION_WALL_DISTANCE_HH
#define DUMUX_DISCRETIZATION_WALL_DISTANCE_HH

#include <vector>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/shared_ptr.hh>
#include <dune/common/reservedvector.hh>
#include <dune/grid/common/partitionset.hh>

#include <dumux/common/tag.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/parallel/parallel_for.hh>
#include <dumux/geometry/distancefield.hh>

namespace Dumux {

/*!
 * \ingroup Discretization
 * \brief Class to calculate the wall distance at every element or vertex of a grid
 * \tparam GridGeometry The grid geometry.
 * \tparam DistanceField The type of distance field to use (parameterized with the geometry type)
 */
template<class GridGeometry, template<class> class DistanceField = AABBDistanceField>
class WallDistance
{
    using GridView = typename GridGeometry::GridView;
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using Scalar = typename GridView::Grid::ctype;
    using GlobalPosition = typename SubControlVolumeFace::GlobalPosition;

    static constexpr int dim = GridView::dimension;
    static constexpr int dimWorld = GridView::dimensionworld;
    static_assert(dim == dimWorld, "Wall distances cannot be computed for embedded surface or network domains.");

    using CornerStorage = Dune::ReservedVector<GlobalPosition, (1<<(GridView::dimension-1))>;

    // We use a simplified geometry type here which allows much easier MPI communication
    // for parallel runs compared to the Dune geometry types (due to fixed-size storage).
    // This type only implements a minimal part of the Geometry interface.
    struct SimpleGeometry
    {
        SimpleGeometry() = default;
        SimpleGeometry(CornerStorage&& corners)
        : corners_(std::move(corners))
        , center_(0.0)
        {
            for (int i = 0; i < corners_.size(); ++i)
                center_ += corners_[i];
            center_ /= corners_.size();
        }

        using GlobalCoordinate = GlobalPosition;
        using ctype = typename GlobalCoordinate::value_type;
        static constexpr int coorddimension = GridView::dimensionworld;
        static constexpr int mydimension = GridView::dimension-1;

        std::size_t corners() const
        { return corners_.size(); }

        const auto& corner(int i) const
        { return corners_[i]; }

        const auto& center() const
        { return center_; }

    private:
        CornerStorage corners_;
        GlobalCoordinate center_;
    };

    struct WallDataImpl
    {
        GridIndexType eIdx;
        GridIndexType scvfIdx;
        GlobalPosition scvfOuterNormal;
        int rank; // for parallel runs
    };

    struct ElementCenters {};
    struct VertexCenters {};

public:
    static constexpr auto atElementCenters = Utility::Tag<ElementCenters>{};
    static constexpr auto atVertexCenters = Utility::Tag<VertexCenters>{};
    using WallData = WallDataImpl;

    /*!
     * \brief Constructs a new wall distance object
     * \param gridGeometry the grid geometry of the grid
     * \param tag the policy (either WallDistance::atElementCenters or WallDistance::atVertexCenters)
     * \param select a selection functor taking an scvf and returning bool signifying whether the scvfs is part of the wall
     */
    template<class LocationTag, class ScvfSelectionFunctor>
    WallDistance(std::shared_ptr<const GridGeometry> gridGeometry, LocationTag tag, const ScvfSelectionFunctor& select)
    : gridGeometry_(gridGeometry)
    {
        initializeWallDistance_(select, tag);
    }

    /*!
     * \brief Constructs a new wall distance object
     * \param gridGeometry the grid geometry of the grid
     * \param tag the policy (either WallDistance::atElementCenters or WallDistance::atVertexCenters)
     * \note selects all boundary scvfs as wall faces
     */
    template<class LocationTag>
    WallDistance(std::shared_ptr<const GridGeometry> gridGeometry, LocationTag tag)
    : WallDistance(gridGeometry, tag, [](const FVElementGeometry& fvGeometry, const SubControlVolumeFace& scvf) { return true; }) {}

    //! caller has to make sure the lifetime of grid geometry exceeds the lifetime of wall distance
    template<class LocationTag, class ScvfSelectionFunctor>
    WallDistance(const GridGeometry& gridGeometry, LocationTag tag, const ScvfSelectionFunctor& select)
    : WallDistance(Dune::stackobject_to_shared_ptr(gridGeometry), tag, select) {}

    //! caller has to make sure the lifetime of grid geometry exceeds the lifetime of wall distance
    template<class LocationTag>
    WallDistance(const GridGeometry& gridGeometry, LocationTag tag)
    : WallDistance(Dune::stackobject_to_shared_ptr(gridGeometry), tag) {}

    /*!
     * \brief Returns a vector storing the distance from each DOF location to the nearest wall.
     *        For the atElementCenter policy, this is the distance from the element center to the nearest wall.
     *        For the atVertexCenter policy, this is the distance from the vertex to the nearest wall.
     */
    const std::vector<Scalar>& wallDistance() const
    { return distance_; }

    /*!
     * \brief Returns a vector storing additional information about the nearest scvf on the wall (element index and scvf index).
     *        For the atElementCenter policy, this information is given for each element
     *        For the atVertexCenter policy, this information is given for each vertex.
     */
    const std::vector<WallData>& wallData() const
    { return wallData_; }

private:
    /*!
     * \brief Perform the actual calculation of wall distances.
     * \param considerFace Function object (e.g. a lambda) that determines whether a certain scvf shall be considered
     *                     for the calculation of wall distances.
     * \param loc Location a which the distance to the wall shall be calculated (elementCenter or vertex)
     */
    template<class ConsiderFaceFunction, class LocationTag>
    void initializeWallDistance_(const ConsiderFaceFunction& considerFace, LocationTag loc)
    {
        const std::size_t numSamplingPoints = (loc == atElementCenters)
                                            ? gridGeometry_->gridView().size(0)
                                            : gridGeometry_->gridView().size(dim);
        // Reset the containers.
        wallData_.resize(numSamplingPoints);
        distance_.resize(numSamplingPoints, std::numeric_limits<Scalar>::max());

        std::vector<SimpleGeometry> wallGeometries;
        wallGeometries.reserve(gridGeometry_->numBoundaryScvf());

        std::vector<WallData> tempWallData;
        tempWallData.reserve(gridGeometry_->numBoundaryScvf());

        // Loop over all elements: find all wall scvfs.
        auto fvGeometry = localView(*gridGeometry_);
        for (const auto& element : elements(gridGeometry_->gridView(), Dune::Partitions::interior))
        {
            fvGeometry.bindElement(element);
            if (!fvGeometry.hasBoundaryScvf())
                continue;

            const auto eIdx = gridGeometry_->elementMapper().index(element);

            for (const auto& scvf : scvfs(fvGeometry))
            {
                if (scvf.boundary() && considerFace(fvGeometry, scvf))
                {
                    const auto& geo = fvGeometry.geometry(scvf);
                    CornerStorage corners;
                    for (int i = 0; i < geo.corners(); ++i)
                        corners.push_back(geo.corner(i));

                    wallGeometries.emplace_back(std::move(corners));
                    tempWallData.push_back(WallData{
                        eIdx, scvf.index(), scvf.unitOuterNormal(), gridGeometry_->gridView().comm().rank()
                    });
                }
            }
        }

#if HAVE_MPI
        // Handle parallel runs. We need to prepare a global vector of wall geometries,
        // containing the wall geometries of each process in order to get a correct distance field.
        const bool isParallel = gridGeometry_->gridView().comm().size() > 1;
        std::vector<SimpleGeometry> globalWallGeometries;
        std::vector<WallData> globalTempWallData;
        const auto distanceField = [&]
        {
            if (isParallel)
            {
                const auto& communication = gridGeometry_->gridView().comm();
                const int totalNumberOfBoundaryGeometries = communication.sum(wallGeometries.size());
                globalWallGeometries.resize(totalNumberOfBoundaryGeometries);
                globalTempWallData.resize(totalNumberOfBoundaryGeometries);

                // prepare a displacement vector
                std::vector<int> numGeosPerProcLocal{static_cast<int>(wallGeometries.size())};
                std::vector<int> numGeosPerProcGlobal(communication.size());
                communication.allgather(numGeosPerProcLocal.data(), 1, numGeosPerProcGlobal.data());

                std::vector<int> disp(communication.size(), 0);
                disp[1] = numGeosPerProcGlobal[0];
                for (int i = 2; i < numGeosPerProcGlobal.size(); ++i)
                    disp[i] = disp[i-1] + numGeosPerProcGlobal[i-1];

                // concatenate the wall geometries and temp scvf data of each process into a global vector
                communication.allgatherv(
                    wallGeometries.data(),
                    wallGeometries.size(),
                    globalWallGeometries.data(),
                    numGeosPerProcGlobal.data(),
                    disp.data()
                );

                communication.allgatherv(
                    tempWallData.data(),
                    tempWallData.size(),
                    globalTempWallData.data(),
                    numGeosPerProcGlobal.data(),
                    disp.data()
                );

                // pass the global vector of wall geometries to the distance field
                return DistanceField<SimpleGeometry>(globalWallGeometries);
            }
            else
                return DistanceField<SimpleGeometry>(wallGeometries);
        }();
#else
        const DistanceField<SimpleGeometry> distanceField(wallGeometries);
#endif

        // compute sampling points
        std::vector<GlobalPosition> points(numSamplingPoints);
        if (loc == atElementCenters)
            for (const auto& element : elements(gridGeometry_->gridView()))
                points[gridGeometry_->elementMapper().index(element)] = element.geometry().center();
        else
            for (const auto& vertex : vertices(gridGeometry_->gridView()))
                points[gridGeometry_->vertexMapper().index(vertex)] = vertex.geometry().corner(0);

        // get the actual distances (this is the most expensive part)
        if (loc == atElementCenters)
        {
            const auto kernel = [&](std::size_t eIdx){
                const auto [d, idx] = distanceField.distanceAndIndex(points[eIdx]);
                distance_[eIdx] = d;
#if HAVE_MPI
                wallData_[eIdx] = isParallel ? globalTempWallData[idx] : tempWallData[idx];
#else
                wallData_[eIdx] = tempWallData[idx];
#endif
            };

            runKernel_(numSamplingPoints, kernel);
        }
        else
        {
            const auto kernel = [&](std::size_t vIdx){
                const auto [d, idx] = distanceField.distanceAndIndex(points[vIdx]);
                distance_[vIdx] = d;
#if HAVE_MPI
                wallData_[vIdx] = isParallel ? globalTempWallData[idx] : tempWallData[idx];
#else
                wallData_[vIdx] = tempWallData[idx];
#endif
            };

            runKernel_(numSamplingPoints, kernel);
        }
    }

    template<class Kernel>
    void runKernel_(std::size_t size, const Kernel& kernel)
    {
        // parallelize, if we have enough work (enough evaluation points)
        if (size > 10000)
            Dumux::parallelFor(size, [&](const std::size_t i){ kernel(i); });
        else
            for (std::size_t i = 0; i < size; ++i) kernel(i);
    }

    std::vector<Scalar> distance_;
    std::vector<WallData> wallData_;
    std::shared_ptr<const GridGeometry> gridGeometry_;
};

} // end namespace Dumux

#endif
