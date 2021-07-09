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
 * \ingroup RANSModel
 * \copydoc Dumux::BoundarySearchWallDistance
 */
#ifndef DUMUX_RANS_BOUNDARY_SEARCH_WALL_DISTANCE_HH
#define DUMUX_RANS_BOUNDARY_SEARCH_WALL_DISTANCE_HH

#include <vector>
#include <array>
#include <dune/common/fvector.hh>
#include <dune/common/reservedvector.hh>
#include <dune/geometry/multilineargeometry.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/discretization/evalgradients.hh>
#include <dumux/discretization/evalsolution.hh>
#include <dumux/geometry/distance.hh>

namespace Dumux {

namespace Detail {

template <class ct>
struct TriangleMLGTraits : public Dune::MultiLinearGeometryTraits<ct>
{
    // we use static vectors to store the corners as we know
    // the number of corners in advance (2^(dim-1) corners (1<<(dim-1))
    template<int mydim, int cdim>
    struct CornerStorage
    {
        using Type = std::array<Dune::FieldVector<ct, cdim>, 3>;
    };

    // we know all scvfs will have the same geometry type
    template<int mydim>
    struct hasSingleGeometryType
    {
        static const bool v = true;
        static const unsigned int topologyId = Dune::GeometryTypes::simplex(mydim).id();
    };
};

} // end namespace Detail

/*!
 * \ingroup RANSModel
 * \brief Class to calculate the wall distance based on a simple search algorithm.
 */
template<class GridGeometry>
class BoundarySearchWallDistance
{
    using GridView = typename GridGeometry::GridView;
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using ScvfGeometry = std::decay_t<decltype(std::declval<SubControlVolumeFace>().geometry())>;
    using Scalar = typename GridView::Grid::ctype;
    using GlobalPosition = typename SubControlVolumeFace::GlobalPosition;

    static constexpr auto dim = GridView::dimension;
    static constexpr auto dimWorld = GridView::dimensionworld;

    struct WallScvfData
    {
        GridIndexType eIdx;
        GridIndexType scvfIdx;
    };

    using Triangle = Dune::MultiLinearGeometry<Scalar, dim-1, dimWorld, Detail::TriangleMLGTraits<Scalar>>;
    using Corners = typename Detail::TriangleMLGTraits<Scalar>::template CornerStorage<dim-1, dimWorld>::Type;

    static_assert(dim > 1, "Wall distances cannot be computed for 1D domains.");

public:

    BoundarySearchWallDistance(const GridGeometry& gridGeometry)
    : gridGeometry_(gridGeometry)  {}

    /*!
     * \brief Perform the actual calculation of wall distances. This function overload
     *        considers all scvfs on the boundary.
     */
    void updateWallDistance()
    {
        updateWallDistance([](const auto& scvf) { return true; });
    }

    /*!
     * \brief Perform the actual calculation of wall distances. This function overload
     *        allows to specify which boundary faces shall be considered for computing the distances.
     *
     * \param considerFace Function object (e.g. a lambda) that determines whether a certain scvf shall be considered
     *                     for the calculation of wall distances.
     */
    template<class ConsiderFaceFunction>
    void updateWallDistance(const ConsiderFaceFunction& considerFace)
    {
        // The idea of the algorithm is to iterate over all DOFs (either element centers
        // for cell-centered methods or vertices for the box method) to find the closest distance
        // to the wall. At each DOF, we iterate over all relevant scvfs at the boundary (which we store beforehand)
        // in order to find the minimum wall distance. In 2D, point-segment distances are computed while in 3D,
        // we consider point-triangle distances. In 3D, quadrilateral scvfs are split into two triangles.
        // To speed up computation in 3D, we precompute circumspheres around each wall scvf which allows to
        // exclude certain wall scvfs faster from further more expensive comparisons.

        // Store some temporary info. In 3D, we need more data.
        struct TempScvfInfo2D
        {
            ScvfGeometry geometry;
            GridIndexType eIdx;
            GridIndexType scvfIdx;
        };

        struct TempScvfInfo3D
        {
            std::array<Triangle, 2> geometry;
            GridIndexType eIdx;
            GridIndexType scvfIdx;
            int numTriangles;

            Dune::ReservedVector<GlobalPosition, 2> circumCenter;
            Dune::ReservedVector<Scalar, 2> sphereRadius;
        };

        using TempScvfInfo = std::conditional_t<(dim == 2), TempScvfInfo2D, TempScvfInfo3D>;

        std::vector<TempScvfInfo> tempScvfInfo;
        tempScvfInfo.reserve(gridGeometry_.numBoundaryScvf());

        // Reset the containers.
        wallScvfData_.clear();
        wallScvfData_.resize(gridGeometry_.numDofs());
        distance_.clear();
        distance_.resize(gridGeometry_.numDofs(), std::numeric_limits<Scalar>::max());

        auto fvGeometry = localView(gridGeometry_);

        // First loop over all elements: find all wall scvfs and store circumspheres in 3D.
        for (const auto& element : elements(gridGeometry_.gridView()))
        {
            fvGeometry.bindElement(element);
            if (!fvGeometry.hasBoundaryScvf())
                continue;

            const auto eIdx = gridGeometry_.elementMapper().index(element);

            for (const auto& scvf : scvfs(fvGeometry))
            {
                if (scvf.boundary() && considerFace(scvf))
                {
                    auto geo = scvf.geometry();

                    if constexpr (dim == 3)
                    {
                        Dune::ReservedVector<GlobalPosition, 2> centers;
                        Dune::ReservedVector<Scalar, 2> radii;

                        Corners corners1;
                        corners1[0] = geo.corner(0);
                        corners1[1] = geo.corner(1);
                        corners1[2] = geo.corner(2);
                        auto triangle1 = Triangle(Dune::GeometryTypes::simplex(2), corners1);

                        // Circumvent the problem that Dune::MultiLinearGeometry is not default constructible.
                        // Only the first entry will be used if the scvf itself is a triangle.
                        std::array<Triangle, 2> triangles{triangle1, std::move(triangle1)};
                        int numTriangles = 1;

                        if (geo.corners() == 3)
                        {
                            // Scvf is already triangle, don't add another one.
                            auto [center, radius] = getCircumSphere_(geo);
                            centers.push_back(std::move(center));
                            radii.push_back(std::move(radius));
                        }
                        else
                        {
                            // Scvf has four corners, split into two triangles.
                            numTriangles = 2;

                            // second triangle
                            Corners corners2;
                            corners2[0] = geo.corner(1);
                            corners2[1] = geo.corner(2);
                            corners2[2] = geo.corner(3);
                            triangles[1] = Triangle(Dune::GeometryTypes::simplex(2), corners2);

                            auto [center1, radius1] = getCircumSphere_(triangles[0]);
                            centers.push_back(std::move(center1));
                            radii.push_back(std::move(radius1));

                            auto [center2, radius2] = getCircumSphere_(triangles[1]);
                            centers.push_back(std::move(center2));
                            radii.push_back(std::move(radius2));
                        }

                        tempScvfInfo.push_back({std::move(triangles), eIdx, scvf.index(), numTriangles, std::move(centers), std::move(radii)});
                    }
                    else
                        tempScvfInfo.push_back({std::move(geo), eIdx, scvf.index()});
                }
            }
        }

        // Second loop over all elements: find distances of all dofs to wall.
        for (const auto& element : elements(gridGeometry_.gridView()))
        {
            fvGeometry.bindElement(element);

            for (const auto& scv : scvs(fvGeometry))
            {
                for (const auto& wallScvfInfo : tempScvfInfo)
                {
                    if constexpr (dim == 3) // 3D
                    {
                        // Computing the distance to the circumsphere is much cheaper than computing the distance to the actual triangle.
                        // Use this to early exclude candidates.
                        Dune::ReservedVector<bool, 2> skip;
                        for (int i = 0; i < wallScvfInfo.numTriangles; ++i)
                        {
                            if (distanceToSphere_(scv.dofPosition(), wallScvfInfo.circumCenter[i], wallScvfInfo.sphereRadius[i]) > distance_[scv.dofIndex()])
                                skip.push_back(true);
                            else
                                skip.push_back(false);
                        }

                        if (std::all_of(skip.begin(), skip.end(), [](const bool s) { return s; }))
                            continue;

                        // Compute the actual distance to the triangle. Remember that there are two triangles for quadrilateral scvfs.
                        for (int i = 0; i < wallScvfInfo.numTriangles; ++i)
                        {
                            const auto d = getDistance_(scv.dofPosition(), wallScvfInfo.geometry[i]);
                            if (d < distance_[scv.dofIndex()])
                            {
                                distance_[scv.dofIndex()] = d;
                                wallScvfData_[scv.dofIndex()] = WallScvfData{wallScvfInfo.eIdx, wallScvfInfo.scvfIdx};
                            }
                        }
                    }
                    else // 2D
                    {
                        const auto d = getDistance_(scv.dofPosition(), wallScvfInfo.geometry);
                        if (d < distance_[scv.dofIndex()])
                        {
                            distance_[scv.dofIndex()] = d;
                            wallScvfData_[scv.dofIndex()] = WallScvfData{wallScvfInfo.eIdx, wallScvfInfo.scvfIdx};
                        }
                    }
                }
            }
        }
    }

    /*!
     * \brief Returns a vector storing the distance from each DOF to the nearest wall.
     *        For cell-centered schemes, this is the distance from the element center to the nearest wall.
     *        For node-centered schemes (box), this is the distance from the vertex to the nearest wall.
     */
    const std::vector<Scalar>& wallDinstance() const
    { return distance_; }

    /*!
     * \brief Returns a vector storing additional information about the nearest scvf on the wall (element index and scvf index).
     *        For cell-centered schemes, this information is given for each element
     *        For node-centered schemes (box), this information is given for each vertex.
     */
    const std::vector<WallScvfData>& wallScvfData() const
    { return wallScvfData_; }

    /*!
     * \brief Clear the data storage to free some memory.
     */
    void clear()
    {
        distance_.clear();
        wallScvfData_.clear();
    }

private:

    template<class Geometry>
    Scalar getDistance_(const GlobalPosition& point, const Geometry& geo) const
    {
        if constexpr (dim == 2)
            return distancePointSegment(point, geo);
        else
            return distancePointTriangle(point, geo);
    }

    Scalar distanceToSphere_(const GlobalPosition& point, const GlobalPosition& sphereCenter, const Scalar sphereRadius) const
    {
        return (point - sphereCenter).two_norm() - sphereRadius;
    }

    // see https://gamedev.stackexchange.com/a/60631 and https://www.ics.uci.edu/~eppstein/junkyard/circumcenter.html
    template<class Geometry>
    auto getCircumSphere_(const Geometry& geo)
    {
        const auto& a = geo.corner(0);
        const auto& b = geo.corner(1);
        const auto& c = geo.corner(2);

        const auto ac = c - a;
        const auto ab = b - a;
        const auto n = crossProduct(ab, ac);

        const auto distCenterToA = (crossProduct(n, ab)*ac.two_norm2() + crossProduct(ac, n)*ab.two_norm2()) / (2.0*n.two_norm2());
        const auto radius = distCenterToA.two_norm();
        const auto circumCenter = a + distCenterToA;

        return std::make_pair(circumCenter, radius);
    }

    std::vector<Scalar> distance_;
    std::vector<WallScvfData> wallScvfData_;
    const GridGeometry& gridGeometry_;
};

} // end namespace Dumux

#endif
