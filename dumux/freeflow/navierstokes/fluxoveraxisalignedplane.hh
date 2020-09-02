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
 * \ingroup NavierStokesModel
 * \copydoc Dumux::FluxOverSurface
 */
#ifndef DUMUX_FLUX_OVER_SURFACE_STAGGERED_HH
#define DUMUX_FLUX_OVER_SURFACE_STAGGERED_HH

#include <algorithm>
#include <type_traits>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/geometry/axisalignedcubegeometry.hh>

#include <dumux/common/parameters.hh>
#include <dumux/geometry/diameter.hh>
#include <dumux/geometry/distance.hh>
#include <dumux/geometry/intersectspointgeometry.hh>
#include <dumux/geometry/geometricentityset.hh>
#include <dumux/geometry/intersectingentities.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesModel
 * \brief  Class used to calculate fluxes over surfaces. This only works for the staggered grid discretization.
 */
template<class GridVariables, class SolutionVector, class LocalResidual>
class FluxOverAxisAlignedPlane
{
    using Scalar = typename GridVariables::Scalar;
    using GridGeometry = typename GridVariables::GridGeometry;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GridView = typename GridGeometry::GridView;
    using VolumeVariables = typename GridVariables::VolumeVariables;
    using Element = typename GridView::template Codim<0>::Entity;
    using NumEqVector = typename LocalResidual::ElementResidualVector::value_type;

    static constexpr auto dim = GridView::dimension;
    static constexpr auto dimWorld = GridView::dimensionworld;

    static_assert(dim > 1);

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    // in 2D, the plane is represeted as a line
    using PlaneT = Dune::AxisAlignedCubeGeometry<Scalar, (dim == 2 ? 1 : 2), dimWorld>;

    struct PlaneData
    {
        PlaneT plane;
        std::size_t normalDirectionIndex;
        NumEqVector flux;
    };

public:

    using Plane = PlaneT;

    enum class EvaluationMode
    {
        snapToClosestFace, average
    };

    /*!
     * \brief The constructor
     */
    FluxOverAxisAlignedPlane(const GridVariables& gridVariables,
                             const SolutionVector& sol)
    : gridVariables_(gridVariables),
      sol_(sol)
    {
        verbose_  = getParamFromGroup<bool>(problem_().paramGroup(), "FluxOverPlane.Verbose", false);
    }

    /*!
     * \brief Add a collection of sub surfaces under a given name
     *
     * \param name The name of the surface
     * \param surfaces The list of sub surfaces
     */
    template<class T>
    void addPlane(const std::string& name, T&& plane)
    {
        static_assert(std::is_same_v<std::decay_t<T>, Plane>);
        planes_.emplace(std::make_pair(name, std::make_pair(plane, NumEqVector{})));
    }

    /*!
     * \brief Add a surface under a given name, specifying the surface's corner points.
     *        This is a specialization for 2D, therefore the surface is actually a line.
     *
     * \param name The name of the surface
     * \param p0 The first corner
     * \param p1 The second corner
     */
     void addPlane(const std::string& name,
                   const GlobalPosition& lowerLeft,
                   const GlobalPosition& upperRight,
                   const std::size_t normalDirectionIndex,
                   const bool performCheck = true)
    {
        if constexpr (dimWorld == 2)
            addPlane2D_(name, lowerLeft, upperRight, normalDirectionIndex, performCheck);
        else
            addPlane3D_(name, lowerLeft, upperRight, normalDirectionIndex, performCheck);
    }

    void calculateAllScalarFluxes(const EvaluationMode evalMode = EvaluationMode::snapToClosestFace)
    {
        auto fluxType = [this](const auto& element,
                               const auto& fvGeometry,
                               const auto& elemVolVars,
                               const auto& scvf,
                               const auto& elemFluxVarsCache)
        {
            LocalResidual localResidual(&problem_());
            return localResidual.evalFlux(problem_(), element, fvGeometry, elemVolVars, elemFluxVarsCache, scvf);
        };

        calculateFluxes(fluxType, evalMode);
    }

    /*!
     * \brief Calculate the fluxes over all surfaces for a given flux type.
     *
     * \param fluxType The flux type. This can be a lambda of the following form:
     *                 [](const auto& element,
                          const auto& fvGeometry,
                          const auto& elemVolVars,
                          const auto& scvf,
                          const auto& elemFluxVarsCache)
                          { return ... ; }
     */
    template<class FluxType>
    void calculateFluxes(const FluxType& fluxType, const EvaluationMode evalMode = EvaluationMode::snapToClosestFace)
    {
        // make sure to reset all the values of the planes, in case this method has been called already before
        for (auto& plane : planes_)
            plane.second.flux = 0.0;

        if (evalMode == EvaluationMode::snapToClosestFace)
            calculateFluxesWithSnapping_(fluxType);
        else
            DUNE_THROW(Dune::NotImplemented, "Averaging of fluxes not implemented yet");
    }

    /*!
     * \brief Return the fluxes of the individual sub surface of a given name.
     *
     * \param name The name of the surface
     */
    const auto& netFlux(const std::string& name) const
    {
        return planes_.at(name).flux;
    }

private:

    /*!
     * \brief Add a surface under a given name, specifying the surface's corner points.
     *        This is a specialization for 2D, therefore the surface is actually a line.
     *
     * \param name The name of the surface
     * \param p0 The first corner
     * \param p1 The second corner
     */
     void addPlane2D_(const std::string& name,
                      const GlobalPosition& lowerLeft,
                      const GlobalPosition& upperRight,
                      const std::size_t normalDirectionIndex,
                      const bool performCheck)
    {
        auto inPlaneAxes = std::move(std::bitset<dimWorld>{}.set());
        inPlaneAxes.set(normalDirectionIndex, false);
        auto plane = Plane(lowerLeft, upperRight, inPlaneAxes);

        if (performCheck)
        {
            auto v01 = upperRight - lowerLeft;
            v01 /= v01.two_norm();
            GlobalPosition testVector(0.0);
            testVector[normalDirectionIndex] = 1.0;

            using std::abs;
            if (abs(v01 * testVector) > 1e-8)
                DUNE_THROW(Dune::InvalidStateException, "Points do not span a axis-parallel plane in 3D");
        }

        planes_.emplace(std::make_pair(name, PlaneData{std::move(plane), normalDirectionIndex, NumEqVector{}}));
    }

    /*!
     * \brief Add a surface under a given name, specifying the surface's corner points.
     *        This is a specialization for 2D, therefore the surface is actually a line.
     *
     * \param name The name of the surface
     * \param p0 The first corner
     * \param p1 The second corner
     */
    void addPlane3D_(const std::string& name,
                     const GlobalPosition& lowerLeft,
                     const GlobalPosition& upperRight,
                     const std::size_t normalDirectionIndex,
                     const bool performCheck)
    {
        static constexpr std::array<std::size_t, 3> dirIndices = {0, 1, 2};
        if (std::none_of(dirIndices.begin(), dirIndices.end(), [&](auto i) { return normalDirectionIndex == i;} ))
            DUNE_THROW(Dune::InvalidStateException, normalDirectionIndex << " is not a valid index. Use 0, 1 or 2.");

        auto inPlaneAxes = std::move(std::bitset<dimWorld>{}.set());
        inPlaneAxes.set(normalDirectionIndex, false);
        auto plane = Plane(lowerLeft, upperRight, inPlaneAxes);

        if (performCheck)
        {
            const auto v01 = plane.corner(1) - plane.corner(0);
            const auto v02 = plane.corner(2) - plane.corner(0);
            auto planeNormal = crossProduct(v01, v02);
            planeNormal /= planeNormal.two_norm();
            GlobalPosition testVector(0.0);
            testVector[normalDirectionIndex] = 1.0;

            using std::abs;
            if (abs(1.0 - abs(planeNormal * testVector)) > 1e-8 * v01.two_norm())
                DUNE_THROW(Dune::InvalidStateException, "Points do not span a axis-parallel plane in 3D");
        }

        planes_.emplace(std::make_pair(name, PlaneData{std::move(plane), normalDirectionIndex, NumEqVector{}}));
    }

    template<class FluxType>
    void calculateFluxesWithSnapping_(const FluxType& fluxType)
    {
        snapPlaneToClosestFace_();

        auto fvGeometry = localView(problem_().gridGeometry());
        auto elemVolVars = localView(gridVariables_.curGridVolVars());
        auto elemFluxVarsCache = localView(gridVariables_.gridFluxVarsCache());

        for (const auto& element : elements(problem_().gridGeometry().gridView()))
        {
            fvGeometry.bindElement(element);
            elemVolVars.bindElement(element, fvGeometry, sol_);
            elemFluxVarsCache.bindElement(element, fvGeometry, elemVolVars);

            for (const auto& scvf : scvfs(fvGeometry))
            {
                // iterate through all surfaces and check if the flux at the given position
                // should be accounted for in the respective surface
                for (auto& [name, planeData] : planes_)
                {
                    if (considerScvf_(scvf, planeData))
                    {
                        const auto result = fluxType(element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);
                        planeData.flux += result;

                        if (verbose_)
                            std::cout << "At element " << problem_().gridGeometry().elementMapper().index(element)
                                      << ": Flux at face " << scvf.ipGlobal() << ": " << result << " (" << name << ")" << std::endl;
                    }
                }
            }
        }
    }

    //! Check whether a scvf should be considered for the flux calculation
    bool considerScvf_(const SubControlVolumeFace& scvf, const PlaneData& planeData) const
    {
        // In order to avoid considering scvfs at the same element intersection (and hence, the corresponding flux) twice,
        // only use those with a unit outer normal pointing towards positive coordinate direction,
        // unless the scvf lies on a boundary (then there is no second scvf).
        if (scvf.boundary() || !std::signbit(scvf.unitOuterNormal()[planeData.normalDirectionIndex]))
            return intersectsPointGeometry(scvf.ipGlobal(), planeData.plane);
        else
            return false;
    }

    void snapPlaneToClosestFace_()
    {
        using GeometriesEntitySet = Dumux::GeometriesEntitySet<Plane>;
        const auto gridView = problem_().gridGeometry().gridView();

        for (auto& [name, planeData] : planes_)
        {
            GeometriesEntitySet entitySet({planeData.plane});
            Dumux::BoundingBoxTree<GeometriesEntitySet> geometriesTree(std::make_shared<GeometriesEntitySet>(entitySet));
            const auto intersectingElements = intersectingEntities(problem_().gridGeometry().boundingBoxTree(), geometriesTree);

            if (intersectingElements.empty())
            {
                std::cout << "Plane boundaries: " << std::endl;
                printPlaneBoundaries_(planeData.plane);

                DUNE_THROW(Dune::InvalidStateException, "Plane " << name << " does not intersect with any element");
            }

            std::vector<std::size_t> sortedResults;
            sortedResults.reserve(gridView.size(0));

            for (const auto& i : intersectingElements)
                sortedResults.push_back(i.first());

            std::sort(sortedResults.begin(), sortedResults.end());
            sortedResults.erase(std::unique(sortedResults.begin(), sortedResults.end()), sortedResults.end());

            // pick the first intersecting element and make sure the plane snaps to the closest face with the same (or opposite facing) normal vector
            GlobalPosition normalVector(0.0);
            normalVector[planeData.normalDirectionIndex] = 1.0;

            const auto& firstIntersectingElement = problem_().gridGeometry().element(sortedResults[0]);
            Scalar distance = std::numeric_limits<Scalar>::max();
            bool snappingOccured = false;

            GlobalPosition planeLowerLeft = planeData.plane.corner(0);
            GlobalPosition planeUpperRight = planeData.plane.corner(3);

            bool planeAlreadyOnFaces = false;
            for (const auto& intersection : intersections(gridView, firstIntersectingElement))
            {
                if (planeAlreadyOnFaces)
                    continue;

                using std::abs;
                if (abs(1.0 - abs(normalVector * intersection.centerUnitOuterNormal())) < 1e-8)
                {
                    const auto& geo = intersection.geometry();
                    const auto& faceCenter = geo.center();
                    const auto d = pointToPlaneDistance_(intersection.geometry().center(), planeData.plane);

                    // no snapping required, face already lies on plane
                    if (d < 1e-8 * diameter(geo))
                    {
                        planeAlreadyOnFaces = true;
                        snappingOccured = false;
                        continue;
                    }

                    if (d < distance)
                    {
                        distance = d;
                        snappingOccured = true;

                        // move the plane boundaries
                        for (int i = 0; i < planeData.plane.corners(); ++i)
                        {
                            planeLowerLeft[planeData.normalDirectionIndex] = faceCenter[planeData.normalDirectionIndex];
                            planeUpperRight[planeData.normalDirectionIndex] = faceCenter[planeData.normalDirectionIndex];
                        }
                    }
                }
            }

            if (snappingOccured)
            {
                std::cout << "\n\nPlane '" << name <<  "' was automatically snapped to the closest faces" << std::endl;
                std::cout << "Old plane boundaries: " << std::endl;
                printPlaneBoundaries_(planeData.plane);

                // overwrite the old plane with the new boundaries
                auto inPlaneAxes = std::move(std::bitset<dimWorld>{}.set());
                inPlaneAxes.set(planeData.normalDirectionIndex, false);
                planeData.plane = Plane{planeLowerLeft, planeUpperRight, inPlaneAxes};

                std::cout << "New plane boundaries: " << std::endl;
                printPlaneBoundaries_(planeData.plane);
                std::cout << std::endl;
            }
        }
    }

    auto getPlaneCoefficients_(const Plane& plane)
    {
        struct Result {Scalar a; Scalar b; Scalar c; Scalar d;} result;
        const auto& p0 = plane.corner(0);
        const auto& p1 = plane.corner(1);
        const auto& p2 = plane.corner(2);

        const auto& v01 = p1 -p0;
        const auto& v02 = p2 -p0;
        const auto c = crossProduct(v01, v02);

        result.a = c[0];
        result.b = c[1];
        result.c = c[2];
        result.d = -(result.a*p0[0] + result.b*p0[1] + result.c*p0[2]);

        return result;
    }

    Scalar pointToPlaneDistance_(const GlobalPosition& p0, const Plane& plane)
    {
        if constexpr (dimWorld == 2)
            return distancePointLine(p0, plane);
        else
        {
            const auto [a, b, c, d] = getPlaneCoefficients_(plane);
            using std::abs; using std::sqrt;
            return abs(a*p0[0] + b*p0[1] + c*p0[2] + d) / sqrt(a*a + b*b + c*c);
        }
    }

    void printPlaneBoundaries_(const Plane& plane) const
    {
        for (int i = 0; i < plane.corners(); ++i)
            std::cout << plane.corner(i) << std::endl;
    }

    const auto& problem_() const { return gridVariables_.curGridVolVars().problem(); }

    std::map<std::string, PlaneData> planes_;
    const GridVariables& gridVariables_;
    const SolutionVector& sol_;
    bool verbose_;
};

} //end namespace

#endif
