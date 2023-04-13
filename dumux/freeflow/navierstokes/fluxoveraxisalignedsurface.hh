// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesModel
 * \copydoc Dumux::FluxOverAxisAlignedSurface
 */
#ifndef DUMUX_FREELOW_NAVIERSTOKES_FLUX_OVER_AXISALIGNED_SURFACE_HH
#define DUMUX_FREELOW_NAVIERSTOKES_FLUX_OVER_AXISALIGNED_SURFACE_HH

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
 * \brief  Class used to calculate fluxes over axis-aligned surfaces.
 */
template<class GridVariables, class SolutionVector, class LocalResidual>
class FluxOverAxisAlignedSurface
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

    static_assert(dim > 1, "Only implemented for dim > 1");

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    // in 2D, the surface is represented as a line
    using SurfaceT = Dune::AxisAlignedCubeGeometry<Scalar, (dim == 2 ? 1 : 2), dimWorld>;

    struct SurfaceData
    {
        SurfaceT surface;
        std::size_t normalDirectionIndex;
        NumEqVector flux;
    };

public:

    using Surface = SurfaceT;

    /*!
     * \brief The constructor
     */
    FluxOverAxisAlignedSurface(const GridVariables& gridVariables,
                               const SolutionVector& sol,
                               const LocalResidual& localResidual)
    : gridVariables_(gridVariables)
    , sol_(sol)
    , localResidual_(localResidual)
    {
        verbose_ = getParamFromGroup<bool>(problem_().paramGroup(), "FluxOverAxisAlignedSurface.Verbose", false);
    }

    /*!
     * \brief Add an axis-aligned surface with a given name
     *
     * \param name The name of the surface
     * \param surface The surface to add
     */
    template<class T>
    void addAxisAlignedSurface(const std::string& name, T&& surface)
    {
        static_assert(std::is_same_v<std::decay_t<T>, Surface>);
        surfaces_.emplace(std::make_pair(
            name, std::make_pair(surface, NumEqVector(0.0))
        ));
    }

    /*!
     * \brief Add an axis-aligned surface (segment in 2D) with a given name, specifying the surface's corner points.
     *
     * \param name The name of the surface
     * \param lowerLeft Lower left corner of surface
     * \param upperRight Upper right corner of surface
     */
     void addAxisAlignedSurface(const std::string& name,
                                const GlobalPosition& lowerLeft,
                                const GlobalPosition& upperRight)
    {
        using std::abs;
        const GlobalPosition v = upperRight - lowerLeft;
        const auto it = std::find_if(v.begin(), v.end(), [](const auto& x){ return abs(x) < 1e-20; });
        if (it == v.end())
            DUNE_THROW(Dune::InvalidStateException, "Surface is not axis-parallel!");

        const std::size_t normalDirectionIndex = std::distance(v.begin(), it);
        auto inSurfaceAxes = std::move(std::bitset<dimWorld>{}.set());
        inSurfaceAxes.set(normalDirectionIndex, false);
        auto surface = Surface(lowerLeft, upperRight, inSurfaceAxes);

        surfaces_.emplace(std::make_pair(
            name,
            SurfaceData{
                std::move(surface), normalDirectionIndex, NumEqVector(0.0)
            }
        ));
    }

    /*!
     * \brief Add an axis-aligned plane (line in 2D) with a given name, specifying the planes's center and normal.
     *
     * \param name The name of the plane
     * \param center Center point of the plane
     * \param normalDirectionIndex Index of the plane's normal axis (0=x, 1=y, 2=z)
     */
     void addAxisAlignedPlane(const std::string& name,
                              const GlobalPosition& center,
                              const std::size_t normalDirectionIndex)
    {
        GlobalPosition lowerLeft = gridVariables_.gridGeometry().bBoxMin();
        GlobalPosition upperRight = gridVariables_.gridGeometry().bBoxMax();

        lowerLeft[normalDirectionIndex] = center[normalDirectionIndex];
        upperRight[normalDirectionIndex] = center[normalDirectionIndex];

        auto inSurfaceAxes = std::move(std::bitset<dimWorld>{}.set());
        inSurfaceAxes.set(normalDirectionIndex, false);
        auto surface = Surface(lowerLeft, upperRight, inSurfaceAxes);

        surfaces_.emplace(std::make_pair(
            name,
            SurfaceData{
                std::move(surface), normalDirectionIndex, NumEqVector(0.0)
            }
        ));
    }

    /*!
     * \brief Calculate the fluxes over all surfaces.
     */
    void calculateAllFluxes()
    {
        auto fluxType = [this](const auto& element,
                               const auto& fvGeometry,
                               const auto& elemVolVars,
                               const auto& scvf,
                               const auto& elemFluxVarsCache)
        {
            return localResidual_.evalFlux(
                problem_(), element, fvGeometry, elemVolVars, elemFluxVarsCache, scvf
            );
        };

        calculateFluxes(fluxType);
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
    void calculateFluxes(const FluxType& fluxType)
    {
        // make sure to reset all the values of the surfaces, in case this method has been called already before
        for (auto& surface : surfaces_)
            surface.second.flux = 0.0;

        snapSurfaceToClosestFace_();
        calculateFluxes_(fluxType);
    }

    /*!
     * \brief Return the flux over given surface
     *
     * \param name The name of the surface
     */
    const auto& flux(const std::string& name) const
    {
        return surfaces_.at(name).flux;
    }

    /*!
     * \brief Provides access to all surfaces.
     */
    const std::map<std::string, SurfaceData>& surfaces() const
    { return surfaces_; }

    /*!
     * \brief Prints all fluxes.
     */
    void printAllFluxes() const
    {
        for (const auto& [name, data] : surfaces_)
            std::cout << "Flux over surface " << name << ": " << data.flux << std::endl;
    }

private:

    template<class FluxType>
    void calculateFluxes_(const FluxType& fluxType)
    {
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
                for (auto& [name, surfaceData] : surfaces_)
                {
                    if (considerScvf_(scvf, surfaceData))
                    {
                        const auto result = fluxType(element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);
                        surfaceData.flux += result;

                        if (verbose_)
                            std::cout << "At element " << problem_().gridGeometry().elementMapper().index(element)
                                      << ": Flux at face " << scvf.ipGlobal() << ": " << result << " (" << name << ")" << std::endl;
                    }
                }
            }
        }
    }

    //! Check whether a scvf should be considered for the flux calculation
    bool considerScvf_(const SubControlVolumeFace& scvf, const SurfaceData& SurfaceData) const
    {
        // In order to avoid considering scvfs at the same element intersection (and hence, the corresponding flux) twice,
        // only use those with a unit outer normal pointing towards positive coordinate direction,
        // unless the scvf lies on a boundary (then there is no second scvf).
        if (scvf.boundary() || !std::signbit(scvf.unitOuterNormal()[SurfaceData.normalDirectionIndex]))
            return intersectsPointGeometry(scvf.ipGlobal(), SurfaceData.surface);
        else
            return false;
    }

    void snapSurfaceToClosestFace_()
    {
        using GeometriesEntitySet = Dumux::GeometriesEntitySet<Surface>;
        const auto gridView = problem_().gridGeometry().gridView();

        for (auto& [name, surfaceData] : surfaces_)
        {
            GeometriesEntitySet entitySet({surfaceData.surface});
            Dumux::BoundingBoxTree<GeometriesEntitySet> geometriesTree(std::make_shared<GeometriesEntitySet>(entitySet));
            const auto intersectingElements = intersectingEntities(
                problem_().gridGeometry().boundingBoxTree(), geometriesTree
            );

            if (intersectingElements.empty())
            {
                std::cout << "surface boundaries: " << std::endl;
                printSurfaceBoundaries_(surfaceData.surface);

                DUNE_THROW(Dune::InvalidStateException, "surface " << name << " does not intersect with any element");
            }

            std::vector<std::size_t> sortedResults;
            sortedResults.reserve(gridView.size(0));

            for (const auto& i : intersectingElements)
                sortedResults.push_back(i.first());

            std::sort(sortedResults.begin(), sortedResults.end());
            sortedResults.erase(std::unique(
                sortedResults.begin(), sortedResults.end()
            ), sortedResults.end());

            // pick the first intersecting element and make sure the surface snaps to the closest face with the same (or opposite facing) normal vector
            GlobalPosition normalVector(0.0);
            normalVector[surfaceData.normalDirectionIndex] = 1.0;

            const auto& firstIntersectingElement = problem_().gridGeometry().element(sortedResults[0]);
            Scalar distance = std::numeric_limits<Scalar>::max();
            bool snappingOcurred = false;

            GlobalPosition surfaceLowerLeft = surfaceData.surface.corner(0);
            GlobalPosition surfaceUpperRight = surfaceData.surface.corner(3);

            bool surfaceAlreadyOnFaces = false;
            for (const auto& intersection : intersections(gridView, firstIntersectingElement))
            {
                if (surfaceAlreadyOnFaces)
                    continue;

                using std::abs;
                if (abs(1.0 - abs(normalVector * intersection.centerUnitOuterNormal())) < 1e-8)
                {

                    const auto getDistance = [](const auto& p, const auto& geo)
                    {
                        if constexpr (dim == 2)
                            return distancePointSegment(p, geo);
                        else
                            return distancePointPolygon(p, geo);
                    };

                    const auto& geo = intersection.geometry();
                    if (const Scalar d = getDistance(geo.center(), surfaceData.surface); d < 1e-8 * diameter(geo))
                    {
                        // no snapping required, face already lies on surface
                        surfaceAlreadyOnFaces = true;
                        snappingOcurred = false;
                    }
                    else if (d < distance)
                    {
                        distance = d;
                        snappingOcurred = true;

                        // move the surface boundaries
                        for (int i = 0; i < surfaceData.surface.corners(); ++i)
                        {
                            const auto& faceCenter = geo.center();
                            surfaceLowerLeft[surfaceData.normalDirectionIndex] = faceCenter[surfaceData.normalDirectionIndex];
                            surfaceUpperRight[surfaceData.normalDirectionIndex] = faceCenter[surfaceData.normalDirectionIndex];
                        }
                    }
                }
            }

            if (snappingOcurred)
            {
                std::cout << "\n\nSurface '" << name <<  "' was automatically snapped to the closest faces" << std::endl;
                std::cout << "Old surface boundaries: " << std::endl;
                printSurfaceBoundaries_(surfaceData.surface);

                // overwrite the old surface with the new boundaries
                auto inSurfaceAxes = std::move(std::bitset<dimWorld>{}.set());
                inSurfaceAxes.set(surfaceData.normalDirectionIndex, false);
                surfaceData.surface = Surface{surfaceLowerLeft, surfaceUpperRight, inSurfaceAxes};

                std::cout << "New surface boundaries: " << std::endl;
                printSurfaceBoundaries_(surfaceData.surface);
                std::cout << std::endl;
            }
        }
    }

    void printSurfaceBoundaries_(const Surface& surface) const
    {
        for (int i = 0; i < surface.corners(); ++i)
            std::cout << surface.corner(i) << std::endl;
    }

    const auto& problem_() const { return gridVariables_.curGridVolVars().problem(); }

    std::map<std::string, SurfaceData> surfaces_;
    const GridVariables& gridVariables_;
    const SolutionVector& sol_;
    const LocalResidual localResidual_; // store a copy of the local residual
    bool verbose_;
};

} // end namespace Dumux

#endif
