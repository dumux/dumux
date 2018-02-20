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
 * \ingroup NavierStokesModel
 * \copydoc Dumux::FluxOverPlane
 */
#ifndef DUMUX_FLUX_OVER_PLANE_STAGGERED_HH
#define DUMUX_FLUX_OVER_PLANE_STAGGERED_HH

#include <numeric>

#include <dune/common/version.hh>
#include <dune/common/fvector.hh>
#include <dune/geometry/type.hh>
#include <dune/geometry/affinegeometry.hh>
#include <dune/geometry/referenceelements.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/geometry/intersectspointgeometry.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesModel
 * \brief  Class used to calculate fluxes over planes. This only works for the staggered grid discretization.
 */
template <class TypeTag>
class FluxOverPlane
{
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using SolutionVector = typename GET_PROP_TYPE(TypeTag, SolutionVector);
    using GridVariables = typename GET_PROP_TYPE(TypeTag, GridVariables);
    using LocalResidual = typename GET_PROP_TYPE(TypeTag, LocalResidual);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using CellCenterPrimaryVariables = typename GET_PROP_TYPE(TypeTag, CellCenterPrimaryVariables);
    using Indices = typename GET_PROP_TYPE(TypeTag, Indices);
    using Element = typename GridView::template Codim<0>::Entity;

    using DofTypeIndices = typename GET_PROP(TypeTag, DofTypeIndices);
    typename DofTypeIndices::CellCenterIdx cellCenterIdx;
    typename DofTypeIndices::FaceIdx faceIdx;

    enum {
        // Grid and world dimension
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    using GlobalPosition = Dune::FieldVector<Scalar, dim>;

    static constexpr auto planeDim = dim - 1;
    using PlaneGeometryType = Dune::AffineGeometry< Scalar, planeDim, dim >;


    /*!
     * \brief Auxiliary class that holds the plane-specific data.
     */
    template<int mydim, int coordDim, class Scalar=double>
    class PlaneData
    {
        using Geo = Dune::AffineGeometry< Scalar, mydim, coordDim >;
    public:

        PlaneData() {}

        PlaneData(Geo&& geo)
        {
            values_.resize(1);
            geometries_.push_back(std::move(geo));
        }

        PlaneData(std::vector<Geo>&& geo)
        {
            values_.resize(geo.size());
            geometries_ = std::forward<decltype(geo)>(geo);
        }

        void addValue(int planeIdx, const CellCenterPrimaryVariables& value)
        {
            values_[planeIdx] += value;
        }

        void addSubPlane(Geo&& geo)
        {
            values_.emplace_back(0.0);
            geometries_.push_back(std::move(geo));
        }

        auto& subPlanes() const
        {
            return geometries_;
        }

        CellCenterPrimaryVariables& value(int planeIdx) const
        {
            return values_[planeIdx];
        }

        auto& values() const
        {
            return values_;
        }

        void printPlaneBoundaries(int planeIdx) const
        {
            const auto& geometry = geometries_[planeIdx];
            for(int i = 0; i < geometry.corners(); ++i)
                std::cout << geometry.corner(i) << "  ";
        }

        void resetValues()
        {
            std::fill(values_.begin(), values_.end(), CellCenterPrimaryVariables(0.0));
        }

    private:

        std::vector<Geo> geometries_;
        std::vector<CellCenterPrimaryVariables> values_;
    };


public:

    using PlaneList = std::vector<PlaneGeometryType>;

    /*!
     * \brief The constructor
     *
     * \param assembler The assembler
     * \param sol The solution vector
     */
    FluxOverPlane(const Problem& problem,
                  const GridVariables& gridVariables,
                  const SolutionVector& sol)
    : problem_(problem)
    , gridVariables_(gridVariables)
    , sol_(sol)
    , localResidual_(LocalResidual(&problem_))
    {
        verbose_  = getParamFromGroup<bool>(GET_PROP_VALUE(TypeTag, ModelParameterGroup), "FluxOverPlane.Verbose", false);
    }

    /*!
     * \brief Add a collection of sub planes under a given name
     *
     * \param name The name of the plane
     * \param planes The list of sub planes
     */
    void addPlane(const std::string& name, PlaneList&& planes )
    {
        planes_[name] = PlaneData<planeDim,dim>(std::forward<decltype(planes)>(planes));
    }

    /*!
     * \brief Add a plane under a given name, specifing the plane's corner points.
     *        This is a specialization for 2D, therefore the plane is actually a line.
     *
     * \param name The name of the plane
     * \param p0 The first corner
     * \param p1 The second corner
     */
    void addPlane(const std::string& name, const GlobalPosition& p0, const GlobalPosition& p1)
    {
        planes_[name].addSubPlane(makePlane(p0, p1));
    }

    /*!
     * \brief Add a plane under a given name, specifing the plane's corner points.
     *        This is a specialization for 3D.
     *
     * \param name The name of the plane
     * \param p0 The first corner
     * \param p1 The second corner
     * \param p2 The third corner
     * \param p3 The fourth corner
     */
    void addPlane(const std::string& name,
                  const GlobalPosition& p0,
                  const GlobalPosition& p1,
                  const GlobalPosition& p2,
                  const GlobalPosition& p3)
    {
        planes_[name].addSubPlane(makePlane(p0, p1, p2, p3));
    }

    /*!
     * \brief Creates a geometrical plane object.
     *        This is a specialization for 2D, therefore the plane is actually a line.
     *
     * \param p0 The first corner
     * \param p1 The second corner
     */
    static PlaneGeometryType makePlane(const GlobalPosition& p0, const GlobalPosition& p1)
    {
        const std::vector< Dune::FieldVector< Scalar, dim > > corners = {p0, p1};
#if DUNE_VERSION_NEWER(DUNE_COMMON,2,6)
        return PlaneGeometryType(Dune::GeometryTypes::line, corners);
#else
        static Dune::GeometryType gt(Dune::GeometryType::simplex, dim-1);
        return PlaneGeometryType(gt, corners);
#endif
    }

    /*!
     * \brief Creates a geometrical plane object.
     *        This is a specialization for 3D.
     *
     * \param p0 The first corner
     * \param p1 The second corner
     * \param p2 The third corner
     * \param p3 The fourth corner
     */
    static PlaneGeometryType makePlane(const GlobalPosition& p0,
                                       const GlobalPosition& p1,
                                       const GlobalPosition& p2,
                                       const GlobalPosition& p3)
    {
        const std::vector< Dune::FieldVector< Scalar, dim > > corners = {p0, p1, p2, p3};
#if DUNE_VERSION_NEWER(DUNE_COMMON,2,6)
        return PlaneGeometryType(Dune::GeometryTypes::quadrilateral, corners);
#else
        static Dune::GeometryType gt(Dune::GeometryType::cube, dim-1);
        return PlaneGeometryType(gt, corners);
#endif
    }

    /*!
     * \brief Calculate the mass or mole fluxes over all planes
     */
    void calculateMassOrMoleFluxes()
    {
        auto fluxType = [this](const auto& problem,
                               const auto& element,
                               const auto& fvGeometry,
                               const auto& elemVolVars,
                               const auto& elemFaceVars,
                               const auto& scvf,
                               const auto& elemFluxVarsCache)
        {
            return localResidual_.computeFluxForCellCenter(problem, element, fvGeometry, elemVolVars, elemFaceVars, scvf, elemFluxVarsCache);
        };

        calculateFluxes(fluxType);
    }

    /*!
     * \brief Calculate the volume fluxes over all planes.
     *        This method simply averages the densities between two adjacent cells.
     */
    void calculateVolumeFluxes()
    {
        const auto isCompositional = std::integral_constant<bool, (GET_PROP_VALUE(TypeTag, NumComponents) > 1) >();
        calculateVolumeFluxesImpl_(isCompositional);
    }

    /*!
     * \brief Calculate the fluxes over all planes for a given flux type.
     *
     * \param fluxType The flux type. This can be a lambda of the following form:
     *                 [](const auto& problem,
                          const auto& element,
                          const auto& fvGeometry,
                          const auto& elemVolVars,
                          const auto& elemFaceVars,
                          const auto& scvf,
                          const auto& elemFluxVarsCache)
                          { return ... ; }
     */
    template<class FluxType>
    void calculateFluxes(const FluxType& fluxType)
    {
        // make sure to reset all the values of the planes, in case this method has been called already before
        for(auto&& plane : planes_)
            plane.second.resetValues();

        // make sure not to iterate over the same dofs twice
        std::vector<bool> dofVisited(problem_.fvGridGeometry().numFaceDofs(), false);

        auto elemVolVars = localView(gridVariables_.curGridVolVars());
        auto elemFluxVarsCache = localView(gridVariables_.gridFluxVarsCache());
        auto elemFaceVars = localView(gridVariables_.curGridFaceVars());

        for(auto&& element : elements(problem_.fvGridGeometry().gridView()))
        {
            auto fvGeometry = localView(problem_.fvGridGeometry());
            fvGeometry.bindElement(element);

            elemVolVars.bind(element, fvGeometry, sol_);
            elemFaceVars.bind(element, fvGeometry, sol_);
            elemFluxVarsCache.bind(element, fvGeometry, elemVolVars);

            for(auto && scvf : scvfs(fvGeometry))
            {
                const auto dofIdx = scvf.dofIndex();
                // do nothing of the dof was already visited
                if(dofVisited[dofIdx])
                    continue;

                dofVisited[dofIdx] = true;

                // iterate through all planes and check if the flux at the given position
                // should be accounted for in the respective plane
                for(auto&& plane : planes_)
                {
                    const auto& subPlanes = plane.second.subPlanes();

                    for(int planeIdx = 0; planeIdx < subPlanes.size(); ++planeIdx)
                    {
                        if(intersectsPointGeometry(scvf.center(), subPlanes[planeIdx]))
                        {
                            const auto result = fluxType(problem_, element, fvGeometry, elemVolVars, elemFaceVars, scvf, elemFluxVarsCache);

                            plane.second.addValue(planeIdx, result);

                            if(verbose_)
                            {
                                std::cout << "Flux at face "  << scvf.center() << " (" << plane.first << "): " << result;
                                std::cout << " (directionIndex: " << scvf.directionIndex() << "; plane boundaries: ";
                                plane.second.printPlaneBoundaries(planeIdx);
                                std::cout << ", planeIdx " << planeIdx << ")" << std::endl;
                            }
                        }
                    }
                }
            }
        }
    }

    /*!
     * \brief Return the fluxes of the individual sub planes of a given name.
     *
     * \param name The name of the plane
     */
    auto& values(const std::string& name) const
    {
        return planes_.at(name).values();
    }

    /*!
     * \brief Return the cumulative net fluxes of a plane of a given name.
     *
     * \param name The name of the plane
     */
    auto netFlux(const std::string& name) const
    {
        const auto& planeResults = values(name);
        return std::accumulate(planeResults.begin(), planeResults.end(), CellCenterPrimaryVariables(0.0));
    }

private:

    /*!
     * \brief Calculate the volume fluxes over all planes for compositional models.
     *        This method simply averages the densities between two adjacent cells.
     */
    void calculateVolumeFluxesImpl_(std::true_type)
    {
        auto fluxType = [this](const auto& problem,
                               const auto& element,
                               const auto& fvGeometry,
                               const auto& elemVolVars,
                               const auto& elemFaceVars,
                               const auto& scvf,
                               const auto& elemFluxVarsCache)
        {
            const auto massOrMoleFlux = localResidual_.computeFluxForCellCenter(problem, element, fvGeometry, elemVolVars, elemFaceVars, scvf, elemFluxVarsCache);

            const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
            const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];

            constexpr bool useMoles = GET_PROP_VALUE(TypeTag, UseMoles);
            const auto density = [useMoles](const auto& volVars)
            {
                return useMoles ? volVars.molarDensity() : volVars.density() ;
            };

            const auto avgDensity = 0.5*density(insideVolVars) + 0.5*density(outsideVolVars);

            constexpr auto replaceCompEqIdx = GET_PROP_VALUE(TypeTag, ReplaceCompEqIdx);
            constexpr auto numComponents = GET_PROP_VALUE(TypeTag, NumComponents);

            const Scalar cumulativeFlux = [replaceCompEqIdx, numComponents, &massOrMoleFlux]()
            {
                Scalar result = 0.0;
                if(replaceCompEqIdx < numComponents)
                    result = massOrMoleFlux[replaceCompEqIdx];
                else
                {
                    for(int i = 0; i < numComponents; ++i)
                        result += massOrMoleFlux[i];
                }

                return result;
            }();

            CellCenterPrimaryVariables tmp(0.0);
            tmp[Indices::conti0EqIdx] = cumulativeFlux / avgDensity;
            return tmp;
        };

        calculateFluxes(fluxType);
    }

    /*!
     * \brief Calculate the volume fluxes over all planes for non-compositional models.
     *        This method simply averages the densities between two adjacent cells.
     */
    void calculateVolumeFluxesImpl_(std::false_type)
    {
        auto fluxType = [this](const auto& problem,
                               const auto& element,
                               const auto& fvGeometry,
                               const auto& elemVolVars,
                               const auto& elemFaceVars,
                               const auto& scvf,
                               const auto& elemFluxVarsCache)
        {
            const Scalar massFlux = localResidual_.computeFluxForCellCenter(problem, element, fvGeometry, elemVolVars, elemFaceVars, scvf, elemFluxVarsCache)[Indices::conti0EqIdx];

            const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
            const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];

            const auto avgDensity = 0.5*insideVolVars.density() + 0.5*outsideVolVars.density();

            CellCenterPrimaryVariables tmp(0.0);
            tmp[Indices::conti0EqIdx] = massFlux / avgDensity;
            return tmp;
        };

        calculateFluxes(fluxType);
    }

    std::map<std::string, PlaneData<planeDim ,dim> > planes_;
    const Problem& problem_;
    const GridVariables& gridVariables_;
    const SolutionVector& sol_;
    bool verbose_;
    LocalResidual localResidual_;
};

} //end namespace

#endif
