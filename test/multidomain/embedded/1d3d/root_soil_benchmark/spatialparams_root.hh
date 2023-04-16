// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup EmbeddedTests
 * \brief The spatial parameters class for the root problem
 */
#ifndef DUMUX_TEST_ROOT_SOIL_BENCHMARK_ROOT_SPATIALPARAMS_HH
#define DUMUX_TEST_ROOT_SOIL_BENCHMARK_ROOT_SPATIALPARAMS_HH

#include <vector>
#include <utility>
#include <memory>

#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>
#include <dumux/io/grid/griddata.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/porousmediumflow/fvspatialparams1p.hh>

namespace Dumux {

/*!
 * \ingroup EmbeddedTests
 * \brief The spatial parameters class for the root problem
 */
template<class GridGeometry, class Scalar>
class RootSpatialParams
: public FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar, RootSpatialParams<GridGeometry, Scalar>>
{
    using ThisType = RootSpatialParams<GridGeometry, Scalar>;
    using ParentType = FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar, ThisType>;
    using Grid = typename GridGeometry::Grid;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    // export permeability type
    using PermeabilityType = Scalar;

    RootSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry,
                      std::shared_ptr<const GridData<Grid>> gridData)
    : ParentType(gridGeometry), gridData_(gridData)
    {
        porosity_ = getParam<Scalar>("Root.SpatialParams.Porosity", 0.4);
        const auto scenarioParamPrefix = "Root.SpatialParams." + getParam<std::string>("Root.SpatialParams.Scenario", "C12a") + ".";

        // read the tabularized root conductivities
        axialRootConductivity_.resize(2);
        radialRootConductivity_.resize(2);
        axialRootConductivity_[0] = { getParam<std::vector<Scalar>>(scenarioParamPrefix + "AxialConductivites.Order0.Age"),
                                      getParam<std::vector<Scalar>>(scenarioParamPrefix + "AxialConductivites.Order0.Kx") };
        axialRootConductivity_[1] = { getParam<std::vector<Scalar>>(scenarioParamPrefix + "AxialConductivites.Order1.Age"),
                                      getParam<std::vector<Scalar>>(scenarioParamPrefix + "AxialConductivites.Order1.Kx") };
        radialRootConductivity_[0] = { getParam<std::vector<Scalar>>(scenarioParamPrefix + "RadialConductivites.Order0.Age"),
                                       getParam<std::vector<Scalar>>(scenarioParamPrefix + "RadialConductivites.Order0.Kr") };
        radialRootConductivity_[1] = { getParam<std::vector<Scalar>>(scenarioParamPrefix + "RadialConductivites.Order1.Age"),
                                       getParam<std::vector<Scalar>>(scenarioParamPrefix + "RadialConductivites.Order1.Kr") };

        // sanity checks
        for (const auto& k : axialRootConductivity_)
            if (k.first.size() != k.second.size())
                DUNE_THROW(Dune::IOError, "AxialConductivites.Age and AxialConductivites.Kx have to have the same length!");

        for (const auto& k : radialRootConductivity_)
            if (k.first.size() != k.second.size())
                DUNE_THROW(Dune::IOError, "RadialConductivites.Age and RadialConductivites.Kr have to have the same length!");

        const auto& gv = gridGeometry->gridView();
        radii_.resize(gv.size(0));
        age_.resize(gv.size(0));
        order_.resize(gv.size(0));
        Kx_.resize(gv.size(0));
        Kr_.resize(gv.size(0));
        const auto paramType = getParam<std::string>("Root.Grid.DGFParam", "Radius");
        for (const auto& element : elements(gv))
        {
            const auto eIdx = gridGeometry->elementMapper().index(element);
            auto level0element = element;
            while (level0element.hasFather())
                level0element = level0element.father();

            if (paramType == "Surface")
            {
                // see e.g. lupine.dgf
                const auto rootLength = element.geometry().volume();
                const auto rootSurface = gridData->parameters(level0element)[2]/(1 << element.level());
                radii_[eIdx] = rootSurface / ( rootLength * 2.0 * M_PI );
            }
            else if (paramType == "Radius")
            {
                static const auto radiusIdx = getParam<int>("Root.Grid.DGFRadiusIndex", 4);
                static const auto radiusScaling = getParam<double>("Root.Grid.DGFRadiusScaling", 1.0);
                radii_[eIdx] = gridData->parameters(level0element)[radiusIdx]*radiusScaling; // read radius from grid file
            }
            else
                DUNE_THROW(Dune::NotImplemented, "Unknown DGF Parameter " << paramType);

            // read age
            static const auto ageIdx = getParam<int>("Root.Grid.DGFAgeIndex", 7);
            static const auto ageScaling = getParam<int>("Root.Grid.DGFAgeScaling", 1.0);
            static const auto currentPlantAge = getParam<int>("Root.Grid.DGFPlantAge", 8.0);
            age_[eIdx] = std::max(0.0, currentPlantAge - gridData->parameters(level0element)[ageIdx]*ageScaling);

            // read order
            static const auto orderIdx = getParam<int>("Root.Grid.DGFOrderIndex", 0);
            order_[eIdx] = gridData->parameters(level0element)[orderIdx];

            // permeability
            const auto order = static_cast<std::size_t>(std::min(std::max(0.0, order_[eIdx]), 1.0));
            Kx_[eIdx] = computeAxialConductivity(order, age_[eIdx]);
            Kr_[eIdx] = computeRadialConductivity(order, age_[eIdx]);
        }
    }

    /*!
     * \brief Returns how much the domain is extruded at a given sub-control volume.
     *
     * The extrusion factor here extrudes the 1d line to a circular tube with
     * cross-section area pi*r^2.
     */
    template<class ElementSolution>
    Scalar extrusionFactor(const Element &element,
                           const SubControlVolume &scv,
                           const ElementSolution& elemSol) const
    {
        const auto eIdx = this->gridGeometry().elementMapper().index(element);
        const auto radius = this->radius(eIdx);
        return M_PI*radius*radius;
    }

    /*!
     * \brief Returns the intrinsic permeability for the current sub-control volume in [m^2].
     *
     * \note Kx has units [m^4/(Pa*s)] so we have to divide by the cross-section area
     *       and multiply with a characteristic viscosity
     */
    template<class ElementSolution>
    PermeabilityType permeability(const Element& element,
                                  const SubControlVolume& scv,
                                  const ElementSolution& elemSol) const
    {
        const auto eIdx = this->gridGeometry().elementMapper().index(element);
        const Scalar r = radius(eIdx);
        return Kx(eIdx) / (M_PI*r*r) * Components::SimpleH2O<Scalar>::liquidViscosity(285.15, 1e5);
    }

    /*!
     * \brief Returns the radius of the circular pipe for the current sub-control volume in [m].
     *
     * \param eIdxGlobal the index of the element
     */
    Scalar radius(std::size_t eIdxGlobal) const
    {
        return radii_[eIdxGlobal];
    }

    /*!
     * \brief Returns the radial permeability.
     *
     * \param eIdxGlobal the index of the element
     */
    Scalar Kr(std::size_t eIdxGlobal) const
    {
        return Kr_[eIdxGlobal];
    }

    /*!
     * \brief Returns the radial permeability.
     *
     * \param eIdxGlobal the index of the element
     */
    Scalar Kx(std::size_t eIdxGlobal) const
    {
        return Kx_[eIdxGlobal];
    }

    // get radii for output
    const std::vector<Scalar>& getRadii() const { return radii_; }
    // get ages for output
    const std::vector<Scalar>& getAges() const { return age_; }
    // get orders for output
    const std::vector<Scalar>& getOrders() const { return order_; }
    // get kx for output
    const std::vector<Scalar>& getKx() const { return Kx_; }
    // get kr for output
    const std::vector<Scalar>& getKr() const { return Kr_; }

    /*!
     * \brief Returns the porosity \f$[-]\f$.
     *
     * \param globalPos the scv center
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    {
        return porosity_;
    }

    //! compute the radial conductivity (m/Pa/s) given the segment age in days
    double computeRadialConductivity(int order, double age) const
    {
        return interpolate<InterpolationPolicy::LinearTable>(age, radialRootConductivity_[order]);
    }

    //! compute the axial conductivity in (m^4/Pa/s) given the segment age in days
    double computeAxialConductivity(int order, double age) const
    {
        return interpolate<InterpolationPolicy::LinearTable>(age, axialRootConductivity_[order]);
    }

    /*!
     * \brief Returns the temperature within the domain in [K].
     */
    Scalar temperatureAtPos(const GlobalPosition& globalPos) const
    { return 273.15 + 10.0; }

private:
    std::shared_ptr<const GridData<Grid>> gridData_;
    Scalar porosity_;

    //! Tabularized root conductivity functions (pairs of x, y) for each order
    std::vector<std::pair<std::vector<double>, std::vector<double>>> axialRootConductivity_;
    std::vector<std::pair<std::vector<double>, std::vector<double>>> radialRootConductivity_;

    std::vector<Scalar> radii_;
    std::vector<Scalar> age_;
    std::vector<Scalar> order_;
    std::vector<Scalar> Kx_;
    std::vector<Scalar> Kr_;
};

} // end namespace Dumux

#endif
