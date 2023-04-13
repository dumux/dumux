// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup EmbeddedTests
 * \brief The spatial parameters class blood flow problem.
 */

#ifndef DUMUX_ROOT_SPATIALPARAMS_HH
#define DUMUX_ROOT_SPATIALPARAMS_HH

#include <dumux/common/parameters.hh>
#include <dumux/io/grid/griddata.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/porousmediumflow/fvspatialparams1p.hh>

namespace Dumux {

/*!
 * \ingroup EmbeddedTests
 * \brief Definition of the spatial parameters for the blood flow problem.
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

    //! Indices to access the parameters in the dgf file
    enum DGFParamIndicesSurfaceModel {
        orderIdx = 0,
        rootIdIdx = 1,
        surfaceIdx = 2,
        massIdx = 3,
        plantIdx = 5
    };

public:
    // export permeability type
    using PermeabilityType = Scalar;

    RootSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry,
                      std::shared_ptr<const GridData<Grid>> gridData)
    : ParentType(gridGeometry), gridData_(gridData)
    {
        porosity_ = getParam<Scalar>("Root.SpatialParams.Porosity", 0.4);
        constantKx_ = getParam<Scalar>("Root.SpatialParams.Kx", 5.0968e-17);
        constantKr_ = getParam<Scalar>("Root.SpatialParams.Kr", 2.04e-13);

        const auto& gv = gridGeometry->gridView();
        radii_.resize(gv.size(0));
        for (const auto& element : elements(gv))
        {
            const auto eIdx = gv.indexSet().index(element);
            auto level0element = element;
            for (auto levelIdx = element.level(); levelIdx != 0; levelIdx--)
                level0element = level0element.father();

            const auto& params = gridData_->parameters(level0element);
            if (params.size() == 1) // assume only radius is given
                radii_[eIdx] = params[0];
            else // assume DGFParamIndicesSurfaceModel
            {
                const Scalar rootLength = element.geometry().volume();
                const Scalar rootSurface = gridData_->parameters(level0element)[DGFParamIndicesSurfaceModel::surfaceIdx]/(1 << element.level());
                radii_[eIdx] = rootSurface / rootLength / 2.0 / M_PI;
            }
        }
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
        const Scalar r = radius(this->gridGeometry().elementMapper().index(element));
        return constantKx_ / (M_PI*r*r) * Components::SimpleH2O<Scalar>::liquidViscosity(285.15, 1e5);
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
        return constantKr_;
    }

    const std::vector<Scalar>& getRadii() const
    { return radii_; }

    /*!
     * \brief Returns the porosity \f$[-]\f$.
     *
     * \param globalPos the scv center
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    {
        return porosity_;
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
        const auto r = radius(eIdx);
        return M_PI*r*r;
    }

private:
    std::shared_ptr<const GridData<Grid>> gridData_;
    Scalar porosity_, constantKx_, constantKr_;
    std::vector<Scalar> radii_;
};

} // end namespace Dumux

#endif
