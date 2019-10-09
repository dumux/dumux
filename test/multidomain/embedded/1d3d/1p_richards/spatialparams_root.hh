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
 * \ingroup EmbeddedTests
 * \brief The spatial parameters class blood flow problem.
 */

#ifndef DUMUX_ROOT_SPATIALPARAMS_HH
#define DUMUX_ROOT_SPATIALPARAMS_HH

#include <dumux/common/parameters.hh>
#include <dumux/io/grid/griddata.hh>
#include <dumux/material/components/simpleh2o.hh>
#include <dumux/material/spatialparams/fv1p.hh>

namespace Dumux {

/*!
 * \ingroup EmbeddedTests
 * \brief Definition of the spatial parameters for the blood flow problem.
 */
template<class GridGeometry, class Scalar>
class RootSpatialParams
: public FVSpatialParamsOneP<GridGeometry, Scalar, RootSpatialParams<GridGeometry, Scalar>>
{
    using ThisType = RootSpatialParams<GridGeometry, Scalar>;
    using ParentType = FVSpatialParamsOneP<GridGeometry, Scalar, ThisType>;
    using Grid = typename GridGeometry::Grid;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    //! Indices to access the parameters in the dgf file
    enum DGFParamIndices {
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
        constantKx_ = getParam<Scalar>("SpatialParams.Kx", 5.0968e-17);
        constantKr_ = getParam<Scalar>("SpatialParams.Kr", 2.04e-13);

        const auto& gv = gridGeometry->gridView();
        radii_.resize(gv.size(0));
        for (const auto& element : elements(gv))
        {
            const auto eIdx = gv.indexSet().index(element);
            auto level0element = element;
            for(auto levelIdx = element.level(); levelIdx != 0; levelIdx--)
                level0element = level0element.father();
            const Scalar rootLength = element.geometry().volume();
            const Scalar rootSurface = gridData_->parameters(level0element)[DGFParamIndices::surfaceIdx]/(1 << element.level());
            radii_[eIdx] = rootSurface / rootLength / 2.0 / M_PI;
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

private:
    std::shared_ptr<const GridData<Grid>> gridData_;
    Scalar porosity_, constantKx_, constantKr_;
    std::vector<Scalar> radii_;
};

} // end namespace Dumux

#endif
