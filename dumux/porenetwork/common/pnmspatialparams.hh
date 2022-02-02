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
 * \ingroup PoreNetworkModels
 * \ingroup SpatialParameters
 * \brief The base class for spatial parameters for pore-network models.
 */
#ifndef DUMUX_PNM_SPATIAL_PARAMS_HH
#define DUMUX_PNM_SPATIAL_PARAMS_HH

#include <type_traits>
#include <memory>

#include <dune/common/fvector.hh>
#include <dumux/common/parameters.hh>
#include <dumux/porousmediumflow/fvspatialparams.hh>
#include <dumux/porenetwork/common/spatialparamstraits_.hh>

namespace Dumux::PoreNetwork {

/*!
 * \ingroup SpatialParameters
 * \ingroup PoreNetworkModels
 * \brief The base class for spatial parameters for pore-network models.
 */
template<class GridGeometry, class Scalar, class Implementation>
class PNMSpatialParams
: public FVPorousMediumFlowSpatialParams<GridGeometry, Scalar, Implementation>
{
    using ParentType = FVPorousMediumFlowSpatialParams<GridGeometry, Scalar, Implementation>;
    using GridView = typename GridGeometry::GridView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    static constexpr auto dimWorld = GridView::dimensionworld;

public:
    using PermeabilityType = Scalar;

    PNMSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {}

    /*!
     * \brief Length of the throat \f$[m]\f$.
     *        Can be solution-dependent.
     *
     *  \param element The finite volume element
     *  \param elemVolVars The element volume variables.
     */
    template<class ElementVolumeVariables>
    Scalar throatLength(const Element& element,
                        const ElementVolumeVariables& elemVolVars) const
    {
        const auto eIdx = this->gridGeometry().elementMapper().index(element);
        return this->gridGeometry().throatLength(eIdx);
    }

    /*!
     * \brief Inscribed radius of the throat \f$[m]\f$.
     *        Can be solution-dependent.
     *
     *  \param element The finite volume element
     *  \param elemVolVars The element volume variables.
     */
    template<class ElementVolumeVariables>
    Scalar throatInscribedRadius(const Element& element,
                                 const ElementVolumeVariables& elemVolVars) const
    {
        const auto eIdx = this->gridGeometry().elementMapper().index(element);
        return this->gridGeometry().throatInscribedRadius(eIdx);
    }

    /*!
     * \brief Cross-sectional area of the throat \f$[m]\f$.
     *        Can be solution-dependent.
     *
     *  \param element The finite volume element
     *  \param elemVolVars The element volume variables.
     */
    template<class ElementVolumeVariables>
    Scalar throatCrossSectionalArea(const Element& element,
                                   const ElementVolumeVariables& elemVolVars) const
    {
        const auto eIdx = this->gridGeometry().elementMapper().index(element);
        return this->gridGeometry().throatCrossSectionalArea(eIdx);
    }

   /*!
    * \brief Inscribed radius of the pore body \f$[m]\f$.
    *        Can be solution-dependent.
    *
    *  \param element The finite volume element
    *  \param scv The sub-control volume
    *  \param elemSol The element solution
    */
    template<class ElementSolutionVector>
    Scalar poreInscribedRadius(const Element& element,
                               const SubControlVolume& scv,
                               const ElementSolutionVector& elemSol) const
    { return this->gridGeometry().poreInscribedRadius(scv.dofIndex()); }

    /*!
     * \brief Returns a reference to the gridview
     */
    const GridView& gridView() const
    { return this->gridGeometry().gridView(); }

    //! Required for compatibility reasons with porous medium-flow models.
    Scalar permeabilityAtPos(const GlobalPosition& globalPos) const
    { return 1.0; }

    //! Required for compatibility reasons with porous medium-flow models.
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return 1.0; }

private:
};

} // namespace Dumux::PoreNetwork

#endif
