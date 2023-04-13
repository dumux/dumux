// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
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

namespace Dumux::PoreNetwork {

/*!
 * \ingroup SpatialParameters
 * \ingroup PoreNetworkModels
 * \brief The base class for spatial parameters for pore-network models.
 */
template<class GridGeometry, class Scalar, class Implementation>
class SpatialParams
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

    SpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
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

};

} // namespace Dumux::PoreNetwork

#endif
