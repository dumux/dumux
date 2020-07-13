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
 *
 * \ingroup SpatialParameters
 * \brief The base class for spatial parameters for pore network models.
 */
#ifndef DUMUX_PNM_SPATIAL_PARAMS_1P_HH
#define DUMUX_PNM_SPATIAL_PARAMS_1P_HH

#include "porenetworkbase.hh"

namespace Dumux
{

/*!
 * \ingroup SpatialParameters
 */

  template<class GridGeometry, class Scalar, class SinglePhaseTransmissibilityLaw, class Implementation>
  class PNMOnePBaseSpatialParams;

  /**
   * \brief The base class for spatial parameters for pore network models.
   */
  template<class GridGeometry, class Scalar, class SinglePhaseTransmissibilityLaw>
  class PNMOnePSpatialParams : public PNMOnePBaseSpatialParams<GridGeometry, Scalar, SinglePhaseTransmissibilityLaw,
                                                               PNMOnePSpatialParams<GridGeometry, Scalar, SinglePhaseTransmissibilityLaw>>
  {
      using ParentType = PNMOnePBaseSpatialParams<GridGeometry, Scalar, SinglePhaseTransmissibilityLaw,
                                                  PNMOnePSpatialParams<GridGeometry, Scalar, SinglePhaseTransmissibilityLaw>>;
  public:
      using ParentType::ParentType;
  };

/**
 * \brief The base class for spatial parameters for pore network models.
 */
template<class GridGeometry, class Scalar, class SinglePhaseTransmissibilityLaw, class Implementation>
class PNMOnePBaseSpatialParams : public PNMBaseSpatialParams<GridGeometry, Scalar, Implementation>
{
    using ParentType = PNMBaseSpatialParams<GridGeometry, Scalar, Implementation>;

    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;

public:
    using PermeabilityType = Scalar;
    using ParentType::ParentType;

    /*!
    * \brief Returns the transmissibility of a throat
    */
   template<class ElementVolumeVariables, class FluxVariablesCache>
   Scalar transmissibility(const Element& element,
                           const FVElementGeometry& fvGeometry,
                           const SubControlVolumeFace& scvf,
                           const ElementVolumeVariables& elemVolVars,
                           const FluxVariablesCache& fluxVarsCache,
                           const int phaseIdx = 0) const
   {
       // forward to specialized function
       return SinglePhaseTransmissibilityLaw::singlePhaseTransmissibility(element, fvGeometry, scvf, fluxVarsCache);
   }

    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return 1.0; }

};

} // namespace Dumux

#endif
