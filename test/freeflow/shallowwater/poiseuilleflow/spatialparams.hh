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
 * \ingroup ShallowWaterTests
 * \brief The spatial parameters for the Poiseuille flow problem.
 */
#ifndef DUMUX_POISEUILLE_FLOW_SPATIAL_PARAMETERS_HH
#define DUMUX_POISEUILLE_FLOW_SPATIAL_PARAMETERS_HH

#include <dumux/material/spatialparams/fv.hh>
#include <dumux/common/parameters.hh>
#include <dumux/material/fluidmatrixinteractions/frictionlaws/frictionlaw.hh>
#include <dumux/material/fluidmatrixinteractions/frictionlaws/manning.hh>
#include <dumux/material/fluidmatrixinteractions/frictionlaws/nikuradse.hh>
// #include <dumux/material/fluidmatrixinteractions/wallfrictionlaws/wallfrictionlaw.hh>
// #include <dumux/material/fluidmatrixinteractions/wallfrictionlaws/wallnoslip.hh>


namespace Dumux {

/*!
 * \ingroup ShallowWaterTests
 * \brief The spatial parameters class for the Poiseuille flow test.
 *
 */
template<class GridGeometry, class Scalar, class VolumeVariables>
class PoiseuilleFlowSpatialParams
: public FVSpatialParams<GridGeometry, Scalar,
                         PoiseuilleFlowSpatialParams<GridGeometry, Scalar, VolumeVariables>>
{
    using ThisType = PoiseuilleFlowSpatialParams<GridGeometry, Scalar, VolumeVariables>;
    using ParentType = FVSpatialParams<GridGeometry, Scalar, ThisType>;
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    PoiseuilleFlowSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        gravity_ = getParam<Scalar>("Problem.Gravity");
        bedSlope_ = getParam<Scalar>("Problem.BedSlope");
        wallFrictionLawType_ = getParam<std::string>("Problem.WallFrictionLaw");
        //initWallFrictionLaw();
    }

    /*!
     * \brief Initialize the WallFrictionLaw
     */
    //void initWallFrictionLaw()
    //{
    //  if (wallFrictionLawType_ == "Noslip")
    //  {
    //      Scalar alphaWall = getParam<Scalar>("Problem.alphaWall");
    //      wallFrictionLaw_ = std::make_unique<WallFrictionLawNoslip<VolumeVariables>>(gravity_, alphaWall);
    //  }
    //  else if (wallFrictionLawType_ == "Nikuradse")
    //  {
    //      Scalar wallKs = getParam<Scalar>("Problem.wallKs"); // equivalent sand roughness of the wall
    //      wallFrictionLaw_ = std::make_unique<WallFrictionLawNikuradse<VolumeVariables>>(wallKs);
    //  }
    //  else
    //  {
    //      std::cout<<"The WallFrictionLaw in params.input is unknown. Valid entries are 'Noslip' and 'Nikuradse'!"<<std::endl;
    //  }
    //}

    /*! \brief Define the gravitation.
    *
    * \return gravity constant
    */
    Scalar gravity(const GlobalPosition& globalPos) const
    {
        return gravity_;
    }

    /*! \brief Define the gravitation.
    *
    * \return gravity constant
    */
    Scalar gravity() const
    {
        return gravity_;
    }

    /*! \brief Get the wallFrictionLaw.
    *
    * Get the wallFrictionLaw, which already includes the wall friction value.
    *
    * \return wallFrictionLaw
    */

    //const WallFrictionLaw<VolumeVariables>& wallFrictionLaw(const Element& element,
    //                                                const SubControlVolumeFace& scvf) const
    //{
    //    return *wallFrictionLaw_;
    //}

    /*! \brief Define the bed surface
    *
    * \param element The current element
    * \param scv The sub-control volume inside the element.
    *
    * \return The bed surface
    */
    Scalar bedSurface(const Element& element,
                      const SubControlVolume& scv) const
    {
        // todo depends on index e.g. eIdx = scv.elementIndex();
        return 9.98 - element.geometry().center()[0] * bedSlope_;
    }

private:
    Scalar gravity_;
    Scalar bedSlope_;
    std::string wallFrictionLawType_;
    //std::unique_ptr<FrictionLaw<VolumeVariables>> wallFrictionLaw_;
};

} // end namespace Dumux

#endif
