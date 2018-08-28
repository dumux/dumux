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
 * \ingroup RichardsTests
 * \brief Definition of the spatial parameters for the non-isothermal Richards problems.
 */
#ifndef DUMUX_RICHARDSNI_SPATIAL_PARAMS_HH
#define DUMUX_RICHARDSNI_SPATIAL_PARAMS_HH

#include <dumux/porousmediumflow/richards/model.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedvangenuchten.hh>
#include <dumux/material/spatialparams/fv.hh>

namespace Dumux {

template<class FVGridGeometry, class Scalar>
class RichardsNISpatialParams
: public FVSpatialParams<FVGridGeometry, Scalar,
                         RichardsNISpatialParams<FVGridGeometry, Scalar>>
{
    using GridView = typename FVGridGeometry::GridView;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;
    using ParentType = FVSpatialParams<FVGridGeometry, Scalar,
                                       RichardsNISpatialParams<FVGridGeometry, Scalar>>;

    enum { dimWorld=GridView::dimensionworld };

    using EffectiveLaw = RegularizedVanGenuchten<Scalar>;

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    // export permeability type
    using PermeabilityType = Scalar;

    using MaterialLaw = EffToAbsLaw<EffectiveLaw>;
    using  MaterialLawParams = typename MaterialLaw::Params;

    RichardsNISpatialParams(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
        : ParentType(fvGridGeometry)
    {
        permeability_ = 1e-10;
        porosity_ = 0.4;


        // residual saturations
        materialParams_.setSwr(0.05);
        materialParams_.setSnr(0.0);

        // parameters for the Van Genuchten law
        // alpha and n

        materialParams_.setVgAlpha(0.0037);
        materialParams_.setVgn(4.7);
    }

    /*!
     * \brief Define the intrinsic permeability \f$\mathrm{[m^2]}\f$.
     *
     * \param element The current finite element
     * \param fvGeometry The current finite volume geometry of the element
     * \param scvIdx The index of the sub-control volume
     */
    PermeabilityType permeabilityAtPos(const GlobalPosition& globalPos) const
    {
        return permeability_;
    }

    /*!
     * \brief Define the porosity \f$\mathrm{[-]}\f$.
     *
     * \param element The finite element
     * \param fvGeometry The finite volume geometry
     * \param scvIdx The local index of the sub-control volume where
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    {
        return porosity_;
    }

        /*!
     * \brief return the parameter object for the Brooks-Corey material law which depends on the position
     *
     * \param element The current finite element
     * \param fvGeometry The current finite volume geometry of the element
     * \param scvIdx The index of the sub-control volume
     */
     const MaterialLawParams& materialLawParamsAtPos(const GlobalPosition &globalPos) const
    {
        return materialParams_;
    }

private:

    MaterialLawParams materialParams_;
    Scalar permeability_;
    Scalar porosity_;
};

} // end namespace Dumux

#endif
