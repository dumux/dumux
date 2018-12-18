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
 *   along with this program. If not, see <http://www.gnu.org/licenses/>.    *
 *****************************************************************************/
/*!
 * \file
 * \ingroup RichardsTests
 * \brief Spatial parameters for the RichardsAnalyticalProblem.
 */

#ifndef DUMUX_RICHARDS_ANALYTICAL_SPATIAL_PARAMETERS_HH
#define DUMUX_RICHARDS_ANALYTICAL_SPATIAL_PARAMETERS_HH

#include <dumux/material/spatialparams/fv.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedvangenuchten.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh>
#include <dumux/material/fluidmatrixinteractions/2p/brookscorey.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>

#include <dumux/porousmediumflow/richards/model.hh>

namespace Dumux {

/*!
 * \ingroup RichardsModel
 * \ingroup ImplicitTestProblems
 * \brief The spatial parameters for the RichardsAnalyticalProblem.
 */
template<class FVGridGeometry, class Scalar>
class RichardsAnalyticalSpatialParams
: public FVSpatialParams<FVGridGeometry, Scalar,
                         RichardsAnalyticalSpatialParams<FVGridGeometry, Scalar>>
{
    using GridView = typename FVGridGeometry::GridView;
    using ParentType = FVSpatialParams<FVGridGeometry, Scalar,
                                       RichardsAnalyticalSpatialParams<FVGridGeometry, Scalar>>;

    enum { dimWorld=GridView::dimensionworld };
    using Element = typename GridView::template Codim<0>::Entity;

    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using EffectiveLaw = LinearMaterial<Scalar>;

public:

    using MaterialLaw = EffToAbsLaw<EffectiveLaw>;
    using MaterialLawParams = typename MaterialLaw::Params;
    // export permeability type
    using PermeabilityType = Scalar;

    RichardsAnalyticalSpatialParams(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
        : ParentType(fvGridGeometry)
    {
        K_ = 5e-12;
        materialParams_.setSwr(0.0);
        materialParams_.setSnr(0.0);
        materialParams_.setEntryPc(0);
        materialParams_.setMaxPc(1e10);
    }

    /*!
     * \brief Returns the intrinsic permeability tensor [m^2] at a given location
     *
     * \param globalPos The global position where we evaluate
     */
    PermeabilityType permeabilityAtPos(const GlobalPosition& globalPos) const
    {
        return K_;
    }

    /*!
     * \brief Returns the porosity [] at a given location
     *
     * \param globalPos The global position where we evaluate
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return 0.4; }

    /*!
     * \brief Returns the parameters for the material law at a given location
     *
     * This method is not actually required by the Richards model, but provided
     * for the convenience of the RichardsLensProblem
     *
     * \param globalPos A global coordinate vector
     */
    const MaterialLawParams& materialLawParamsAtPos(const GlobalPosition &globalPos) const
    {
        return materialParams_;
    }

private:
    // intrinsic permeability
    Scalar K_;

    MaterialLawParams materialParams_;
};
} // end namespace

#endif
