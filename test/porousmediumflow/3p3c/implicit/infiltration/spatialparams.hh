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
 * \ingroup ThreePThreeCTests
 * \brief Definition of the spatial parameters for the kuevette problem
 *        which uses the isothermal two-phase two component
 *        fully implicit model.
 */
#ifndef DUMUX_INFILTRATION_THREEPTHREEC_SPATIAL_PARAMETERS_HH
#define DUMUX_INFILTRATION_THREEPTHREEC_SPATIAL_PARAMETERS_HH

#include <dumux/porousmediumflow/properties.hh>
#include <dumux/material/spatialparams/fv.hh>
#include <dumux/material/fluidmatrixinteractions/3p/regularizedparkervangen3p.hh>
#include <dumux/material/fluidmatrixinteractions/3p/regularizedparkervangen3pparams.hh>
#include <dumux/material/fluidmatrixinteractions/3p/efftoabslaw.hh>
#include <dumux/io/plotmateriallaw3p.hh>

namespace Dumux {

/*!
 * \ingroup ThreePThreeCTests
 * \brief Definition of the spatial parameters for the infiltration problem
 */
template<class FVGridGeometry, class Scalar>
class InfiltrationThreePThreeCSpatialParams
: public FVSpatialParams<FVGridGeometry, Scalar,
                         InfiltrationThreePThreeCSpatialParams<FVGridGeometry, Scalar>>
{
    using GridView = typename FVGridGeometry::GridView;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;
    using ParentType = FVSpatialParams<FVGridGeometry, Scalar,
                                       InfiltrationThreePThreeCSpatialParams<FVGridGeometry, Scalar>>;

    using EffectiveLaw = RegularizedParkerVanGen3P<Scalar>;

    using GlobalPosition = typename SubControlVolume::GlobalPosition;

public:
   using MaterialLaw = EffToAbsLaw<EffectiveLaw>;
   using MaterialLawParams = typename MaterialLaw::Params;
   using PermeabilityType = Scalar;

    /*!
     * \brief The constructor
     *
     * \param fvGridGeometry The finite volume grid geometry
     */
    InfiltrationThreePThreeCSpatialParams(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
    {
        // intrinsic permeabilities
        fineK_ = 1.e-11;
        coarseK_ = 1.e-11;

        // porosities
        porosity_ = 0.40;

        // residual saturations
        materialParams_.setSwr(0.12);
        materialParams_.setSnr(0.07);
        materialParams_.setSgr(0.03);

        // parameters for the 3phase van Genuchten law
        materialParams_.setVgAlpha(0.0005);
        materialParams_.setVgn(4.);
        materialParams_.setKrRegardsSnr(false);

        // parameters for adsorption
        materialParams_.setKdNAPL(0.);
        materialParams_.setRhoBulk(1500.);
    }

    /*!
     * \brief Function for defining the (intrinsic) permeability \f$[m^2]\f$
     * \note  It is possibly solution dependent.
     *
     * \param element The current element
     * \param scv The sub-control volume inside the element.
     * \param elemSol The solution at the dofs connected to the element.
     * \return permeability
     */
    template<class ElementSolution>
    PermeabilityType permeability(const Element& element,
                                  const SubControlVolume& scv,
                                  const ElementSolution& elemSol) const
    {
        const auto& globalPos = scv.dofPosition();
        if (isFineMaterial_(globalPos))
            return fineK_;
        return coarseK_;
    }

    /*!
     * \brief Returns the porosity \f$[-]\f$
     *
     * \param globalPos The position for which the porosity is evaluated
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    {
        return porosity_;
    }


    /*!
     * \brief return the parameter object for the material law which depends on the position
     *
     * \param globalPos The position for which the material law should be evaluated
     */
    const MaterialLawParams& materialLawParamsAtPos(const GlobalPosition& globalPos) const
    {
        return materialParams_;
    }

private:
    bool isFineMaterial_(const GlobalPosition &globalPos) const
    { return
            70.0 <= globalPos[0] && globalPos[0] <= 85.0 &&
            7.0 <= globalPos[1] && globalPos[1] <= 7.50;
    }

    Scalar fineK_;
    Scalar coarseK_;

    Scalar porosity_;

    MaterialLawParams materialParams_;
};

} // end namespace Dumux

#endif
