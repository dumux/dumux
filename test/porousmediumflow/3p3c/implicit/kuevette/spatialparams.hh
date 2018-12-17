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
 * \brief Definition of the spatial parameters for the kuevette problem.
 */
#ifndef DUMUX_KUEVETTE3P3CNI_SPATIAL_PARAMS_HH
#define DUMUX_KUEVETTE3P3CNI_SPATIAL_PARAMS_HH

#include <dune/common/float_cmp.hh>

#include <dumux/porousmediumflow/properties.hh>
#include <dumux/material/spatialparams/fv.hh>
#include <dumux/material/fluidmatrixinteractions/3p/regularizedparkervangen3p.hh>
#include <dumux/material/fluidmatrixinteractions/3p/regularizedparkervangen3pparams.hh>
#include <dumux/material/fluidmatrixinteractions/3p/efftoabslaw.hh>

namespace Dumux {

/*!
 * \ingroup ThreePThreeCTests
 * \brief Definition of the spatial parameters for the kuevette problem
 */
template<class FVGridGeometry, class Scalar>
class KuevetteSpatialParams
: public FVSpatialParams<FVGridGeometry, Scalar,
                         KuevetteSpatialParams<FVGridGeometry, Scalar>>
{
    using GridView = typename FVGridGeometry::GridView;
    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;

    using Element = typename GridView::template Codim<0>::Entity;
    using ParentType = FVSpatialParams<FVGridGeometry, Scalar,
                                       KuevetteSpatialParams<FVGridGeometry, Scalar>>;

    using GlobalPosition = typename SubControlVolume::GlobalPosition;

    using EffectiveLaw = RegularizedParkerVanGen3P<Scalar>;

public:
    using MaterialLaw = EffToAbsLaw<EffectiveLaw>;
    using MaterialLawParams = typename MaterialLaw::Params;
    using PermeabilityType = Scalar;

    /*!
     * \brief The constructor
     *
     * \param fvGridGeometry The finite volume grid geometry
     */
    KuevetteSpatialParams(std::shared_ptr<const FVGridGeometry> fvGridGeometry)
    : ParentType(fvGridGeometry)
    {
        // intrinsic permeabilities
        fineK_ = 6.28e-12;
        coarseK_ = 9.14e-10;

        // porosities
        finePorosity_ = 0.42;
        coarsePorosity_ = 0.42;

        // residual saturations
        fineMaterialParams_.setSwr(0.12);
        fineMaterialParams_.setSnr(0.07);
        fineMaterialParams_.setSgr(0.01);
        coarseMaterialParams_.setSwr(0.12);
        coarseMaterialParams_.setSnr(0.07);
        coarseMaterialParams_.setSgr(0.01);

        // parameters for the 3phase van Genuchten law
        fineMaterialParams_.setVgAlpha(0.0005);
        coarseMaterialParams_.setVgAlpha(0.005);
        fineMaterialParams_.setVgn(4.0);
        coarseMaterialParams_.setVgn(4.0);

        coarseMaterialParams_.setKrRegardsSnr(true);
        fineMaterialParams_.setKrRegardsSnr(true);

        // parameters for adsorption
        coarseMaterialParams_.setKdNAPL(0.);
        coarseMaterialParams_.setRhoBulk(1500.);
        fineMaterialParams_.setKdNAPL(0.);
        fineMaterialParams_.setRhoBulk(1500.);
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
     * \brief Define the porosity \f$[-]\f$ of the spatial parameters
     *
     * \param globalPos The position for which the bc type should be evaluated
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    {
        if (isFineMaterial_(globalPos))
            return finePorosity_;
        else
            return coarsePorosity_;
    }


    /*!
     * \brief Function for defining the parameters needed by constitutive relationships (kr-sw, pc-sw, etc.).
     *
     * \param element The current element
     * \param scv The sub-control volume inside the element.
     * \param elemSol The solution at the dofs connected to the element.
     * \return the material parameters object
     */
    template<class ElementSolution>
    const MaterialLawParams& materialLawParams(const Element& element,
                                               const SubControlVolume& scv,
                                               const ElementSolution& elemSol) const
    {
        const auto& globalPos = scv.dofPosition();
        if (isFineMaterial_(globalPos))
            return fineMaterialParams_;
        else
            return coarseMaterialParams_;
    }

private:
    bool isFineMaterial_(const GlobalPosition &globalPos) const
    {
        return ((Dune::FloatCmp::ge<Scalar>(globalPos[0], 0.13)
                    && Dune::FloatCmp::le<Scalar>(globalPos[0], 1.24)
                    && Dune::FloatCmp::ge<Scalar>(globalPos[1], 0.32)
                    && Dune::FloatCmp::le<Scalar>(globalPos[1], 0.60))
                || (Dune::FloatCmp::ge<Scalar>(globalPos[0], 1.20)
                    && Dune::FloatCmp::le<Scalar>(globalPos[1], 0.15)));
    }

    Scalar fineK_;
    Scalar coarseK_;

    Scalar finePorosity_;
    Scalar coarsePorosity_;

    MaterialLawParams fineMaterialParams_;
    MaterialLawParams coarseMaterialParams_;

};

} // end namespace Dumux

#endif
