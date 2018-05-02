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
 * \ingroup TwoPTests
 * \brief The spatial params for the incompressible 2p test
 */
#ifndef DUMUX_INCOMPRESSIBLE_TWOPBOXDFM_TEST_SPATIAL_PARAMS_HH
#define DUMUX_INCOMPRESSIBLE_TWOPBOXDFM_TEST_SPATIAL_PARAMS_HH

#include <dumux/discretization/methods.hh>

#include <dumux/material/spatialparams/fv.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedvangenuchten.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>

#include <dumux/porousmediumflow/2p/boxmaterialinterfaceparams.hh>

namespace Dumux {

/*!
 * \ingroup TwoPTests
 * \brief The spatial params for the incompressible 2p test
 */
template<class FVGridGeometry, class Scalar>
class TwoPTestSpatialParams : public FVSpatialParams< FVGridGeometry, Scalar, TwoPTestSpatialParams<FVGridGeometry, Scalar> >
{
    using ThisType = TwoPTestSpatialParams<FVGridGeometry, Scalar>;
    using ParentType = FVSpatialParams<FVGridGeometry, Scalar, ThisType>;

    using GridView = typename FVGridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using FVElementGeometry = typename FVGridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;

    static constexpr int dimWorld = GridView::dimensionworld;

public:
    using MaterialLaw = EffToAbsLaw< RegularizedBrooksCorey<Scalar> >;
    using MaterialLawParams = typename MaterialLaw::Params;
    using PermeabilityType = Scalar;

    TwoPTestSpatialParams(std::shared_ptr<const FVGridGeometry> fvGridGeometry) : ParentType(fvGridGeometry)
    {
        // residual saturations
        matrixMaterialParams_.setSwr(0.18);
        matrixMaterialParams_.setSnr(0.0);
        lensMaterialParams_.setSwr(0.28);
        lensMaterialParams_.setSnr(0.0);
        fractureMaterialParams_.setSwr(0.05);
        fractureMaterialParams_.setSnr(0.0);

        // parameters for the Brooks-Corey law
        matrixMaterialParams_.setPe(1e4);
        matrixMaterialParams_.setLambda(2);
        lensMaterialParams_.setPe(5e4);
        lensMaterialParams_.setLambda(2);
        fractureMaterialParams_.setPe(1e3);
        fractureMaterialParams_.setLambda(2);
    }

    /*!
     * \brief Function for defining the (intrinsic) permeability \f$[m^2]\f$.
     *        In this test, we use element-wise distributed permeabilities.
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
        if (scv.isOnFracture())
            return 5e-8;
        else if (isInLens_(element.geometry().center()))
            return 1e-13;
        else
            return 1e-12;
    }

    /*!
     * \brief Returns the porosity \f$[-]\f$
     *
     * \param globalPos The global position
     */
    template<class ElementSolution>
    Scalar porosity(const Element& element,
                    const SubControlVolume& scv,
                    const ElementSolution& elemSol) const
    {
        if (scv.isOnFracture())
            return 0.6;
        else
            return 0.15;
    }

    /*!
     * \brief Returns the parameter object for the Brooks-Corey material law.
     *        In this test, we use element-wise distributed material parameters.
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
        if (scv.isOnFracture())
            return fractureMaterialParams_;
        else if (isInLens_(element.geometry().center()))
            return lensMaterialParams_;
        else
            return matrixMaterialParams_;
    }

    /*!
     * \brief Function for defining which phase is to be considered as the wetting phase.
     *
     * \return the wetting phase index
     * \param globalPos The global position
     */
    template<class FluidSystem>
    int wettingPhaseAtPos(const GlobalPosition& globalPos) const
    { return FluidSystem::phase0Idx; }

    //! Updates the map of which material parameters are associated with a nodal dof.
    template<class SolutionVector>
    void updateMaterialInterfaceParams(const SolutionVector& x)
    {
        materialInterfaceParams_.update(this->fvGridGeometry(), *this, x);
    }

    //! Returns the material parameters associated with a nodal dof
    const BoxMaterialInterfaceParams<ThisType>& materialInterfaceParams() const
    { return materialInterfaceParams_; }

private:
    bool isInLens_(const GlobalPosition& globalPos) const
    {
        bool isInLens = true;
        if ( globalPos[dimWorld-1] > 0.5 || globalPos[dimWorld-1] < 0.25 )
            isInLens = false;

        for (int dir = 0; dir < dimWorld-1; ++dir)
            if (globalPos[dir] < 0.25 || globalPos[dir] > 0.75)
                isInLens = false;

        return isInLens;
    }

    MaterialLawParams matrixMaterialParams_;
    MaterialLawParams lensMaterialParams_;
    MaterialLawParams fractureMaterialParams_;

    // Determines the parameters associated with the dofs at material interfaces
    BoxMaterialInterfaceParams<ThisType> materialInterfaceParams_;

    static constexpr Scalar eps_ = 1.5e-7;
};

} // end namespace Dumux

#endif
