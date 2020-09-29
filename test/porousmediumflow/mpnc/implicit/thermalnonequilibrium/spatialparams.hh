/*****************************************************************************
 *   Institute for Modelling Hydraulic and Environmental Systems             *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup MPNCTests
 * \brief Spatialparameters for the combustionproblem1c.
 *
 * Parameters for the actual simulation domain and an outflow region are provided.
 */

#ifndef DUMUX_COMBUSTION_SPATIALPARAMS_HH
#define DUMUX_COMBUSTION_SPATIALPARAMS_HH

#include <dune/common/parametertreeparser.hh>

#include <dumux/material/spatialparams/fvnonequilibrium.hh>
#include <dumux/material/fluidmatrixinteractions/2p/heatpipelaw.hh>
#include <dumux/material/fluidmatrixinteractions/1pia/fluidsolidinterfacialareashiwang.hh>
#include <dumux/porousmediumflow/properties.hh>
#include <dumux/material/spatialparams/fv.hh>

namespace Dumux {

/**
 * \brief Definition of the spatial parameters for the one component combustion problem
 */
template<class GridGeometry, class Scalar>
class CombustionSpatialParams
: public FVNonEquilibriumSpatialParams<GridGeometry, Scalar,
                                       CombustionSpatialParams<GridGeometry, Scalar>>
{
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;
    using ThisType = CombustionSpatialParams<GridGeometry, Scalar>;
    using ParentType = FVNonEquilibriumSpatialParams<GridGeometry, Scalar, ThisType>;

    enum {dimWorld = GridView::dimensionworld};
    using GlobalPosition = typename SubControlVolume::GlobalPosition;

public:
    //! Export the type used for the permeability
    using PermeabilityType = Scalar;

    using FluidSolidInterfacialAreaFormulation = FluidSolidInterfacialAreaShiWang<Scalar>;

    using FluidMatrixInteraction = FluidMatrix::HeatPipeLaw<Scalar>;

    CombustionSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry) : ParentType(gridGeometry)
    {
        // this is the parameter value from file part
        porosity_ = getParam<Scalar>("SpatialParams.PorousMedium.porosity");
        intrinsicPermeabilityOutFlow_ = getParam<Scalar>("SpatialParams.Outflow.permeabilityOutFlow");
        porosityOutFlow_                = getParam<Scalar>("SpatialParams.Outflow.porosityOutFlow");
        interfacialTension_  = getParam<Scalar>("Constants.interfacialTension");

        characteristicLength_ =getParam<Scalar>("SpatialParams.PorousMedium.meanPoreSize");

        using std::pow;
        intrinsicPermeability_  =  (pow(characteristicLength_,2.0)  * pow(porosity_, 3.0)) / (150.0 * pow((1.0-porosity_),2.0)); // 1.69e-10 ; //

        factorEnergyTransfer_ = getParam<Scalar>("SpatialParams.PorousMedium.factorEnergyTransfer");
        lengthPM_ = getParam<Scalar>("Grid.lengthPM");

        using std::sqrt;
        typename FluidMatrixInteraction::Params params;
        params.p0 = sqrt(porosity_/intrinsicPermeability_);
        params.gamma = interfacialTension_; // interfacial tension of water-air at 100°C

        typename FluidMatrixInteraction::EffToAbsParams effToAbsParams;
        effToAbsParams.swr = getParam<Scalar>("SpatialParams.soil.Swr");
        effToAbsParams.snr = getParam<Scalar>("SpatialParams.soil.Snr");

        fluidMatrixInteraction_ = std::make_unique<FluidMatrixInteraction>(params, effToAbsParams);
    }

    template<class ElementSolution>
    PermeabilityType permeability(const Element& element,
                                  const SubControlVolume& scv,
                                  const ElementSolution& elemSol) const
    {
        const auto& globalPos =  scv.dofPosition();
        if ( inOutFlow(globalPos) )
            return intrinsicPermeabilityOutFlow_ ;
        else
            return intrinsicPermeability_ ;
    }

    /*!
     * \brief Function for defining the porosity which is possibly solution dependent.
     *
     * \param element The current element
     * \param scv The sub-control volume inside the element.
     * \param elemSol The solution at the dofs connected to the element.
     * \return The porosity
     */
    template<class ElementSolution>
    Scalar porosity(const Element& element,
                    const SubControlVolume& scv,
                    const ElementSolution& elemSol) const
    {
        if (inOutFlow(scv.dofPosition()))
            return porosityOutFlow_;
        else
            return porosity_;
    }

    /*!
     * \brief Function for defining the porosity which is possibly solution dependent.
     *
     * \param element The current element
     * \param scv The sub-control volume inside the element.
     * \param elemSol The solution at the dofs connected to the element.
     * \param compIdx The index of the component
     * \return The porosity
     */
    template<class SolidSystem, class ElementSolution>
    Scalar inertVolumeFraction(const Element& element,
                               const SubControlVolume& scv,
                               const ElementSolution& elemSol,
                               int compIdx) const
    {
        if (compIdx == SolidSystem::comp0Idx)
        {
            if (inOutFlow(scv.dofPosition()))
                return 1.0-porosityOutFlow_;
            else
                return 0;
        }
        else
        {
            if (inOutFlow(scv.dofPosition()))
                return 0;
            else
                return 1.0-porosity_;
        }
    }

    /*!
     * \brief Function for defining which phase is to be considered as the wetting phase.
     *
     * \param globalPos The global position
     * \return The wetting phase index
     */
    template<class FluidSystem>
    int wettingPhaseAtPos(const GlobalPosition& globalPos) const
    {
        return FluidSystem::phase0Idx;
    }

    /*!
     * \brief Returns the characteristic length for the mass transfer.
     * \param globalPos The position in global coordinates.
     */
    const Scalar characteristicLengthAtPos(const  GlobalPosition & globalPos) const
    { return characteristicLength_ ; }

    /*!
     * \brief Returns the pre factor the the energy transfer
     * \param globalPos The position in global coordinates.
     */
    const Scalar factorEnergyTransferAtPos(const  GlobalPosition & globalPos) const
    { return factorEnergyTransfer_; }

    //! Returns if the tested position is at the right end of the porous medium.
    bool inOutFlow(const GlobalPosition & globalPos) const { return globalPos[0] > (lengthPM_ - eps_) ;    }
    //! Returns the length of the porous medium domain
    Scalar lengthPM() const { return lengthPM_ ; }
    //! Returns the interfacial tension
    Scalar interfacialTension() const { return interfacialTension_ ; }

    /*!
     * \brief Returns the parameters for the material law at a given location
     * \param globalPos A global coordinate vector
     */
    const FluidMatrixInteraction& fluidMatrixInteractionAtPos(const GlobalPosition &globalPos) const
    {
        return *fluidMatrixInteraction_;
    }

private:
    static constexpr Scalar eps_ = 1e-6;

    // Porous Medium Domain
    Scalar intrinsicPermeability_ ;
    Scalar porosity_ ;
    Scalar factorEnergyTransfer_ ;
    Scalar characteristicLength_ ;

    // Outflow Domain
    Scalar intrinsicPermeabilityOutFlow_ ;
    Scalar porosityOutFlow_ ;

    // solid parameters
    Scalar interfacialTension_ ;

    // grid
    Scalar lengthPM_ ;

    std::unique_ptr<FluidMatrixInteraction> fluidMatrixInteraction_;
};

}

#endif
