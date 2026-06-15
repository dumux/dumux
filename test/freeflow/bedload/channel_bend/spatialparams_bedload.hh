// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup BedloadTests
 * \copydoc Dumux::ChannelBendTestSpatialParamsBedload
 */
#ifndef DUMUX_CHANNEL_BEND_SPATIAL_PARAMETERS_BEDLOAD_HH
#define DUMUX_CHANNEL_BEND_SPATIAL_PARAMETERS_BEDLOAD_HH

#include <dumux/freeflow/spatialparams.hh>
#include <dumux/common/parameters.hh>
#include <dumux/material/sediment/bedloadtransport/bedloadformula.hh>
#include <dumux/material/sediment/bedloadtransport/meyerpetermueller.hh>
#include <dumux/material/sediment/bedloadtransport/hidingexposure.hh>
#include <dumux/material/sediment/bedloadtransport/criticalshieldsstress.hh>

#include <dumux/freeflow/spatialparams.hh>

namespace Dumux {

/*!
 * \ingroup BedloadTests
 * \brief The spatial parameter class for the channel bend test (bedload part)
 *
 */
template<class GridGeometry, class Scalar, class VolumeVariables>
class ChannelBendTestSpatialParamsBedload
: public FreeFlowSpatialParams<GridGeometry, Scalar,
                               ChannelBendTestSpatialParamsBedload<GridGeometry, Scalar, VolumeVariables>>
{
    using ThisType = ChannelBendTestSpatialParamsBedload<GridGeometry, Scalar, VolumeVariables>;
    using ParentType = FreeFlowSpatialParams<GridGeometry, Scalar, ThisType>;
    using GridView = typename GridGeometry::GridView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    ChannelBendTestSpatialParamsBedload(std::shared_ptr<const GridGeometry> gridGeometry,
                                        std::vector<Scalar>& initialBedSurface,
                                        std::vector<Scalar>& initialWaterDepth,
                                        std::vector<Scalar>& initialVelocityX,
                                        std::vector<Scalar>& initialVelocityY,
                                        std::vector<Scalar>& curvatureRadius)
        : ParentType(gridGeometry)
    {
        int nElems = this->gridGeometry().gridView().size(0);
        initialBedSurface_ = initialBedSurface;
        waterDepth_ = initialWaterDepth;
        velocityX_ = initialVelocityX;
        velocityY_ = initialVelocityY;
        curvatureRadius_ = curvatureRadius;
        // use updateCoupledVariables() to fill the bottom shear stress with the actual values
        bottomShearStressX_.assign(nElems, 0.0);
        bottomShearStressY_.assign(nElems, 0.0);

        // read parameter
        int nGrainClasses = getParam<int>("Sediment.NumberGrainClasses");
        gravity_ = getParam<Scalar>("Problem.Gravity", 9.81);
        porosity_ = getParam<Scalar>("Sediment.Porosity");
        bedloadFormulaType_ = getParam<std::string>("Sediment.BedloadTransportFormula");
        meyerPeterMuellerCoefficient_ = getParam<Scalar>("Sediment.MeyerPeterMuellerCoefficient");
        fluxLimiterUpperSedimentThickness_ = getParam<Scalar>("FluxLimiterLET.UpperSedimentThickness");
        fluxLimiterLowerSedimentThickness_ = getParam<Scalar>("FluxLimiterLET.LowerSedimentThickness");

        // calculate the critical Shields stress
        criticalShieldsStress_.resize(nGrainClasses);
        for (int i = 0; i < nGrainClasses; i++) {
            Scalar dimensionlessGrainDiameter = calculateDimensionlessGrainDiameter(representativeGrainDiameter_[i],
                                                                                    gravity_,
                                                                                    kinematicViscosity_,
                                                                                    grainDensity_ / waterDensity_);
            criticalShieldsStress_[i] = calculateCriticalShieldsStress(dimensionlessGrainDiameter);
        }

        // calculate the hiding and exposure coefficient following Karim, Holly, Yang (1987)
        // "Ialluvial Numerical Simulation of Mobile-Bed Rivers: Part I.Theoretical and Numerical Principles."
        Scalar medianGrainDiameter;
        if (representativeGrainDiameter_.size() % 2 == 1)
        {
            medianGrainDiameter = representativeGrainDiameter_[representativeGrainDiameter_.size() / 2];
        }
        else
        {
            medianGrainDiameter = 0.5 * (representativeGrainDiameter_[representativeGrainDiameter_.size() / 2 - 1]
                                         + representativeGrainDiameter_[representativeGrainDiameter_.size() / 2]);
        }
        auto hidingExposureCoefficientsKarimHollyYang = calculateHidingExposureCoefficientsKarimHollyYang(representativeGrainDiameter_,
                                                                                                          medianGrainDiameter);

        // set the bedload transport formula
        if (bedloadFormulaType_ == "MeyerPeterMueller" && meyerPeterMuellerCoefficient_ == 8.0)
        {
            std::cout << "\nUse BedloadTransportFormula MeyerPeterMueller " << std::endl;
            for (int i = 0; i < nGrainClasses; i++)
            {
                bedloadFormula_.push_back(std::make_unique<BedloadFormulaMeyerPeterMueller<VolumeVariables>>(criticalShieldsStress_[i],
                                                                                                             grainDensity_,
                                                                                                             gravity_,
                                                                                                             meyerPeterMuellerCoefficient_,
                                                                                                             representativeGrainDiameter_[i],
                                                                                                             waterDensity_,
                                                                                                             hidingExposureCoefficientsKarimHollyYang[i]));
            }
        }
        else {
            DUNE_THROW(Dune::InvalidStateException, "Please set 'Sediment.BedloadTransportFormula = MeyerPeterMueller' and 'MeyerPeterMuellerCoefficient = 8.0'!");
        }
    }

    /*! \brief Get the bedload formula of a certain grain class.
    *
    * \param grainClassIdx the index of the grain class (starting with 0)
    *
    * \return bedload formula
    */
    const BedloadFormula<VolumeVariables>& bedloadFormula(const Scalar grainClassIdx) const
    {
        return *bedloadFormula_[grainClassIdx];
    }

    /*! \brief Get the bottom shear stress in x direction.
    *
    * \param element The current element
    * \param scv The sub-control volume inside the element.
    *
    * \return bottom shear stress in x direction
    */
    Scalar bottomShearStressX(const Element& element,
                              const SubControlVolume& scv) const
    {
        auto eIdx = scv.elementIndex();
        return bottomShearStressX_[eIdx];
    }

    /*! \brief Get the bottom shear stress in y direction.
    *
    * \param element The current element
    * \param scv The sub-control volume inside the element.
    *
    * \return bottom shear stress in y direction
    */
    Scalar bottomShearStressY(const Element& element,
                              const SubControlVolume& scv) const
    {
        auto eIdx = scv.elementIndex();
        return bottomShearStressY_[eIdx];
    }

    /*! \brief Get the curvature radius.
    *
    * Needed for the secondary currents approach.
    *
    * \param element The current element
    * \param scv The sub-control volume inside the element.
    *
    * \return curvature radius
    */
    Scalar curvatureRadius(const SubControlVolume& scv) const
    {
        auto eIdx = scv.elementIndex();
        return curvatureRadius_[eIdx];
    }

    /*! \brief Get the fluxLimiterLowerSedimentThickness
    *
    * \return fluxLimiterLowerSedimentThickness
    */
    Scalar fluxLimiterLowerSedimentThickness() const
    {
        return fluxLimiterLowerSedimentThickness_;
    }

    /*! \brief Get the fluxLimiterUpperSedimentThickness
    *
    * \return fluxLimiterUpperSedimentThickness
    */
    Scalar fluxLimiterUpperSedimentThickness() const
    {
        return fluxLimiterUpperSedimentThickness_;
    }

    /*! \brief Get the friction parameter.
    *
    * \param scv The sub-control volume inside the element.
    *
    * \return friction parameter
    */
    Scalar frictionParameter(const SubControlVolume& scv) const
    {
        return manningFrictionCoefficient_;
    }

    /*! \brief Get the gravitation.
    *
    * \return gravity constant
    */
    Scalar gravity() const
    {
        return gravity_;
    }

    /*! \brief Get the grain density of a certain grain class.
    *
    * \param grainClassIdx the index of the grain class (starting with 0)
    *
    * \return grain density
    */
    Scalar grainDensity(int grainClassIdx) const
    {
        return grainDensity_;
    }

    /*!
    * \brief Return the fixed ground level of an element.
    *
    * \param element The current element
    * \param scv The sub-control volume inside the element.
    *
    * \return fixed ground level
    */
    const Scalar fixedGroundLevel(const Element& element,
                                  const SubControlVolume& scv) const
    {
        return fixedGroundLevel_;
    }

    /*! \brief Get the initial bed surface.
    *
    * \param element The current element
    * \param scv The sub-control volume inside the element.
    *
    * \return initial bed surface
    */
    Scalar initialBedSurface(const Element& element,
                             const SubControlVolume& scv) const
    {
        auto eIdx = scv.elementIndex();
        return initialBedSurface_[eIdx];
    }

    /*! \brief Get the initial mass fraction for a certain grain class and layer
    *
    * \param element The current element
    * \param scv The sub-control volume inside the element
    * \param layerId The ID of the layer
    * \param grainClassIdx The index of the grain class (starting with zero)
    *
    * \return initial mass fraction
    */
    Scalar initialMassFraction(const Element& element,
                               const SubControlVolume& scv,
                               const int layerId,
                               const int grainClassIdx) const
    {
        // in this test all layers have the same initial mass fractions
        return initialMassFraction_[grainClassIdx];
    }

    /*! \brief Get the initial upper layer limit.
    *
    * This test uses only two layer. Therefore, only the upper limit of the lower layer
    * must be specified, as the upper limit of the upper layer is the bed surface.
    *
    * \param element The current element
    * \param scv The sub-control volume inside the element
    * \param layerId The ID of the layer
    *
    * \return initial upper layer limit
    */
    Scalar initialUpperLayerLimit(const Element& element,
                                  const SubControlVolume& scv,
                                  const int layerId) const
    {
        // the upper layer limit is constant throughout the whole domain
        return initialUpperLayerLimit_.at(layerId);
    }

    /*! \brief Get the porosity.
    *
    * \return porosity
    */
    Scalar porosity() const
    {
        return porosity_;
    }

    /*! \brief Get the representative grain diameter of a certain grain class
    *
    * \param grainClassIdx index of the grain class starting with 0
    *
    * \return representative grain diameter
    */
    Scalar representativeGrainDiameter(int grainClassIdx) const
    {
        return representativeGrainDiameter_[grainClassIdx];
    }
    /*! \brief Get the velocity in x direction.
    *
    * \param element The current element
    * \param scv The sub-control volume inside the element.
    *
    * \return velocity in x direction
    */
    Scalar velocityX(const Element& element,
                     const SubControlVolume& scv) const
    {
        auto eIdx = scv.elementIndex();
        return velocityX_[eIdx];
    }

    /*! \brief Get the velocity in y direction.
    *
    * \param element The current element
    * \param scv The sub-control volume inside the element.
    *
    * \return velocity in y direction
    */
    Scalar velocityY(const Element& element,
                     const SubControlVolume& scv) const
    {
        auto eIdx = scv.elementIndex();
        return velocityY_[eIdx];
    }

    /*! \brief Get the water depth.
    *
    * \param element The current element
    * \param scv The sub-control volume inside the element.
    *
    * \return water depth
    */
    Scalar waterDepth(const Element& element,
                      const SubControlVolume& scv) const
    {
        auto eIdx = scv.elementIndex();
        return waterDepth_[eIdx];
    }

    /*! \brief Update the coupled variables.
    *
    * Update the water depth, both velocity components and the bottom shear stress.
    *
    * \param std::vector<Scalar> waterDepth
    * \param std::vector<Scalar> velocityX
    * \param std::vector<Scalar> velocityY
    * \param std::vector<Scalar> bottomShearStressX
    * \param std::vector<Scalar> bottomShearStressY
    */
    void updateCoupledVariables(std::vector<Scalar> waterDepth,
                                std::vector<Scalar> velocityX,
                                std::vector<Scalar> velocityY,
                                std::vector<Scalar> bottomShearStressX,
                                std::vector<Scalar> bottomShearStressY)
    {
        waterDepth_ = waterDepth;
        velocityX_ = velocityX;
        velocityY_ = velocityY;
        bottomShearStressX_ = bottomShearStressX;
        bottomShearStressY_ = bottomShearStressY;
    }
private:
    std::string bedloadFormulaType_;
    std::vector<std::unique_ptr<BedloadFormula<VolumeVariables>>> bedloadFormula_;
    std::vector<Scalar> bottomShearStressX_;
    std::vector<Scalar> bottomShearStressY_;
    std::vector<Scalar> criticalShieldsStress_;
    std::vector<Scalar> curvatureRadius_;
    Scalar fixedGroundLevel_ = 0.0;  // [m]
    Scalar fluxLimiterUpperSedimentThickness_;
    Scalar fluxLimiterLowerSedimentThickness_;
    Scalar grainDensity_ = 2650.0; // [kg/m^3]
    Scalar gravity_;
    std::vector<Scalar> initialBedSurface_;
    std::vector<Scalar> initialMassFraction_ = { 0.33, 0.33, 0.34 };  // [-]
    // we have two layers. The upper limit of layer 1 is the bed surface,
    // therefore we have to specify only the upper limit of layer 2
    std::map<int, Scalar> initialUpperLayerLimit_ = { {2, 9.0} };  // [m]
    Scalar kinematicViscosity_ = 1e-6;  // [m^2/s]
    Scalar manningFrictionCoefficient_ = 0.0167; // The same value like in the shallow water part [s/m^(1/3)]
    Scalar meyerPeterMuellerCoefficient_;
    Scalar porosity_;
    std::vector<Scalar> representativeGrainDiameter_ {0.00038, 0.00097, 0.00246}; // [m]
    std::vector<Scalar> velocityX_;
    std::vector<Scalar> velocityY_;
    Scalar waterDensity_ = 1000.0;  // [kg/m^3]
    std::vector<Scalar> waterDepth_;
};

} // end namespace Dumux

#endif
