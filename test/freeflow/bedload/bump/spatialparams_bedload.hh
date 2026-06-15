// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup BedloadTests
 * \copydoc Dumux::BumpTestSpatialParamsBedload
 */
#ifndef DUMUX_BUMP_SPATIAL_PARAMETERS_BEDLOAD_HH
#define DUMUX_BUMP_SPATIAL_PARAMETERS_BEDLOAD_HH

#include <dumux/freeflow/spatialparams.hh>
#include <dumux/common/parameters.hh>
#include <dumux/material/sediment/bedloadtransport/bedloadformula.hh>
#include <dumux/material/sediment/bedloadtransport/grass.hh>

#include <dumux/freeflow/spatialparams.hh>

namespace Dumux {

/*!
 * \ingroup BedloadTests
 * \brief The spatial parameter class for the bump test (bedload part)
 *
 */
template<class GridGeometry, class Scalar, class VolumeVariables>
class BumpTestSpatialParamsBedload
: public FreeFlowSpatialParams<GridGeometry, Scalar, BumpTestSpatialParamsBedload<GridGeometry, Scalar, VolumeVariables>>
{
    using ThisType = BumpTestSpatialParamsBedload<GridGeometry, Scalar, VolumeVariables>;
    using ParentType = FreeFlowSpatialParams<GridGeometry, Scalar, ThisType>;
    using GridView = typename GridGeometry::GridView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    BumpTestSpatialParamsBedload(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        int nElems = this->gridGeometry().gridView().size(0);
        gravity_ = getParam<Scalar>("Problem.Gravity", 9.81);
        porosity_ = getParam<Scalar>("Sediment.Porosity");
        bedloadFormulaType_ = getParam<std::string>("Sediment.BedloadTransportFormula");
        grassAlpha_ = getParam<Scalar>("Sediment.GrassAlpha");
        fluxLimiterUpperSedimentThickness_ = getParam<Scalar>("FluxLimiterLET.UpperSedimentThickness");
        fluxLimiterLowerSedimentThickness_ = getParam<Scalar>("FluxLimiterLET.LowerSedimentThickness");

        if (bedloadFormulaType_ == "Grass" && grassAlpha_ == 0.001)
        {
            std::cout << "\nUse BedloadTransportFormula Grass " << std::endl;
            bedloadFormula_.push_back(std::make_unique<BedloadFormulaGrass<VolumeVariables>>(grassAlpha_));
        }
        else{
                DUNE_THROW(Dune::InvalidStateException, "Please set 'Sediment.BedloadTransportFormula = Grass' and 'GrassAlpha = 0.001'!");
        }

        // the following variables must be filled to use updateCoupledVariables()
        waterDepth_.assign(nElems, 0.0);
        velocityX_.assign(nElems, 0.0);
        velocityY_.assign(nElems, 0.0);
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
    * The shear stress is zero as bottom friction is not considered in this test.
    *
    * \param element The current element
    * \param scv The sub-control volume inside the element.
    *
    * \return bottom shear stress in x direction
    */
    Scalar bottomShearStressX(const Element& element,
                              const SubControlVolume& scv) const
    {
        return 0.0;
    }

    /*! \brief Get the bottom shear stress in y direction.
    *
    * The shear stress is zero as bottom friction is not considered in this test.
    *
    * \param element The current element
    * \param scv The sub-control volume inside the element.
    *
    * \return bottom shear stress in y direction
    */
    Scalar bottomShearStressY(const Element& element,
                              const SubControlVolume& scv) const
    {
        return 0.0;
    }

    /*! \brief Get the curvature radius.
    *
    *  Needed for the secondary currents approach, which is not used in this test.
    *
    * \return a dummy curvature radius
    */
    Scalar curvatureRadius(const SubControlVolume& scv) const
    {
        return 0.0;
    }

    /*!
     * \brief Return the fixed ground level of an element.
     *
     * \return fixed ground level
     */
    const Scalar fixedGroundLevel(const Element& element) const
    {
        return fixedGroundLevel_;
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
    * The friction parameter is zero as bottom friction is not considered in this test.
    *
    * \param scv The sub-control volume inside the element.
    *
    * \return friction parameter
    */
    Scalar frictionParameter(const SubControlVolume& scv) const
    {
        return 0.0;
    }

    /*! \brief Get the grain density of a certain grain class.
    *
    * \return grain density
    */
    Scalar grainDensity(int grainClassIdx) const
    {
        return grainDensity_[grainClassIdx];
    }

    /*! \brief Get the Grass parameter.
    *
    * \return Grass parameter
    */
    Scalar grassAlpha() const
    {
        return grassAlpha_;
    }

    /*! \brief Define the gravitation.
    *
    * \return gravity constant
    */
    Scalar gravity() const
    {
        return gravity_;
    }

    /*! \brief Define the gravitation.
    *
    * \param globalPos The global position
    *
    * \return gravity constant
    */
    Scalar gravity(const GlobalPosition& globalPos) const
    {
        return gravity_;
    }

    /*! \brief Get the porosity.
    *
    * \return porosity
    */
    Scalar porosity() const
    {
        return porosity_;
    }

    /*! \brief Get a dummy representative grain diameter
    *
    * This parameter is required by the BedloadFlux class for a functionality,
    * which is not used in this test.
    *
    * param grainClassIdx index of the grain class starting with 0
    *
    * \return a dummy representative grain diameter
    */
    Scalar representativeGrainDiameter(int grainClassIdx) const
    {
        return 1.0;
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
    */
    void updateCoupledVariables(std::vector<Scalar> waterDepth,
                                std::vector<Scalar> velocityX,
                                std::vector<Scalar> velocityY)
    {
        waterDepth_ = waterDepth;
        velocityX_ = velocityX;
        velocityY_ = velocityY;
    }
private:
    std::vector<std::unique_ptr<BedloadFormula<VolumeVariables>>> bedloadFormula_;
    std::string bedloadFormulaType_;
    Scalar fixedGroundLevel_ = -1.0;  // [m]
    Scalar fluxLimiterUpperSedimentThickness_;
    Scalar fluxLimiterLowerSedimentThickness_;
    std::vector<Scalar> grainDensity_{ 2650.0 };  // [kg/m^3]
    Scalar grassAlpha_;
    Scalar gravity_;
    Scalar porosity_;
    std::vector<Scalar> velocityX_;
    std::vector<Scalar> velocityY_;
    std::vector<Scalar> waterDepth_;
};

} // end namespace Dumux

#endif
