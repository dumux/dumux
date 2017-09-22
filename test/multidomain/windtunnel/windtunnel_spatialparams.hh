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
 *
 * \brief The spatial parameters class for the windtunnel test problem (taken from Fetzer2017c)
 */
#ifndef DUMUX_WINDTUNNELSPATIALPARAMETERS_HH
#define DUMUX_WINDTUNNELSPATIALPARAMETERS_HH

#include <dumux/material/spatialparams/implicit.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedvangenuchten.hh>
#include <dumux/io/plotmateriallaw.hh>
#include <dumux/io/plotthermalconductivitymodel.hh>

// TODO !!
//#include <appl/staggeredgrid/multidomain/navierstokes2ctdarcy2p2ct/properties.hh>

namespace Dumux
{

//forward declaration
template<class TypeTag>
class WindTunnelSpatialParams;

namespace Properties // TODO check properties (except "NEW_TYPE_TAG...")
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(WindTunnelSpatialParams);

// Set the spatial parameters
SET_TYPE_PROP(DarcySpatialParams, SpatialParams, WindTunnelSpatialParams<TypeTag>);

// Set the material Law
NEW_PROP_TAG(EffMaterialLaw);
SET_PROP(DarcySpatialParams, EffMaterialLaw)
{
private:
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
public:
//     typedef VanGenuchten<Scalar> type;
//     typedef RegularizedBrooksCorey<Scalar> type;
    typedef RegularizedVanGenuchten<Scalar> type;
//     typedef RegularizedVanGenuchten<Scalar, LinearizedRegVanGenuchtenParams<Scalar, TypeTag> > type;
};

// Set the material Law
SET_TYPE_PROP(DarcySpatialParams, MaterialLaw, EffToAbsLaw<typename GET_PROP_TYPE(TypeTag, EffMaterialLaw)>);

// Decide whether to plot the porous-medium properties or not
NEW_PROP_TAG(OutputPlotMaterialLaw);
SET_BOOL_PROP(DarcySpatialParams, OutputPlotMaterialLaw, true);

// Set properties of the porous medium
NEW_PROP_TAG(SpatialParamsRandomField);
NEW_PROP_TAG(SpatialParamsVerticalLayers);
NEW_PROP_TAG(Soil1RegularizationLowSw);
NEW_PROP_TAG(Soil1RegularizationHighSw);
NEW_PROP_TAG(Soil2RegularizationLowSw);
NEW_PROP_TAG(Soil2RegularizationHighSw);
SET_BOOL_PROP(DarcySpatialParams, SpatialParamsRandomField, false);
SET_BOOL_PROP(DarcySpatialParams, SpatialParamsVerticalLayers, false);
SET_SCALAR_PROP(DarcySpatialParams, Soil1RegularizationLowSw, 0.0);
SET_SCALAR_PROP(DarcySpatialParams, Soil1RegularizationHighSw, 1.0);
SET_SCALAR_PROP(DarcySpatialParams, Soil2RegularizationLowSw, 0.0);
SET_SCALAR_PROP(DarcySpatialParams, Soil2RegularizationHighSw, 1.0);
}

/*!
 * \ingroup TwoPTwoCModel
 * \ingroup ImplicitTestProblems
 * \brief Definition of the spatial parameters
 */
template<class TypeTag>
class WindTunnelSpatialParams : public ImplicitSpatialParams<TypeTag>
{
    using ParentType = ImplicitSpatialParams<TypeTag>;

    using Grid = typename GET_PROP_TYPE(TypeTag, Grid);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using CoordScalar = typename GridView::ctype;
    enum { dimWorld=GridView::dimensionworld };

    using GlobalPosition = Dune::FieldVector<Scalar,dimWorld>;
    using DimVector = Dune::FieldVector<CoordScalar,dimWorld>;
    using IntVector =  Dune::FieldVector<int,dimWorld>;
    using ScalarVector = std::vector<Scalar>;

    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using Element = typename GridView::template Codim<0>::Entity;

    using EffMaterialLaw = typename GET_PROP_TYPE(TypeTag, EffMaterialLaw);
    using MaterialLaw = EffToAbsLaw<EffMaterialLaw>;
    using MaterialLawParams = typename MaterialLaw::Params;

public:
    // export permeability type
    using PermeabilityType = Scalar;

    /*!
     * \brief The constructor
     *
     * \param gridView The grid view
     */
    WindTunnelSpatialParams(const GridView &gridView)
        : ParentType(gridView)
    {
        // domain extents
        Scalar noDarcyX1 = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, NoDarcyX1);
        Scalar noDarcyX2 = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, NoDarcyX2);
        std::vector<Scalar> positions0 = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::vector<Scalar>, Grid, Positions0);
        std::vector<Scalar> positions1 = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::vector<Scalar>, Grid, Positions1);

        bBoxMin_[0] = std::max(positions0.front(),noDarcyX1);
        bBoxMax_[0] = std::min(positions0.back(),noDarcyX2);
        bBoxMin_[1] = positions1.front();
        bBoxMax_[1] = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Grid, InterfacePosY);
        lengthPM_ = bBoxMax_[0] - bBoxMin_[0];
        heightPM_ = bBoxMax_[1] - bBoxMin_[1];

        // soil properties
        plotMaterialLaw_ = GET_PARAM_FROM_GROUP(TypeTag, bool, Output, PlotMaterialLaw);
        randomField_ = GET_PARAM_FROM_GROUP(TypeTag, bool, SpatialParams, RandomField);
        verticalLayers_ = GET_PARAM_FROM_GROUP(TypeTag, bool, SpatialParams, VerticalLayers);
        numberOfSoilLayers_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, SpatialParams, NumberOfSoilLayers);
        firstSoilLayerIdx_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, SpatialParams, SoilLayerIdx1);
        secondSoilLayerIdx_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, int, SpatialParams, SoilLayerIdx2);

        if (firstSoilLayerIdx_ == 1 || secondSoilLayerIdx_ == 1)
        {
            permeability1_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Soil1, Permeability);
            porosity1_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Soil1, Porosity);
            alphaBJ1_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Soil1, AlphaBJ);
            solidThermalConductivity1_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Soil1, ThermalConductivitySolid);
            spatialParams1_.setSwr(GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Soil1, Swr));
            spatialParams1_.setSnr(GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Soil1, Snr));
            spatialParams1_.setVgAlpha(GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Soil1, VgAlpha));
            spatialParams1_.setVgn(GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Soil1, VgN));
            spatialParams1_.setPcLowSw(std::max(0.01,GET_PARAM_FROM_GROUP(TypeTag, Scalar, Soil1, RegularizationLowSw)));
            spatialParams1_.setPcHighSw(std::min(0.99,GET_PARAM_FROM_GROUP(TypeTag, Scalar, Soil1, RegularizationHighSw)));
        }

        if (firstSoilLayerIdx_ == 2 || secondSoilLayerIdx_ == 2)
        {
            permeability2_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Soil2, Permeability);
            porosity2_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Soil2, Porosity);
            alphaBJ2_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Soil2, AlphaBJ);
            solidThermalConductivity2_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Soil2, ThermalConductivitySolid);
            spatialParams2_.setSwr(GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Soil2, Swr));
            spatialParams2_.setSnr(GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Soil2, Snr));
            spatialParams2_.setVgAlpha(GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Soil2, VgAlpha));
            spatialParams2_.setVgn(GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, Scalar, Soil2, VgN));
            spatialParams2_.setPcLowSw(std::max(0.01,GET_PARAM_FROM_GROUP(TypeTag, Scalar, Soil2, RegularizationLowSw)));
            spatialParams2_.setPcHighSw(std::min(0.99,GET_PARAM_FROM_GROUP(TypeTag, Scalar, Soil2, RegularizationHighSw)));
        }
    }

    /*!
     * \brief Function for defining the intrinsic (absolute) permeability.
     *
     * \return intrinsic (absolute) permeability
     * \param globalPos The position of the center of the element
     */
    const Scalar intrinsicPermeabilityAtPos(const GlobalPosition& globalPos) const
    {
        int curSoilType = soilType(globalPos);
        if (curSoilType == 0)
            DUNE_THROW(Dune::NotImplemented, "This soil type is not implemented");//randomPermeability_[indexSet_.index(element.template subEntity<dim> (scvIdx))];
        else if (curSoilType == 1)
            return permeability1_;
        else if (curSoilType == 2)
            return permeability2_;
        else
            DUNE_THROW(Dune::NotImplemented, "This soil type is not implemented");
    }

    /*!
     * \brief Apply the intrinsic permeability tensor to a pressure
     *        potential gradient.
     *
     * \param element The current finite element
     * \param fvGeometry The current finite volume geometry of the element
     * \param scvIdx The index of the sub-control volume
     */
    const Scalar intrinsicPermeability(const Element &element,
                                       const FVElementGeometry &fvGeometry,
                                       const int scvIdx) const
    {
        return intrinsicPermeabilityAtPos(element.geometry().center());
    }

    /*!
     * \brief Define the porosity \f$[-]\f$ of the spatial parameters
     *
     * \param globalPos The position of the center of the element
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    {
        int curSoilType = soilType(globalPos);
        if (curSoilType == 0)
            DUNE_THROW(Dune::NotImplemented, "This soil type is not implemented");//randomPorosity_[indexSet_.index(element.template subEntity<dim> (scvIdx))];
        else if (curSoilType == 1)
            return porosity1_;
        else if (curSoilType == 2)
            return porosity2_;
        else
            DUNE_THROW(Dune::NotImplemented, "This soil type is not implemented");
    }

    /*!
     * \brief Define the porosity \f$[-]\f$ of the spatial parameters
     *
     * \param element The finite element
     * \param fvGeometry The finite volume geometry
     * \param scvIdx The local index of the sub-control volume where
     *                    the porosity needs to be defined
     */
    Scalar porosity(const Element &element,
                    const FVElementGeometry &fvGeometry,
                    const int scvIdx) const
    {
        return porosityAtPos(element.geometry().center());
    }

    /*!
     * \brief return the parameter object for the material law
     *
     * \param globalPos The position of the center of the element
     */
    const MaterialLawParams& materialLawParamsAtPos(const GlobalPosition& globalPos) const
    {
        int curSoilType = soilType(globalPos);
        if (curSoilType == 0)
            DUNE_THROW(Dune::NotImplemented, "This soil type is not implemented");//randomSpatialParams_[indexSet_.index(element.template subEntity<dim> (scvIdx))];
        else if (curSoilType == 1)
            return spatialParams1_;
        else if (curSoilType == 2)
            return spatialParams2_;
        else
            DUNE_THROW(Dune::NotImplemented, "This soil type is not implemented");
    }

    /*!
     * \brief return the parameter object for the material law
     *
     * \param element The current finite element
     * \param fvGeometry The current finite volume geometry of the element
     * \param scvIdx The index of the sub-control volume
     */
    const MaterialLawParams& materialLawParams(const Element &element,
                                               const FVElementGeometry &fvGeometry,
                                               const int scvIdx) const
    {
        return materialLawParamsAtPos(element.geometry().center());
    }

    /*!
     * \brief Returns the heat capacity \f$[J / (kg K)]\f$ of the rock matrix.
     *
     * This is only required for non-isothermal models.
     *
     * \param element The finite element
     * \param fvGeometry The finite volume geometry
     * \param scvIdx The local index of the sub-control volume
     */
    Scalar solidHeatCapacity(const Element &element,
                             const FVElementGeometry &fvGeometry,
                             const int scvIdx) const
    {
        return 790;
    }

    /*!
     * \brief Returns the mass density \f$[kg / m^3]\f$ of the rock matrix.
     *
     * This is only required for non-isothermal models.
     *
     * \param globalPos The position of the center of the element
     */
    Scalar solidDensityAtPos(const GlobalPosition& globalPos) const
    {
        return 2700;
    }

    /*!
     * \brief Returns the mass density \f$[kg / m^3]\f$ of the rock matrix.
     *
     * This is only required for non-isothermal models.
     *
     * \param element The finite element
     * \param fvGeometry The finite volume geometry
     * \param scvIdx The local index of the sub-control volume
     */
    Scalar solidDensity(const Element &element,
                        const FVElementGeometry &fvGeometry,
                        const int scvIdx) const
    {
        return solidDensityAtPos(element.geometry().center());
    }

    /*!
     * \brief Returns the thermal conductivity \f$\mathrm{[W/(m K)]}\f$ of the porous material.
     *
     * \param globalPos The position of the center of the element
     */
    Scalar solidThermalConductivityAtPos(const GlobalPosition& globalPos) const
    {
        int curSoilType = soilType(globalPos);
        if (curSoilType == 0)
            DUNE_THROW(Dune::NotImplemented, "This soil type is not implemented");//randomSolidThermalConductivity_[indexSet_.index(element.template subEntity<dim> (scvIdx))];
        else if (curSoilType == 1)
            return solidThermalConductivity1_;
        else if (curSoilType == 2)
            return solidThermalConductivity2_;
        else
            DUNE_THROW(Dune::NotImplemented, "This soil type is not implemented");
    }

    /*!
     * \brief Returns the thermal conductivity \f$\mathrm{[W/(m K)]}\f$ of the porous material.
     *
     * \param element The finite element
     * \param fvGeometry The finite volume geometry
     * \param scvIdx The local index of the sub-control volume where
     *                    the heat capacity needs to be defined
     */
    Scalar solidThermalConductivity(const Element &element,
                                    const FVElementGeometry &fvGeometry,
                                    const int scvIdx) const
    {
        return solidThermalConductivityAtPos(element.geometry().center());
    }

    /*!
     * \brief Evaluate the Beavers-Joseph coefficient at given position
     *
     * \param globalPos The global position
     *
     * \return Beavers-Joseph coefficient
     */
    Scalar beaversJosephCoeffAtPos(const GlobalPosition &globalPos) const
    {
        int curSoilType = soilType(globalPos);
        if (curSoilType == 0)
            return 1.0;
        else if (curSoilType == 1)
            return alphaBJ1_;
        else if (curSoilType == 2)
            return alphaBJ2_;
        else
            DUNE_THROW(Dune::NotImplemented, "This soil type is not implemented");
    }

    /*!
     * \brief Returns the relative position between 0 and 1 of a global coordinate
     */
    GlobalPosition relativePosition(const GlobalPosition &globalPos) const
    {
        GlobalPosition relativePosition(0.0);
        relativePosition[0] = (globalPos[0] - bBoxMin_[0]) / lengthPM_;
        relativePosition[1] = (globalPos[1] - bBoxMin_[1]) / heightPM_;
        return relativePosition;
    }

    /*!
     * \brief Returns the layerIdx of a global coordinate
     */
    IntVector layerIdx(const GlobalPosition &globalPos) const
    {
        GlobalPosition curRelativePosition = relativePosition(globalPos);
        curRelativePosition *= numberOfSoilLayers_ * 1.0;

        IntVector layerIdx(0);
        layerIdx[0] = int(floor(curRelativePosition[0]));
        layerIdx[0] = std::max(layerIdx[0], 0);
        layerIdx[0] = std::min(layerIdx[0], numberOfSoilLayers_-1);
        layerIdx[1] = int(floor(curRelativePosition[1]));
        layerIdx[1] = std::max(layerIdx[1], 0);
        layerIdx[1] = std::min(layerIdx[1], numberOfSoilLayers_-1);
        return layerIdx;
    }

    /*!
     * \brief Returns global position of the heterogeneity locations
     */
    std::vector<Scalar> heterogeneityPos() const
    {
        std::vector<Scalar> heterogeneityPos(numberOfSoilLayers_-1, bBoxMin_[0]);
        heterogeneityPos[0] = bBoxMin_[0] + lengthPM_ / numberOfSoilLayers_;
        for (int idx = 1; idx < heterogeneityPos.size(); idx++)
        {
            heterogeneityPos[idx] = heterogeneityPos[idx-1] + lengthPM_ / numberOfSoilLayers_;
        }
        return heterogeneityPos;
    }

    /*!
     * \brief Returns the the soil typ at a global coordinate
     */
    int soilType(const GlobalPosition &globalPos) const
    {
        if (randomField_)
            return 0;

        int curLayerID = layerIdx(globalPos)[0];
        if (verticalLayers_)
            curLayerID = layerIdx(globalPos)[1];
        if (curLayerID % 2 == 0)
            return firstSoilLayerIdx_;
        return secondSoilLayerIdx_;
    }

    /*!
     * \brief This is called from the coupled problem and creates
     *        a gnuplot output of the Pc-Sw curve
     */
    void plotMaterialLaw()
    {
        if (plotMaterialLaw_ && !randomField_)
        {
            GnuplotInterface<Scalar> gnuplot;
            PlotMaterialLaw<TypeTag> plotMaterialLaw;
            PlotThermalConductivityModel<TypeTag> plotThermalConductivityModel(293.0, 1e5);

            if (firstSoilLayerIdx_ == 1 || secondSoilLayerIdx_ == 1)
                plotMaterialLaw.addpcswcurve(gnuplot, spatialParams1_, 0.0, 1.0, "pcsw_soil1");
            if (firstSoilLayerIdx_ == 2 || secondSoilLayerIdx_ == 2)
                plotMaterialLaw.addpcswcurve(gnuplot, spatialParams2_, 0.0, 1.0, "pcsw_soil2");
            gnuplot.plot("pc-Sw");

            gnuplot.resetAll();
            if (firstSoilLayerIdx_ == 1 || secondSoilLayerIdx_ == 1)
                plotThermalConductivityModel.addlambdaeffcurve(gnuplot, porosity1_, 2700.0, solidThermalConductivity1_, 0.0, 1.0, "lambda_eff_soil1");
            if (firstSoilLayerIdx_ == 2 || secondSoilLayerIdx_ == 2)
                plotThermalConductivityModel.addlambdaeffcurve(gnuplot, porosity2_, 2700.0, solidThermalConductivity2_, 0.0, 1.0, "lambda_eff_soil2");
            gnuplot.plot("lambdaeff");
        }
    }

private:
    GlobalPosition bBoxMin_;
    GlobalPosition bBoxMax_;
    Scalar lengthPM_;
    Scalar heightPM_;

    bool plotMaterialLaw_;
    bool randomField_;
    int numberOfSoilLayers_;
    int firstSoilLayerIdx_;
    int secondSoilLayerIdx_;
    bool verticalLayers_;

    Scalar permeability1_;
    Scalar porosity1_;
    Scalar alphaBJ1_;
    Scalar solidThermalConductivity1_;
    MaterialLawParams spatialParams1_;

    Scalar permeability2_;
    Scalar porosity2_;
    Scalar alphaBJ2_;
    Scalar solidThermalConductivity2_;
    MaterialLawParams spatialParams2_;
};

} // end namespace

#endif // DUMUX_WINDTUNNELSPATIALPARAMETERS_HH
