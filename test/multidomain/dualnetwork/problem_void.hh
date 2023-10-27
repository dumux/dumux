// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief A test problem for the one-phase pore network model.
 */
#ifndef DUMUX_TEST_MULTIDOMAIN_DUALNETWORK_PROBLEM_VOID_HH
#define DUMUX_TEST_MULTIDOMAIN_DUALNETWORK_PROBLEM_VOID_HH

#include <dumux/porousmediumflow/problem.hh>
#include <dumux/porenetwork/1p/model.hh>

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/numeqvector.hh>

namespace Dumux {

template <class TypeTag>
class VoidSubProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using Labels = GetPropType<TypeTag, Properties::Labels>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;

    using Element = typename GridView::template Codim<0>::Entity;
    using Vertex = typename GridView::template Codim<GridView::dimension>::Entity;

    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

    enum class SourceMode {conduction, convection, max};

public:
    template<class SpatialParams>
    VoidSubProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                   std::shared_ptr<SpatialParams> spatialParams,
                   std::shared_ptr<const CouplingManager> couplingManager)
    : ParentType(gridGeometry, spatialParams, "Void"), couplingManager_(couplingManager)
    {
        problemName_ = getParam<std::string>("Vtk.OutputName") + "_" + getParamFromGroup<std::string>(this->paramGroup(), "Problem.Name");
        pressureIn_ = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.InletPressure");
        pressureOut_ = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.OutletPressure");
        temperatureInitial_ = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.InitialTemperature");
        pressureInitial_ = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.InitialPressure");
        temperatureIn_ = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.InletTemperature");
        temperatureBottom_ = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.BottomTemperature");
        enableCoupling_ = getParamFromGroup<bool>(this->paramGroup(), "Problem.EnableCoupling", true);
        heatingOn_ = getParamFromGroup<bool>(this->paramGroup(), "Problem.EnableHeating", true);
        fixedOutletTemperature_ = getParamFromGroup<bool>(this->paramGroup(), "Problem.FixedOutletTemperature", false);
        useRobinInlet_ = getParamFromGroup<bool>(this->paramGroup(), "Problem.UseRobinInlet", true);

        if (fixedOutletTemperature_)
            temperatureOut_ = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.OutletTemperature");

        inletIndex_ = getParamFromGroup<int>(this->paramGroup(), "Problem.InletIndex");
        outletIndex_ = getParamFromGroup<int>(this->paramGroup(), "Problem.OutletIndex");
        heaterIndex_ = getParamFromGroup<int>(this->paramGroup(), "Problem.HeaterIndex");

        const auto mode = getParamFromGroup<std::string>(this->paramGroup(), "Problem.DualNetworkSourceMode", "max");
        if (mode == "conduction")
            sourceMode_ = SourceMode::conduction;
        else if (mode == "convection")
            sourceMode_ = SourceMode::convection;
        else
            sourceMode_ = SourceMode::max;
    }

    const std::string& name() const
    { return problemName_; }

    void setGridVariables(std::shared_ptr<GridVariables> gridVars)
    {
        gridVars_ = gridVars;
    }

    BoundaryTypes boundaryTypes(const Element& element, const SubControlVolume& scv) const
    {
        BoundaryTypes bcTypes;
        if (onInletBoundary_(scv))
        {
            bcTypes.setDirichlet(Indices::pressureIdx);
            if (!useRobinInlet_)
                bcTypes.setDirichlet(Indices::temperatureIdx);
            else
                bcTypes.setNeumann(Indices::energyEqIdx);
        }
        else if (onOutletBoundary(scv))
        {
            // bcTypes.setAllNeumann();
            bcTypes.setDirichlet(Indices::pressureIdx);
            if (fixedOutletTemperature_)
                bcTypes.setDirichlet(Indices::temperatureIdx);
            else
                bcTypes.setNeumann(Indices::energyEqIdx);
        }
        else if (onHeaterBoundary_(scv) && heatingOn_)
        {
            bcTypes.setDirichlet(Indices::temperatureIdx);
            bcTypes.setNeumann(Indices::conti0EqIdx);
        }
        else // neuman for the remaining boundaries
            bcTypes.setAllNeumann();

        // treat insular pores connected only to the solid phase
        const auto poreLabel = this->gridGeometry().poreLabel(scv.dofIndex());
        if (poreLabel == 22)
            bcTypes.setDirichlet(Indices::pressureIdx);

        return bcTypes;
    }

    PrimaryVariables dirichlet(const Element& element,
                               const SubControlVolume& scv) const
    {
        PrimaryVariables values(0.0);
        if (onInletBoundary_(scv))
        {
            values[Indices::pressureIdx] = pressureIn_;
            values[Indices::temperatureIdx] = temperatureIn_;

            if (!dirichletValuesForOutput_.empty())
                values[Indices::temperatureIdx] = dirichletValuesForOutput_[scv.dofIndex()];
        }
        else if (onHeaterBoundary_(scv))
            values[Indices::temperatureIdx] = temperatureBottom_;
        else
        {
            values[Indices::pressureIdx] = pressureOut_;
            values[Indices::temperatureIdx] = temperatureOut_;

            if (!dirichletValuesForOutput_.empty())
                values[Indices::temperatureIdx] = dirichletValuesForOutput_[scv.dofIndex()];
        }

        return values;
    }


    template<class ElementVolumeVariables>
    NumEqVector source(const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolume& scv) const
    {
        NumEqVector value(0.0);

        if (enableCoupling_ && couplingManager_->isCoupledPore(CouplingManager::voidDomainIdx, scv.dofIndex()))
        {
            if (sourceMode_ == SourceMode::conduction)
            {
                value[Indices::energyEqIdx] = couplingManager_->conductionSource(
                    CouplingManager::voidDomainIdx, element, fvGeometry, elemVolVars, scv
                );
            }
            else if (sourceMode_ == SourceMode::convection)
            {
                value[Indices::energyEqIdx] = couplingManager_->convectionSource(
                    CouplingManager::voidDomainIdx, element, fvGeometry, elemVolVars, scv
                );
            }
            else
            {
                const Scalar condSource = couplingManager_->conductionSource(
                    CouplingManager::voidDomainIdx, element, fvGeometry, elemVolVars, scv
                );
                const Scalar convSource = couplingManager_->convectionSource(
                    CouplingManager::voidDomainIdx, element, fvGeometry, elemVolVars, scv
                );
                using std::abs;
                if (abs(condSource) > abs(convSource))
                    value[Indices::energyEqIdx] = condSource;
                else
                    value[Indices::energyEqIdx] = convSource;
            }
         }

        value[Indices::energyEqIdx] += robinInletHeatFlux(element, fvGeometry, elemVolVars, scv);

        // outflow condition for heat
        if (onOutletBoundary(scv) && !fixedOutletTemperature_)
            value[Indices::energyEqIdx] += heatOutFlowCondition(element, fvGeometry, elemVolVars, scv);

        return value;
    }


    template<class ElementVolumeVariables>
    Scalar heatOutFlowCondition(const Element& element,
                                const FVElementGeometry& fvGeometry,
                                const ElementVolumeVariables& elemVolVars,
                                const SubControlVolume& scv) const
    {
        Scalar value = 0.0;
        using FluxVariables = GetPropType<TypeTag, Properties::FluxVariables>;
        auto elemFluxVarsCache = localView(gridVars_->gridFluxVarsCache());
        elemFluxVarsCache.bindElement(element, fvGeometry, elemVolVars);

        for (auto&& scvf : scvfs(fvGeometry))
        {
            FluxVariables fluxVars;
            fluxVars.init(*this, element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);
            const Scalar flux = fluxVars.advectiveFlux(0,
                [&](const auto& volVars)
                {
                    return elemVolVars[scv].mobility(0)
                            * elemVolVars[scv].density(0)
                            * elemVolVars[scv].enthalpy(0);
                }
            );

            const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());

            if (insideScv.dofIndex() == scv.dofIndex())
                value += flux / scv.volume();
            else
                value -= flux / scv.volume();
        }
        return value;
    }

    template<class ElementVolumeVariables>
    Scalar robinInletHeatFlux(const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolume& scv) const
    {
        Scalar flux = 0.0;

        if (useRobinInlet_ && onInletBoundary_(scv))
        {
            flux += robinInletAdvectiveHeatFlux(element, fvGeometry, elemVolVars, scv);
            flux += robinInletConductiveHeatFlux(element, fvGeometry, elemVolVars, scv);
        }

        return flux;
    }

    template<class ElementVolumeVariables>
    Scalar robinInletAdvectiveHeatFlux(const Element& element,
                                    const FVElementGeometry& fvGeometry,
                                    const ElementVolumeVariables& elemVolVars,
                                    const SubControlVolume& scv) const
    {
        using FluxVariables = GetPropType<TypeTag, Properties::FluxVariables>;
        auto elemFluxVarsCache = localView(gridVars_->gridFluxVarsCache());
        elemFluxVarsCache.bindElement(element, fvGeometry, elemVolVars);
        const auto& volVars = elemVolVars[scv];

        for (auto&& scvf : scvfs(fvGeometry))
        {
            FluxVariables fluxVars;
            fluxVars.init(*this, element, fvGeometry, elemVolVars, scvf, elemFluxVarsCache);
            const Scalar enthalypy = ElementVolumeVariables::VolumeVariables::FluidSystem::enthalpy(
                temperatureIn_, volVars.pressure(0)
            );

            const Scalar result = fluxVars.advectiveFlux(0,
                [&](const auto& v)
                {
                    return volVars.mobility(0)
                            * volVars.density(0)*enthalypy;
                }
            );
            const auto& insideScv = fvGeometry.scv(scvf.insideScvIdx());

            if (insideScv.dofIndex() == scv.dofIndex())
                return result / scv.volume();
            else
                return  -result / scv.volume();
        }

        DUNE_THROW(Dune::InvalidStateException, "Flux failed");
    }

    template<class ElementVolumeVariables>
    Scalar robinInletConductiveHeatFlux(const Element& element,
                                        const FVElementGeometry& fvGeometry,
                                        const ElementVolumeVariables& elemVolVars,
                                        const SubControlVolume& scv) const
    {
        const auto& volVars = elemVolVars[scv];
        const Scalar r = volVars.poreInscribedRadius();
        static const Scalar lambdaFluid = getParam<Scalar>("2.Component.LiquidThermalConductivity");
        static const Scalar dPadding = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.PaddingThickness");
        static const Scalar robinShapeFactor = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.RobinShapeFactor", 1.0);
        const Scalar A = M_PI * r*r * robinShapeFactor;

        return (temperatureIn_ - volVars.temperature())
            * A * lambdaFluid / (r + dPadding) / (scv.volume() * fvGeometry.gridGeometry().coordinationNumber()[scv.dofIndex()]);
    }

    PrimaryVariables initial(const Vertex& vertex) const
    {
        PrimaryVariables values(0.0);
        values[Indices::pressureIdx] = pressureInitial_;
        values[Indices::temperatureIdx] = temperatureInitial_;
        return values;
    }

    int outletPoreLabel() const
    { return outletIndex_; }

    int inletPoreLabel() const
    { return inletIndex_; }

    int heaterPoreLabel() const
    { return heaterIndex_; }

    template<class Sol>
    void setOutletToDirichletForOutput(const Sol& sol)
    {
        fixedOutletTemperature_ = true;
        useRobinInlet_ = false;
        dirichletValuesForOutput_.resize(sol.size());

        for (const auto& vertex : vertices(this->gridGeometry().gridView()))
        {
            const auto vIdx = this->gridGeometry().vertexMapper().index(vertex);
            dirichletValuesForOutput_[vIdx] = sol[vIdx][Indices::temperatureIdx];
        }
    }

    bool onOutletBoundary(const SubControlVolume& scv) const
    { return this->gridGeometry().poreLabel(scv.dofIndex()) == outletIndex_; }

    auto sourceMode() const
    { return sourceMode_; }

private:

    bool onInletBoundary_(const SubControlVolume& scv) const
    { return this->gridGeometry().poreLabel(scv.dofIndex()) == inletIndex_; }

    bool onHeaterBoundary_(const SubControlVolume& scv) const
    { return this->gridGeometry().poreLabel(scv.dofIndex()) == heaterIndex_; }

    std::shared_ptr<const CouplingManager> couplingManager_;
    std::shared_ptr<GridVariables> gridVars_;

    std::string problemName_;
    Scalar pressureIn_;
    Scalar pressureOut_;
    Scalar temperatureInitial_;
    Scalar pressureInitial_;
    Scalar temperatureIn_;
    Scalar temperatureOut_;
    Scalar temperatureBottom_;
    bool enableCoupling_;
    bool heatingOn_;
    bool fixedOutletTemperature_;
    int inletIndex_;
    int outletIndex_;
    int heaterIndex_;
    SourceMode sourceMode_;
    std::vector<Scalar> dirichletValuesForOutput_;
    bool useRobinInlet_;
};

} // end namespace Dumux

#endif
