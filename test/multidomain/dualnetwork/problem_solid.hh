// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_TEST_MULTIDOMAIN_DUALNETWORK_PROBLEM_SOLID_HH
#define DUMUX_TEST_MULTIDOMAIN_DUALNETWORK_PROBLEM_SOLID_HH

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/porousmediumflow/problem.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/common/numeqvector.hh>

namespace Dumux {

/*!
 * \brief Heat problem with multiple solid spheres
 */
template <class TypeTag>
class SolidSubProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;

    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

    enum class SourceMode {conduction, convection, max};

public:
    template<class SpatialParams>
    SolidSubProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                    std::shared_ptr<SpatialParams> spatialParams,
                    std::shared_ptr<const CouplingManager> couplingManager)
    : ParentType(gridGeometry, spatialParams,  "Solid"), couplingManager_(couplingManager)
    {
        problemName_  =  getParam<std::string>("Vtk.OutputName") + "_" + getParamFromGroup<std::string>(this->paramGroup(), "Problem.Name");
        initialTemperature_ = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.InitialTemperature");
        temperatureBottom_ = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.BottomTemperature");
        temperatureIn_ = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.InletTemperature");
        enableCoupling_ = getParamFromGroup<bool>(this->paramGroup(), "Problem.EnableCoupling", true);
        heatingOn_ = getParamFromGroup<bool>(this->paramGroup(), "Problem.EnableHeating", true);
        useRobinInlet_ = getParamFromGroup<bool>(this->paramGroup(), "Problem.UseRobinInlet", true);

        inletIndex_ = getParamFromGroup<int>(this->paramGroup(), "Problem.InletIndex");
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

    BoundaryTypes boundaryTypes(const Element &element, const SubControlVolume& scv) const
    {
        BoundaryTypes values;
        values.setAllNeumann();

        if ((!useRobinInlet_ && onInletBoundary_(scv)) || (onHeaterBoundary_(scv) && heatingOn_))
            values.setAllDirichlet();

        return values;
    }

    template<class ElementVolumeVariables>
    NumEqVector source(const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolume& scv) const
    {
        NumEqVector value = 0.0;

        if (enableCoupling_ && couplingManager_->isCoupledPore(CouplingManager::solidDomainIdx, scv.dofIndex()))
        {
            if (sourceMode_ == SourceMode::conduction)
            {
                value[Indices::energyEqIdx] = couplingManager_->conductionSource(CouplingManager::solidDomainIdx,
                                                                                 element, fvGeometry, elemVolVars, scv);
            }
            else if (sourceMode_ == SourceMode::convection)
            {
                value[Indices::energyEqIdx] = couplingManager_->convectionSource(CouplingManager::solidDomainIdx,
                                                                                 element, fvGeometry, elemVolVars, scv);
            }
            else
            {
                const Scalar condSource = couplingManager_->conductionSource(CouplingManager::solidDomainIdx,
                                                                             element, fvGeometry, elemVolVars, scv);
                const Scalar convSource = couplingManager_->convectionSource(CouplingManager::solidDomainIdx,
                                                                             element, fvGeometry, elemVolVars, scv);

                using std::abs;
                if (abs(condSource) > abs(convSource))
                    value[Indices::energyEqIdx] = condSource;
                else
                    value[Indices::energyEqIdx] = convSource;
            }
        }

        value[Indices::energyEqIdx] += robinInletHeatFlux(element, fvGeometry, elemVolVars, scv);

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
            const Scalar r = this->spatialParams().poreExtendedRadius(scv.dofIndex());
            const auto& volVars = elemVolVars[scv];
            static const Scalar lambdaSolid = getParam<Scalar>("1.Component.SolidThermalConductivity");
            static const Scalar lambdaFluid = getParam<Scalar>("2.Component.LiquidThermalConductivity");
            static const Scalar dPadding = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.PaddingThickness");
            static const Scalar robinShapeFactor = getParamFromGroup<Scalar>(this->paramGroup(), "Problem.RobinShapeFactor", 1.0);

            const Scalar lambdaHarmonic = (r + dPadding)/(r/lambdaSolid + dPadding/lambdaFluid);
            const Scalar A = M_PI * r*r * robinShapeFactor;

            flux = (temperatureIn_ - volVars.temperature())
                   * A * lambdaHarmonic / (r + dPadding) / (scv.volume() * fvGeometry.gridGeometry().coordinationNumber()[scv.dofIndex()]);
        }

        return flux;
    }

    PrimaryVariables dirichlet(const Element& element, const SubControlVolume& scv) const
    {
        auto values = initialAtPos(scv.dofPosition());

        if (onInletBoundary_(scv))
        {
            values = temperatureIn_;

            if (!dirichletValuesForOutput_.empty())
                values = dirichletValuesForOutput_[scv.dofIndex()];
        }
        else if (onHeaterBoundary_(scv))
            values = temperatureBottom_;

        return values;
    }

    template<class ElementVolumeVariables, class ElementFluxVariablesCache>
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElementFluxVariablesCache& elemFluxVarsCache,
                        const SubControlVolumeFace& scvf) const
    {
        return NumEqVector(0.0);
    }

    PrimaryVariables initialAtPos(const GlobalPosition& pos) const
    {
        PrimaryVariables values(initialTemperature_);
        return values;
    }

    template<class Sol>
    void setInletToDirichletForOutput(const Sol& sol)
    {
        useRobinInlet_ = false;
        dirichletValuesForOutput_.resize(sol.size());

        for (const auto& vertex : vertices(this->gridGeometry().gridView()))
        {
            const auto vIdx = this->gridGeometry().vertexMapper().index(vertex);
            dirichletValuesForOutput_[vIdx] = sol[vIdx][Indices::temperatureIdx];
        }
    }

private:

    bool onInletBoundary_(const SubControlVolume& scv) const
    { return this->gridGeometry().poreLabel(scv.dofIndex()) == inletIndex_; }

    bool onHeaterBoundary_(const SubControlVolume& scv) const
    { return this->gridGeometry().poreLabel(scv.dofIndex()) == heaterIndex_; }

    std::shared_ptr<const CouplingManager> couplingManager_;

    std::string problemName_;
    Scalar initialTemperature_;
    Scalar temperatureIn_;
    Scalar temperatureBottom_;
    bool heatingOn_;
    bool enableCoupling_;
    int inletIndex_;
    int heaterIndex_;
    SourceMode sourceMode_;
    bool useRobinInlet_;
    std::vector<Scalar> dirichletValuesForOutput_;
};
} // end namespace Dumux

#endif
