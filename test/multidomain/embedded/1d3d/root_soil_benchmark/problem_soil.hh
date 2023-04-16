// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup EmbeddedTests
 * \brief The soil problem for the benchmark cases C1.2a/b from Schnepf et al (2020) https://doi.org/10.3389/fpls.2020.00316
 */
#ifndef DUMUX_TEST_ROOT_SOIL_BENCHMARK_SOIL_PROBLEM_HH
#define DUMUX_TEST_ROOT_SOIL_BENCHMARK_SOIL_PROBLEM_HH

#include <dumux/common/math.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/common/numeqvector.hh>

#include <dumux/porousmediumflow/problem.hh>

namespace Dumux {

/*!
 * \ingroup EmbeddedTests
 * \brief The soil problem for the benchmark cases C1.2a/b from Schnepf et al (2020) https://doi.org/10.3389/fpls.2020.00316
 */
template <class TypeTag>
class SoilProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using GlobalPosition = typename GridGeometry::GlobalCoordinate;
    using Element = typename GridView::template Codim<0>::Entity;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<PrimaryVariables::size()>;
    using PointSource = GetPropType<TypeTag, Properties::PointSource>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

public:
    SoilProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry, "Soil")
    , couplingManager_(couplingManager)
    {
        // read parameters from input file
        name_ = getParamFromGroup<std::string>(this->paramGroup(), "Problem.Name");

        if (hasParam("BoundaryConditions.InitialSoilSaturationTop"))
        {
            const auto sw = getParam<Scalar>("BoundaryConditions.InitialSoilSaturationTop");
            pcTop_ = this->spatialParams().fluidMatrixInteractionAtPos(this->gridGeometry().bBoxMax()).pc(sw);
        }
        else
        {
            const auto pw = getParam<Scalar>("BoundaryConditions.InitialSoilPressureTop");
            pcTop_ = nonwettingReferencePressure() - pw;
        }

        kernelWidthFactor_ = getParam<Scalar>("MixedDimension.KernelWidthFactor");
    }

    const std::string& name() const
    { return name_; }

    Scalar nonwettingReferencePressure() const
    { return 1.0e5; }

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;
        values.setAllNeumann();
        return values;
    }

    template<class ElementVolumeVariables>
    NumEqVector source(const Element &element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolume &scv) const
    {
        NumEqVector values(0.0);

        const auto eIdx = this->gridGeometry().elementMapper().index(element);
        const auto& sourceIds = this->couplingManager().bulkSourceIds(eIdx);
        const auto& sourceWeights = this->couplingManager().bulkSourceWeights(eIdx);

        for (int i = 0; i < sourceIds.size(); ++i)
        {
            const auto id = sourceIds[i];
            const auto weight = sourceWeights[i];

            // compute source at every integration point
            const Scalar p0 = this->couplingManager().bulkPriVars(id)[Indices::pressureIdx];
            const Scalar pRoot = this->couplingManager().lowDimPriVars(id)[Indices::pressureIdx];
            const Scalar radius = this->couplingManager().radius(id);
            const Scalar kernelRadius = kernelWidthFactor_*radius;
            const Scalar Kr = this->couplingManager().Kr(id);
            const Scalar density = 1000;
            const Scalar viscosity = 1e-3;
            const Scalar delta = this->couplingManager().averageDistance(id);

            // sink defined as radial flow Jr * density [m^2 s-1]* [kg m-3]
            const auto lEIdx = this->couplingManager().pointSourceData(id).lowDimElementIdx();
            const auto bEIdx = this->couplingManager().pointSourceData(id).bulkElementIdx();
            const auto sourceValue = this->couplingManager().reconstruction().computeSource(
                bEIdx, lEIdx, p0, pRoot, kernelRadius, radius, Kr, density, viscosity, id, delta
            );
            values[Indices::conti0EqIdx] += sourceValue*weight;
        }

        const auto volume = scv.volume()*elemVolVars[scv].extrusionFactor();
        values[Indices::conti0EqIdx] /= volume;

        return values;
    }

    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables priVars(0.0);
        static const bool initialGravityProfile = getParam<bool>("BoundaryConditions.InitialGravityProfile");
        if (initialGravityProfile)
            priVars[0] = nonwettingReferencePressure() - pcTop_ - 9.81*1000*(globalPos[2] - this->gridGeometry().bBoxMax()[2]);
        else
            priVars[0] = nonwettingReferencePressure() - pcTop_;
        return priVars;
    }

    //! Called after every time step
    //! Output the total global exchange term
    void computeSourceIntegral(const SolutionVector& sol, const GridVariables& gridVars)
    {
        PrimaryVariables source(0.0);
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            auto fvGeometry = localView(this->gridGeometry());
            fvGeometry.bindElement(element);

            auto elemVolVars = localView(gridVars.curGridVolVars());
            elemVolVars.bindElement(element, fvGeometry, sol);

            for (auto&& scv : scvs(fvGeometry))
            {
                auto pointSources = this->source(element, fvGeometry, elemVolVars, scv);
                pointSources *= scv.volume()*elemVolVars[scv].extrusionFactor();
                source += pointSources;
            }
        }

        std::cout << "Global integrated source (soil): " << source << " (kg/s) / "
                  <<                           source*3600*24*1000 << " (g/day)" << '\n';

    }

    //! Get the coupling manager
    const CouplingManager& couplingManager() const
    { return *couplingManager_; }

private:
    Scalar pcTop_;
    Scalar kernelWidthFactor_;

    static constexpr Scalar eps_ = 1.5e-7;
    std::string name_;

    std::shared_ptr<CouplingManager> couplingManager_;
};

} // end namespace Dumux

#endif
