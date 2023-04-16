// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup EmbeddedTests
 * \brief The root problem for the benchmark cases C1.2a/b from Schnepf et al (2020) https://doi.org/10.3389/fpls.2020.00316
 */
#ifndef DUMUX_TEST_ROOT_SOIL_BENCHMARK_ROOT_PROBLEM_HH
#define DUMUX_TEST_ROOT_SOIL_BENCHMARK_ROOT_PROBLEM_HH

#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/porousmediumflow/problem.hh>
#include <dumux/discretization/evalsolution.hh>
#include <dumux/discretization/elementsolution.hh>

#include "couplingreconstruction.hh"

namespace Dumux {

/*!
 * \ingroup EmbeddedTests
 * \brief The root problem for the benchmark cases C1.2a/b from Schnepf et al (2020) https://doi.org/10.3389/fpls.2020.00316
 */
template <class TypeTag>
class RootProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PointSource = GetPropType<TypeTag, Properties::PointSource>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<PrimaryVariables::size()>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using GlobalPosition = typename GridGeometry::GlobalCoordinate;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using Element = typename GridView::template Codim<0>::Entity;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

public:
    enum BCType
    { constCollarPressure, constTranspiration, cyclicTranspiration };

    template<class SpatialParams>
    RootProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                std::shared_ptr<SpatialParams> spatialParams,
                std::shared_ptr<CouplingManager> couplingManager,
                const GlobalPosition& domainSize)
    : ParentType(gridGeometry, spatialParams, "Root")
    , couplingManager_(couplingManager)
    {
        // read parameters from input file
        name_ = getParamFromGroup<std::string>(this->paramGroup(), "Problem.Name");
        criticalCollarPressure_ = getParam<Scalar>("BoundaryConditions.CriticalCollarPressure", -1.4e6); // in Pa
        initialRootPressure_ = getParam<Scalar>("BoundaryConditions.InitialRootPressure", -1e5); // in Pa
        dailyTranspirationRate_ = getParam<Scalar>("BoundaryConditions.DailyTranspirationRate", 1.0); // mm/day
        dailyTranspirationRate_ *= domainSize[0]*domainSize[1]/86400.0; // kg/s (domain size in m, kg->g=mm->m)

        const auto bcType = getParam<std::string>("BoundaryConditions.BoundaryType", "CyclicTranspiration");
        if (bcType == "CyclicTranspiration") bcType_ = BCType::cyclicTranspiration;
        else if (bcType == "ConstCollarPressure") bcType_ = BCType::constCollarPressure;
        else if (bcType == "ConstTranspiration") bcType_ = BCType::constTranspiration;
        else DUNE_THROW(Dune::InvalidStateException, "Unknown BCType: " << bcType);

        bcSmoothingParam_ = getParam<Scalar>("BoundaryConditions.SmoothingK", 0.1);
        kernelWidthFactor_ = getParam<Scalar>("MixedDimension.KernelWidthFactor");
    }

    const std::string& name() const
    { return name_; }

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;
        values.setAllNeumann();

        if (bcType_ == BCType::constCollarPressure)
            if (globalPos[2] + eps_ > this->gridGeometry().bBoxMax()[2])
                values.setAllDirichlet();

        return values;
    }

    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(criticalCollarPressure_); }

    template<class ElementVolumeVariables, class ElementFluxVarsCache>
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElementFluxVarsCache& elemFluxVarsCache,
                        const SubControlVolumeFace& scvf) const
    {
        NumEqVector values(0.0);

        // set transpiration rate at the root collar
        const auto globalPos = scvf.center();
        if (globalPos[2] + eps_ > this->gridGeometry().bBoxMax()[2])
        {
            // if the plant has water stress reduce the transpiration rate (imposing Dirichlet boundary condition weakly)
            const auto& volVars = elemVolVars[scvf.insideScvIdx()];
            const auto dist = fvGeometry.scv(scvf.insideScvIdx()).volume();
            const auto elemSol = elementSolution(element, elemVolVars, fvGeometry);
            const auto p = evalSolution(element, element.geometry(), elemSol, globalPos)[0];
            const auto eIdx = this->gridGeometry().elementMapper().index(element);
            const Scalar Kx = this->spatialParams().Kx(eIdx);
            const auto rho = volVars.density(0);
            const Scalar criticalTranspiration = rho*Kx*(p - criticalCollarPressure_ - rho*9.81*dist)/dist;
            values[Indices::conti0EqIdx] = smoothMin(potentialTranspirationRate(), criticalTranspiration, dailyTranspirationRate_*bcSmoothingParam_);
            values /= volVars.extrusionFactor() * scvf.area(); // convert from kg/s to kg/(s*m^2)
        }

        return values;
    }

    //! Smooth minimum function
    Scalar smoothMin(const Scalar a, const Scalar b, const Scalar k) const
    {
        using std::max; using std::min; using std::abs;
        const auto h = max(k-abs(a-b), 0.0 )/k;
        return min(a, b) - h*h*h*k*(1.0/6.0);
    }

    //! Compute potential transpiration rate in kg/s
    Scalar potentialTranspirationRate() const
    {
        // possibly make the transpiration rate dependent on the current root length for growth
        if (bcType_ == BCType::cyclicTranspiration)
            return dailyTranspirationRate_*std::sin(time_*2.0*M_PI / 86400.0 - M_PI/2.0) + dailyTranspirationRate_;
        else if (bcType_ == BCType::constTranspiration)
            return dailyTranspirationRate_;
        else
            return std::numeric_limits<Scalar>::max();
    }

    void addPointSources(std::vector<PointSource>& pointSources) const
    { pointSources = this->couplingManager().lowDimPointSources(); }

    template<class ElementVolumeVariables>
    void pointSource(PointSource& source,
                     const Element &element,
                     const FVElementGeometry& fvGeometry,
                     const ElementVolumeVariables& elemVolVars,
                     const SubControlVolume &scv) const
    {
        // compute source at every integration point
        const auto id = source.id();

        const Scalar p0 = this->couplingManager().bulkPriVars(id)[Indices::pressureIdx];
        const Scalar pRoot = this->couplingManager().lowDimPriVars(id)[Indices::pressureIdx];
        const Scalar radius = this->couplingManager().radius(id);
        const auto lEIdx = this->couplingManager().pointSourceData(id).lowDimElementIdx();
        const auto bEIdx = this->couplingManager().pointSourceData(id).bulkElementIdx();
        const Scalar Kr = this->spatialParams().Kr(lEIdx);
        const Scalar kernelRadius = kernelWidthFactor_*radius;
        const Scalar density = 1000;
        const Scalar viscosity = 1e-3;
        const Scalar delta = this->couplingManager().averageDistance(id);

        // sink defined as radial flow Jr * density [m^2 s-1]* [kg m-3]
        const auto sourceValue = -1.0*this->couplingManager().reconstruction().computeSource(
            bEIdx, lEIdx, p0, pRoot, kernelRadius, radius, Kr, density, viscosity, id, delta
        );
        source = sourceValue*source.quadratureWeight()*source.integrationElement();
    }

    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        if (bcType_ == BCType::constCollarPressure)
            return PrimaryVariables(criticalCollarPressure_);
        else
            return PrimaryVariables(initialRootPressure_);
    }

    //! Called after every time step
    //! Output the total global exchange term
    //! Output the individual source terms into sourceTerms
    void computeSourceIntegral(const SolutionVector& sol, const GridVariables& gridVars,
                               std::vector<Scalar>& sourceTerms)
    {
        PrimaryVariables source(0.0);
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            auto fvGeometry = localView(this->gridGeometry());
            fvGeometry.bindElement(element);

            auto elemVolVars = localView(gridVars.curGridVolVars());
            elemVolVars.bindElement(element, fvGeometry, sol);

            PrimaryVariables localSource(0.0);
            for (const auto& scv : scvs(fvGeometry))
            {
                auto pointSources = this->scvPointSources(element, fvGeometry, elemVolVars, scv);
                pointSources *= scv.volume()*elemVolVars[scv].extrusionFactor();
                localSource += pointSources;
            }

            const auto eIdx = this->gridGeometry().elementMapper().index(element);
            sourceTerms[eIdx] = localSource[0];
            source += localSource;
        }

        std::cout << "Global integrated source (root): " << source << " (kg/s) / "
                  <<                           source*3600*24*1000 << " (g/day)" << '\n';
    }

    //! compute the actual transpiration rate
    Scalar computeActualTranspirationRate(const SolutionVector& sol, const GridVariables& gridVars, bool verbose = true) const
    {
        NumEqVector transpirationRate(0.0);
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            auto fvGeometry = localView(this->gridGeometry());
            fvGeometry.bindElement(element);

            auto elemVolVars = localView(gridVars.curGridVolVars());
            elemVolVars.bindElement(element, fvGeometry, sol);

            for (const auto& scvf : scvfs(fvGeometry))
                if (scvf.boundary())
                    transpirationRate += this->neumann(element, fvGeometry, elemVolVars, 0.0, scvf)
                                         *scvf.area()*elemVolVars[scvf.insideScvIdx()].extrusionFactor();
        }

        if (verbose)
        {
            std::cout << "Actual transpiration rate:       " << transpirationRate << " (kg/s) / "
                      << transpirationRate[0]*86400*1000 << " (g/day)";
            if (bcType_ != BCType::constCollarPressure)
                std::cout << " / Potential transpiration rate:    " << potentialTranspirationRate() << " (kg/s) / "
                          << potentialTranspirationRate()*86400*1000 << " (g/day) / ";
            std::cout << std::endl;
        }

        return transpirationRate[0];
    }

    //! Get the coupling manager
    const CouplingManager& couplingManager() const
    { return *couplingManager_; }

    //! set the current time for evaluation of time-dependent boundary conditions
    void setTime(const Scalar t)
    { time_= t; }

    //! get the boundary condition type
    BCType bcType() const
    { return bcType_; }

private:
    Scalar initialRootPressure_, criticalCollarPressure_, dailyTranspirationRate_;
    Scalar time_;
    Scalar kernelWidthFactor_;

    BCType bcType_;
    Scalar bcSmoothingParam_;

    static constexpr Scalar eps_ = 1.5e-7;
    std::string name_;

    std::shared_ptr<CouplingManager> couplingManager_;
};

} // end namespace Dumux

#endif
