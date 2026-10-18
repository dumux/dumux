// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup EmbeddedTests
 * \brief The tissue problem of the vessel tip benchmark.
 */
#ifndef DUMUX_TEST_EMBEDDED_1D3D_VESSEL_TIP_TISSUE_PROBLEM_HH
#define DUMUX_TEST_EMBEDDED_1D3D_VESSEL_TIP_TISSUE_PROBLEM_HH

#include <cassert>
#include <cmath>
#include <memory>
#include <vector>

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/numeqvector.hh>

#include <dumux/porousmediumflow/problem.hh>

#include "coupling.hh"
#include "exactsolution.hh"

namespace Dumux {

/*!
 * \ingroup EmbeddedTests
 * \brief The tissue problem of the vessel tip benchmark.
 *
 * The tissue pressure solves \f$ -\Delta p_t = q \delta \f$ in \f$ [-1, 1]^2 \times [0, 1] \f$. The plane
 * \f$ z = 0 \f$ is a symmetry plane (no flow), the exact solution is prescribed on the other boundaries.
 */
template <class TypeTag>
class TissueProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using GridView = typename GridGeometry::GridView;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using PointSource = GetPropType<TypeTag, Properties::PointSource>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename GridGeometry::GlobalCoordinate;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;
    using ExactSolution = VesselTipExactSolution<Scalar>;

    static constexpr auto distribution = vesselSourceDistribution<CouplingManager>();

public:
    TissueProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                  std::shared_ptr<CouplingManager> couplingManager,
                  std::shared_ptr<const ExactSolution> exact)
    : ParentType(gridGeometry, "Tissue")
    , couplingManager_(couplingManager)
    , exactSolution_(exact)
    {
        name_ = getParam<std::string>("Vtk.OutputName") + "_" + getParamFromGroup<std::string>(this->paramGroup(), "Problem.Name");

        // several sub-control volumes share a dof (box), the semi-analytic solution is evaluated once per dof
        exactPressure_.resize(this->gridGeometry().numDofs());
        std::vector<bool> evaluated(this->gridGeometry().numDofs(), false);
        auto fvGeometry = localView(this->gridGeometry());
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            fvGeometry.bindElement(element);
            for (auto&& scv : scvs(fvGeometry))
            {
                if (!evaluated[scv.dofIndex()])
                {
                    exactPressure_[scv.dofIndex()] = exactSolution(scv.dofPosition());
                    evaluated[scv.dofIndex()] = true;
                }
            }
        }
    }

    const std::string& name() const
    { return name_; }

    //! Evaluate the coupling coefficient for each coupling source (call after the coupling manager is initialized)
    void computeCouplingCoefficients()
    { couplingCoefficients_ = vesselTipCouplingCoefficients(*couplingManager_, *exactSolution_); }

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        if (globalPos[2] < this->gridGeometry().bBoxMin()[2] + eps_)
            values.setAllNeumann();
        else
            values.setAllDirichlet();
        return values;
    }

    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(exactSolution(globalPos)); }

    //! The coupling sources for the line and surface source methods
    void addPointSources(std::vector<PointSource>& pointSources) const
    { pointSources = this->couplingManager().bulkPointSources(); }

    template<class ElementVolumeVariables>
    void pointSource(PointSource& source,
                     const Element& element,
                     const FVElementGeometry& fvGeometry,
                     const ElementVolumeVariables& elemVolVars,
                     const SubControlVolume& scv) const
    {
        if constexpr (distribution != VesselSourceDistribution::kernel)
        {
            const Scalar pressure3D = this->couplingManager().bulkPriVars(source.id())[Indices::pressureIdx];
            const Scalar pressure1D = this->couplingManager().lowDimPriVars(source.id())[Indices::pressureIdx];
            const Scalar sourceValue = couplingCoefficients_[source.id()]*(pressure1D - pressure3D);
            source = sourceValue*source.quadratureWeight()*source.integrationElement();
        }
    }

    //! The coupling sources for the kernel method
    template<class ElementVolumeVariables>
    NumEqVector source(const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolume& scv) const
    {
        NumEqVector source(0.0);

        if constexpr (distribution == VesselSourceDistribution::kernel)
        {
            const auto eIdx = this->gridGeometry().elementMapper().index(element);
            const auto& sourceIds = this->couplingManager().bulkSourceIds(eIdx, scv.indexInElement());
            const auto& sourceWeights = this->couplingManager().bulkSourceWeights(eIdx, scv.indexInElement());

            for (int i = 0; i < sourceIds.size(); ++i)
            {
                const auto id = sourceIds[i];
                const Scalar xi = this->couplingManager().fluxScalingFactor(id);
                const Scalar pressure3D = this->couplingManager().bulkPriVars(id)[Indices::pressureIdx];
                const Scalar pressure1D = this->couplingManager().lowDimPriVars(id)[Indices::pressureIdx];
                source[Indices::conti0EqIdx] += couplingCoefficients_[id]*(pressure1D - pressure3D)*xi*sourceWeights[i];
            }

            source[Indices::conti0EqIdx] /= scv.volume()*elemVolVars[scv].extrusionFactor();
        }

        return source;
    }

    //! The flux scaling factor of the kernel method with the coupling coefficient at the coupling source
    Scalar fluxScalingFactor(const Scalar r, const Scalar R, const Scalar rho, std::size_t id) const
    {
        const auto& source = this->couplingManager().lowDimPointSources()[id];
        assert(source.id() == id);
        const Scalar beta = exactSolution_->couplingCoefficient(distribution, source.position()[2]);
        return kernelFluxScalingFactor(r, R, rho, beta);
    }

    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(0.0); }

    //! The exact pressure at the degrees of freedom
    const std::vector<Scalar>& exactSolutionAtDofs() const
    { return exactPressure_; }

    //! The exact tissue pressure
    Scalar exactSolution(const GlobalPosition& globalPos) const
    {
        using std::hypot;
        return exactSolution_->tissuePressure(distribution, hypot(globalPos[0], globalPos[1]), globalPos[2]);
    }

    //! The total exchange rate added to the tissue
    Scalar computeSourceIntegral(const SolutionVector& sol, const GridVariables& gridVars) const
    {
        Scalar source = 0.0;
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            const auto fvGeometry = localView(this->gridGeometry()).bindElement(element);
            const auto elemVolVars = localView(gridVars.curGridVolVars()).bindElement(element, fvGeometry, sol);

            for (auto&& scv : scvs(fvGeometry))
            {
                auto sources = this->scvPointSources(element, fvGeometry, elemVolVars, scv);
                sources += this->source(element, fvGeometry, elemVolVars, scv);
                source += sources[Indices::conti0EqIdx]*scv.volume()*elemVolVars[scv].extrusionFactor();
            }
        }
        return source;
    }

    template<class VtkOutputModule>
    void addVtkOutputFields(VtkOutputModule& vtk) const
    { vtk.addField(exactPressure_, "exact pressure"); }

    const CouplingManager& couplingManager() const
    { return *couplingManager_; }

private:
    std::vector<Scalar> exactPressure_;
    std::vector<Scalar> couplingCoefficients_;

    static constexpr Scalar eps_ = 1.5e-7;
    std::string name_;

    std::shared_ptr<CouplingManager> couplingManager_;
    std::shared_ptr<const ExactSolution> exactSolution_;
};

} // end namespace Dumux

#endif
