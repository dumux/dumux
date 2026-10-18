// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup EmbeddedTests
 * \brief The vessel problem of the vessel tip benchmark.
 */
#ifndef DUMUX_TEST_EMBEDDED_1D3D_VESSEL_TIP_BLOODFLOW_PROBLEM_HH
#define DUMUX_TEST_EMBEDDED_1D3D_VESSEL_TIP_BLOODFLOW_PROBLEM_HH

#include <cmath>
#include <memory>
#include <vector>

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

#include <dumux/porousmediumflow/problem.hh>

#include "coupling.hh"
#include "exactsolution.hh"

namespace Dumux {

/*!
 * \ingroup EmbeddedTests
 * \brief The vessel problem of the vessel tip benchmark.
 *
 * The vessel pressure solves \f$ -\partial_s (\kappa \partial_s p_v) = -q \f$ on \f$ [0, L] \f$ with the
 * exact pressure at \f$ s = 0 \f$ and zero axial flux at the closed tip \f$ s = L \f$.
 */
template <class TypeTag>
class BloodFlowProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PointSource = GetPropType<TypeTag, Properties::PointSource>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename GridGeometry::GlobalCoordinate;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;
    using ExactSolution = VesselTipExactSolution<Scalar>;

    static constexpr auto distribution = vesselSourceDistribution<CouplingManager>();

public:
    BloodFlowProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                     std::shared_ptr<CouplingManager> couplingManager,
                     std::shared_ptr<const ExactSolution> exact)
    : ParentType(gridGeometry, "Vessel")
    , couplingManager_(couplingManager)
    , exactSolution_(exact)
    {
        name_ = getParam<std::string>("Vtk.OutputName") + "_" + getParamFromGroup<std::string>(this->paramGroup(), "Problem.Name");

        exactPressure_.resize(this->gridGeometry().numDofs());
        std::vector<bool> evaluated(this->gridGeometry().numDofs(), false);
        meanExchangeRate_.assign(this->gridGeometry().numDofs(), 0.0);
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
            values.setAllDirichlet();
        else
            values.setAllNeumann();
        return values;
    }

    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(exactSolution(globalPos)); }

    void addPointSources(std::vector<PointSource>& pointSources) const
    { pointSources = this->couplingManager().lowDimPointSources(); }

    template<class ElementVolumeVariables>
    void pointSource(PointSource& source,
                     const Element& element,
                     const FVElementGeometry& fvGeometry,
                     const ElementVolumeVariables& elemVolVars,
                     const SubControlVolume& scv) const
    {
        const Scalar pressure1D = this->couplingManager().lowDimPriVars(source.id())[Indices::pressureIdx];
        const Scalar pressure3D = this->couplingManager().bulkPriVars(source.id())[Indices::pressureIdx];

        Scalar sourceValue = couplingCoefficients_[source.id()]*(pressure3D - pressure1D);
        if constexpr (distribution == VesselSourceDistribution::kernel)
            sourceValue *= this->couplingManager().fluxScalingFactor(source.id());

        source = sourceValue*source.quadratureWeight()*source.integrationElement();
    }

    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(0.0); }

    //! The exact pressure at the degrees of freedom
    const std::vector<Scalar>& exactSolutionAtDofs() const
    { return exactPressure_; }

    //! The exact vessel pressure
    Scalar exactSolution(const GlobalPosition& globalPos) const
    { return exactSolution_->vesselPressure(globalPos[2]); }

    /*!
     * \brief The total exchange rate and the relative error of the exchange rate
     *
     * With \f$ Q_K = \int_K q \,\mathrm{d}s \f$ the exchange over the control volume \f$ K \f$, the error is
     * \f$ \left( \sum_K (Q_K - Q_{K,e})^2 / |K| \right)^{1/2} / \left( \sum_K Q_{K,e}^2 / |K| \right)^{1/2} \f$,
     * the relative discrete \f$ L^2 \f$ error of the mean exchange rate per control volume.
     * The mean exchange rate \f$ Q_K / |K| \f$ is added to the VTK output.
     */
    std::pair<Scalar, Scalar> computeSourceIntegral(const SolutionVector& sol, const GridVariables& gridVars)
    {
        std::vector<Scalar> exchange(this->gridGeometry().numDofs(), 0.0);
        std::vector<Scalar> length(this->gridGeometry().numDofs(), 0.0);
        auto fvGeometry = localView(this->gridGeometry());
        auto elemVolVars = localView(gridVars.curGridVolVars());
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            fvGeometry.bindElement(element);
            elemVolVars.bindElement(element, fvGeometry, sol);

            for (auto&& scv : scvs(fvGeometry))
            {
                // the vessel loses q per unit length
                const auto sources = this->scvPointSources(element, fvGeometry, elemVolVars, scv);
                exchange[scv.dofIndex()] -= sources[Indices::conti0EqIdx]*scv.volume()*elemVolVars[scv].extrusionFactor();
                length[scv.dofIndex()] += scv.volume();
            }
        }

        Scalar total = 0.0;
        Scalar errorSquared = 0.0;
        Scalar normSquared = 0.0;
        const Scalar q0 = exactSolution_->exchangeRate();
        for (std::size_t i = 0; i < exchange.size(); ++i)
        {
            const Scalar exactExchange = q0*length[i];
            total += exchange[i];
            meanExchangeRate_[i] = exchange[i]/length[i];
            errorSquared += (exchange[i] - exactExchange)*(exchange[i] - exactExchange)/length[i];
            normSquared += exactExchange*exactExchange/length[i];
        }

        using std::sqrt;
        return {total, sqrt(errorSquared/normSquared)};
    }

    template<class VtkOutputModule>
    void addVtkOutputFields(VtkOutputModule& vtk) const
    {
        vtk.addField(exactPressure_, "exact pressure");
        vtk.addField(meanExchangeRate_, "exchange rate");
    }

    const CouplingManager& couplingManager() const
    { return *couplingManager_; }

private:
    std::vector<Scalar> exactPressure_;
    std::vector<Scalar> meanExchangeRate_;
    std::vector<Scalar> couplingCoefficients_;

    static constexpr Scalar eps_ = 1.5e-7;
    std::string name_;

    std::shared_ptr<CouplingManager> couplingManager_;
    std::shared_ptr<const ExactSolution> exactSolution_;
};

} // end namespace Dumux

#endif
