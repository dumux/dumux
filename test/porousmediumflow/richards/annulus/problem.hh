// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//

#ifndef DUMUX_RICHARDS_ANNULUS_PROBLEM_HH
#define DUMUX_RICHARDS_ANNULUS_PROBLEM_HH

#include <iostream>
#include <cmath>

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/numeqvector.hh>

#include <dumux/common/integrate.hh>
#include <dumux/nonlinear/findscalarroot.hh>
#include <dumux/discretization/extrusion.hh>

#include <dumux/porousmediumflow/problem.hh>

namespace Dumux {

template<class TypeTag, class FluidSystem>
class RichardsAnnulusProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;

    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<PrimaryVariables::size()>;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    RichardsAnnulusProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        // The default parameters correspond to the benchmark cases C1.1
        // of Schnepf et al (2020) https://doi.org/10.3389/fpls.2020.00316

        // material parameters
        k_ = getParam<Scalar>(this->spatialParams().paramGroup() + ".Permeability");
        static_assert(!FluidSystem::isCompressible(0), "Assumes incompressible fluid");
        rho_ = FluidSystem::H2O::liquidDensity(0, 0);
        static_assert(FluidSystem::viscosityIsConstant(0), "Assumes constant viscosity fluid");
        mu_ = FluidSystem::H2O::liquidViscosity(0, 0);

        // geometry
        innerRadius_ = gridGeometry->bBoxMin()[0];
        outerRadius_ = gridGeometry->bBoxMax()[0];

        // boundary conditions (volumetric flow rate in m^2/s, pressure/psi in Pa)
        const auto innerFlux = getParam<Scalar>("Problem.InnerFlux", 0.1); // positive means outflow (cm/day)
        innerFlowRate_ = innerFlux*2*M_PI*innerRadius_*0.01/(24*60*60);
        const auto outerFlux = getParam<Scalar>("Problem.OuterFlux", 0.0); // positive means outflow (cm/day)
        outerFlowRate_ = outerFlux*2*M_PI*outerRadius_*0.01/(24*60*60);

        const auto innerPressureHead = getParam<Scalar>("Problem.InnerPressureHeadInCm", -15000);
        const auto innerPressure = headToPressure(innerPressureHead);
        innerPsi_ = applyKirchhoffTransformation_(innerPressure);
        source_ = (innerFlowRate_ + outerFlowRate_)/(M_PI*(innerRadius_*innerRadius_ - outerRadius_*outerRadius_));

        // boundary condition types
        enableOuterNeumann_ = getParam<bool>("Problem.EnableOuterNeumann", true);
        enableInnerNeumann_ = getParam<bool>("Problem.EnableInnerNeumann", false);

        std::cout << "Inner flow rate: " << innerFlowRate_ << " m^2/s" << std::endl;
        std::cout << "Outer flow rate: " << outerFlowRate_ << " m^2/s" << std::endl;
        std::cout << "Critical inner pressure: " << innerPressure << " Pa "
                  << "(-> head: " << innerPressureHead << " cm)" << std::endl;

        // initial condition
        const auto initialHead = getParam<Scalar>("Problem.InitialPressureHeadInCm", -100);
        initialPressure_ = headToPressure(initialHead);
    }

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        if (enableOuterNeumann_ && globalPos[0] > this->gridGeometry().bBoxMax()[0] - 1e-10)
            values.setAllNeumann();
        else if (enableInnerNeumann_ && globalPos[0] < this->gridGeometry().bBoxMin()[0] + 1e-10)
            values.setAllNeumann();
        else
            values.setAllDirichlet();
        return values;
    }

    // negative means extraction
    NumEqVector sourceAtPos(const GlobalPosition& globalPos) const
    {
        if (enableSource_)
            return { -rho_*source_ };
        else
            return { 0.0 };
    }

    // positive means outflow
    NumEqVector neumannAtPos(const GlobalPosition& globalPos) const
    {
        if (globalPos[0] > this->gridGeometry().bBoxMax()[0] - 1e-10)
            return { rho_*outerFlowRate_/(2*M_PI*outerRadius_) };
        else
            return { rho_*innerFlowRate_/(2*M_PI*innerRadius_) };
    }

    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    { return exactSolution(globalPos); }

    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    {
        if (useStationaryAnalyticInitialSolution_)
            return exactSolution(globalPos);
        else
        {
            return { initialPressure_ };
        }
    }

    /*! \brief Analytical solution

        Analytical solution of stationary incompressible
        Richards equations with source term B:

          ψ = -µ/k A/(2π) ln(r) + r^2/4 µ/k B + C

        Three constants to fix with three boundary conditions:
         - flux at inner boundary (r=r_i): q_i = 2πr_i k/µ ∂ψ/∂r
         - flux at outer boundary (r=r_o): q_o = -2πr_o k/µ ∂ψ/∂r
         - ψ = T(p) at inner boundary -> ψ = ψ_ri

        Yields the solution in terms of q_i, q_o, ψ_ri:

          ψ = ψ_ri - µ/k A/(2π) ln(r/r_i) + (r^2 - r_i^2)/4 µ/k B
        where
          A = - q_i + πr_i^2 B
          B = (q_i + q_o)/(π(r_i^2 - r_o^2)) [source]
    */
    PrimaryVariables exactSolution(const GlobalPosition& globalPos, const Scalar innerPsi) const
    {
        const auto r = globalPos[0];
        const auto ri2 = innerRadius_*innerRadius_;
        const auto A = -innerFlowRate_ + M_PI*ri2*source_;
        const auto psi = innerPsi
            - A*mu_/(2*M_PI*k_)*std::log(r/innerRadius_)
            + source_*0.25*(r*r - ri2)*mu_/k_;

        return { applyInverseKirchhoffTransformation_(psi) };
    }

    PrimaryVariables exactSolution(const GlobalPosition& globalPos) const
    { return exactSolution(globalPos, innerPsi_); }

    Scalar nonwettingReferencePressure() const
    { return 1.0e5; };

    template<class SolutionVector, class GridVariables>
    Scalar computeTime(const SolutionVector& p, const GridVariables& vars, const Scalar initialSaturation) const
    {
        const auto& gg = this->gridGeometry();
        auto fvGeometry = localView(gg);
        auto elemVolVars = localView(vars.curGridVolVars());
        Scalar totalWaterVolume = 0.0;
        Scalar refWaterVolume = 0.0;
        using Extrusion = Extrusion_t<GridGeometry>;
        for (const auto& element : elements(gg.gridView()))
        {
            fvGeometry.bindElement(element);
            elemVolVars.bindElement(element, fvGeometry, p);
            for (const auto& scv : scvs(fvGeometry))
            {
                const auto& volVars = elemVolVars[scv];
                totalWaterVolume += Extrusion::volume(fvGeometry, scv)*volVars.porosity()*volVars.saturation(0);
                refWaterVolume += Extrusion::volume(fvGeometry, scv)*volVars.porosity()*initialSaturation;
            }
        }

        return (refWaterVolume - totalWaterVolume)/(innerFlowRate_ + outerFlowRate_);
    }

    template<class SolutionVector, class GridVariables>
    void updateStorageDerivative(const SolutionVector& p, const SolutionVector& pOld,
                                 const GridVariables& vars, const Scalar dt,
                                 std::vector<Scalar>& storageDerivative) const
    {
        const auto& gg = this->gridGeometry();
        auto fvGeometry = localView(gg);
        auto elemVolVars = localView(vars.curGridVolVars());
        auto oldElemVolVars = localView(vars.prevGridVolVars());
        storageDerivative.assign(storageDerivative.size(), 0.0);
        auto volumes = storageDerivative;
        using Extrusion = Extrusion_t<GridGeometry>;
        for (const auto& element : elements(gg.gridView()))
        {
            fvGeometry.bindElement(element);
            elemVolVars.bindElement(element, fvGeometry, p);
            oldElemVolVars.bindElement(element, fvGeometry, pOld);
            for (const auto& scv : scvs(fvGeometry))
            {
                const auto& volVars = elemVolVars[scv];
                const auto& oldVolVars = oldElemVolVars[scv];
                storageDerivative[scv.dofIndex()] += Extrusion::volume(fvGeometry, scv)*(volVars.saturation(0) - oldVolVars.saturation(0))/dt;
                volumes[scv.dofIndex()] += Extrusion::volume(fvGeometry, scv);
            }
        }

        for (int i = 0; i < storageDerivative.size(); ++i)
            storageDerivative[i] /= volumes[i];
    }

    // return saturation provided pressure in Pa
    const Scalar saturation(Scalar const pw) const
    {
        const auto& pcKrSwCurve = this->spatialParams().pcKrSwCurve();
        const Scalar pc = this->nonwettingReferencePressure() - pw;
        return pcKrSwCurve.sw(pc);
    }

    /*!
     * \brief convert pressure head in cm to pressure in Pa
     */
    Scalar headToPressure(const Scalar head) const
    {
        return this->nonwettingReferencePressure() + head / 100.0 * rho_ * 9.81;
    }

    /*!
     * \brief convert pressure in Pa to potential psi in Pa (Kirchhoff transformation)
     */
    Scalar pressureToPsi(const Scalar pw) const
    {
        return applyKirchhoffTransformation_(pw);
    }

    /*!
     * \brief Disable the source term (for the instationary problem)
     */
    void disableSource()
    { enableSource_ = false; }

    /*!
     * \brief Disable using the stationary initial solution (for the instationary problem)
     */
    void disableStationaryInitialSolution()
    { useStationaryAnalyticInitialSolution_ = false; }

    /*!
     * \brief Enable using Neumann boundary conditions on the inner boundary
     */
    void enableInnerNeumannBC()
    { enableInnerNeumann_ = true; }

private:

    /*!
     * \brief relative wetting phase permeability
     */
    Scalar kr_(Scalar pw) const
    {
        const auto& pcKrSwCurve = this->spatialParams().pcKrSwCurve();
        const Scalar pc = this->nonwettingReferencePressure() - pw;
        const Scalar sw = pcKrSwCurve.sw(pc);
        return pcKrSwCurve.krw(sw);
    }

    /*!
     * \brief transformed variable p->psi
     */
    Scalar applyKirchhoffTransformation_(Scalar pw) const
    {
        return integrateScalarFunction(
            [&](const double p){ return kr_(p); }, -1.5e6, pw, 1e-20
        );
    }

    /*!
     * \brief inverse transformation psi->p
     */
    Scalar applyInverseKirchhoffTransformation_(Scalar psi) const
    {
        const auto residual = [&, psi](const auto& p){
            return psi - applyKirchhoffTransformation_(p);
        };

        return findScalarRootBrent(-1.5e6, 1.0e6, residual, 1e-14, 200000);
    }

    Scalar k_, mu_, rho_;
    Scalar innerRadius_, outerRadius_, source_;
    Scalar innerFlowRate_, outerFlowRate_, innerPsi_;
    Scalar initialPressure_;
    bool enableOuterNeumann_, enableInnerNeumann_;
    bool useStationaryAnalyticInitialSolution_ = true;
    bool enableSource_ = true;
};

} // end namespace Dumux

#endif
