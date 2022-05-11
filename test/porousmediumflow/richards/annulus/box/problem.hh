// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
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
        rho_ = FluidSystem::H2O::liquidDensity(283.15, 1e5);
        mu_ = FluidSystem::H2O::liquidViscosity(283.15, 1e5);

        // geometry
        innerRadius_ = 0.0002;
        outerRadius_ = 0.006;

        // boundary conditions (volumetric flow rate in m^2/s, pressure/psi in Pa)
        const auto fi = getParam<Scalar>("Problem.InnerFlux", 0.1); // positive means outflow (cm/day)
        innerFlowRate_ = fi*2*M_PI*innerRadius_*0.01/(24*60*60);
        const auto fo = getParam<Scalar>("Problem.OuterFlux", 0.0); // positive means outflow (cm/day)
        outerFlowRate_ = fo*2*M_PI*outerRadius_*0.01/(24*60*60);

        const auto innerPressureHead = getParam<Scalar>("Problem.InnerPressureHeadInCm", -15000);
        const auto innerPressure = headToPressure(innerPressureHead);
        innerPsiCritical_ = applyKirchhoffTransformation_(innerPressure);
        source_ = (innerFlowRate_ + outerFlowRate_)/(M_PI*(innerRadius_*innerRadius_ - outerRadius_*outerRadius_));

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
        values.setAllNeumann();
        return values;
    }

    // positive means outflow
    NumEqVector neumannAtPos(const GlobalPosition& globalPos) const
    {
        if (globalPos[0] > this->gridGeometry().bBoxMax()[0] - 1e-10
            || globalPos[1] > this->gridGeometry().bBoxMax()[1] - 1e-10
            || globalPos[0] < this->gridGeometry().bBoxMin()[0] + 1e-10
            || globalPos[1] < this->gridGeometry().bBoxMin()[1] + 1e-10)
            return 0.0;
        else
            return { rho_*innerFlowRate_/(2*M_PI*innerRadius_) };
    }

    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    { return exactSolution(globalPos); }

    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    {
        PrimaryVariables values(initialPressure_);
        values.setState(Indices::bothPhases);
        return values;
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

        PrimaryVariables priVars(applyInverseKirchhoffTransformation_(psi));
        priVars.setState(Indices::bothPhases);
        return priVars;
    }

    PrimaryVariables exactSolution(const GlobalPosition& globalPos) const
    { return exactSolution(globalPos, innerPsiCritical_); }

    Scalar nonwettingReferencePressure() const
    { return 1.0e5; };

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
    Scalar innerFlowRate_, outerFlowRate_, innerPsiCritical_;
    Scalar initialPressure_;
};

} // end namespace Dumux

#endif
