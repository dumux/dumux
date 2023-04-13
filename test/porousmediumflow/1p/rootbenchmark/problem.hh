// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePTests
 * \brief Root benchmark case Schnepf et al 2020 M3.1 doi: 10.3389/fpls.2020.00316
 */

#ifndef DUMUX_ROOTBENCHMARK_TEST_PROBLEM_HH
#define DUMUX_ROOTBENCHMARK_TEST_PROBLEM_HH

#include <dune/common/fvector.hh>
#include <dune/common/fmatrix.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/common/numeqvector.hh>

#include <dumux/io/format.hh>

#include <dumux/porousmediumflow/problem.hh>

#include <dumux/material/components/constant.hh>

namespace Dumux {

/*!
 * \ingroup OnePTests
 * \brief Root benchmark case Schnepf et al 2020 M3.1 doi: 10.3389/fpls.2020.00316
 */
template <class TypeTag>
class RootBenchmarkProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    static constexpr int dimWorld = GridView::dimensionworld;

    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;

public:
    RootBenchmarkProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        name_ = getParam<std::string>("Problem.Name");
        rootCollarPressure_ = getParam<double>("Problem.RootCollarPressure");
        soilPressure_ = getParam<double>("Problem.SoilPressure");
        radius_ = getParam<double>("Problem.Radius");
        density_ = getParam<double>("Component.LiquidDensity");
        Kr_ = this->spatialParams().radialConductivity();
        Kx_ = this->spatialParams().axialConductivity();
        const auto enableGravity = getParam<bool>("Problem.EnableGravity");
        const Scalar gravity = enableGravity ? 9.81 : 0.0;

        // cf. Schnepf et al 2020 doi: 10.3389/fpls.2020.00316
        // Equation is 2*pi*R*Kr(ps - p0) = d/ds(-Kx (dp/ds + rho*g))
        // Analytical solution is p = ps + A*exp(Cs) + B*exp(-Cs)
        // where C = sqrt(2*pi*R*Kr/Kx)
        // Two boundary conditions are needed to determine A and B
        // Dirichlet at root collar (s=0, p=p0)
        // (1) A + B = p0 - ps
        // Neumann no-flow at tip (s=-L, dp/ds + rho*g = 0)
        // (2) A*C*exp(-CL) - B*C*exp(CL) = -rho*g
        // Assemble (1) and (2) and solve for A and B
        Dune::FieldMatrix<double, 2, 2> M;
        Dune::FieldVector<double, 2> b;
        using std::exp; using std::sqrt;
        const auto C = sqrt(2.0*M_PI*radius_*Kr_/Kx_); C_ = C;
        const auto L = this->gridGeometry().bBoxMax()[dimWorld-1] - this->gridGeometry().bBoxMin()[dimWorld-1];
        M[0][0] = 1.0; M[0][1] = 1.0; b[0] = rootCollarPressure_ - soilPressure_;
        M[1][0] = C*exp(-C*L); M[1][1] = -C*exp(C*L); b[1] = -density_*gravity;
        M.solve(d_, b);

        std::cout << "Computed analytical solution with:\n";
        std::cout << Fmt::format("A or d1 = {}, B or d2 = {}, c = {}, √c = {}\n", d_[0], d_[1], C*C, C);
        std::cout << Fmt::format("Kx = {}, Kr = {}\n", Kx_, Kr_);

        const auto psiCollar = 100*(rootCollarPressure_ - 1e5)/(1000*9.81);
        const auto psiSoil = 100*(soilPressure_ - 1e5)/(1000*9.81);
        const auto psi0 = 100*(analyticalSolution(gridGeometry->bBoxMax()) - 1e5)/(1000*9.81);
        const auto psiTip = 100*(analyticalSolution(gridGeometry->bBoxMin()) - 1e5)/(1000*9.81);
        std::cout << Fmt::format("ψ(0) = {}, ψ(-0.5) = {}, ψ0 = {}, ψs = {} cm\n", psi0, psiTip, psiCollar, psiSoil);

        analyticalPressures_.resize(gridGeometry->gridView().size(0));
        analyticalHead_.resize(gridGeometry->gridView().size(0));
        for (const auto& element : elements(gridGeometry->gridView()))
        {
            const auto eIdx = gridGeometry->elementMapper().index(element);
            const auto& globalPos = element.geometry().center();
            analyticalPressures_[eIdx] = analyticalSolution(globalPos);
            analyticalHead_[eIdx] = 100*(analyticalPressures_[eIdx] - 1e5)/(1000*9.81);
        }
    }

    /*!
     * \brief The problem name.
     */
    const std::string& name() const
    { return name_; }

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary control volume.
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes bcTypes;
        if (globalPos[dimWorld-1] > this->gridGeometry().bBoxMax()[dimWorld-1] - eps_)
            bcTypes.setAllDirichlet(); // Dirichlet at root collar
        else
            bcTypes.setAllNeumann(); // No-flow at root tips
        return bcTypes;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet control volume.
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(rootCollarPressure_); }

    /*!
     * \brief Evaluates the source term for all phases within a given sub-control volume.
     * Positive values mean that the conserved quantity is created, negative ones mean that it vanishes.
     * E.g. for the mass balance that would be a mass rate in \f$ [ kg / (m^3 \cdot s)] \f$.
     */
    template<class ElementVolumeVariables>
    NumEqVector source(const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolume& scv) const
    {
        NumEqVector values(0.0);
        const auto& volVars = elemVolVars[scv];
        // source term with static soil pressure
        values[0] = density_*2.0*M_PI*radius_*Kr_*(soilPressure_ - volVars.pressure());
        // make this volume-specific source
        values[0] /= volVars.extrusionFactor();
        return values;
    }

    template<class SolutionVector>
    void outputL2Norm(const SolutionVector& solution) const
    {
        // calculate the discrete L2-Norm and maximum element size
        Scalar l2Norm = 0.0;
        Scalar l2NormNormalization = 0.0;
        Scalar hMax = 0.0;

        const auto& gg = this->gridGeometry();
        for (const auto& element : elements(gg.gridView()))
        {
            const auto eIdx = gg.elementMapper().index(element);
            const auto geo = element.geometry();
            const auto globalPos = geo.center();
            const auto pe = analyticalSolution(globalPos);
            const auto p = solution[eIdx][0];
            const auto dV = geo.volume();

            l2Norm += (p - pe)*(p - pe)*dV;
            l2NormNormalization += pe*pe*dV;
            hMax = std::max(dV, hMax);
        }

        l2Norm = std::sqrt(l2Norm/l2NormNormalization);

        // write the norm into a log file
        std::ofstream logFile(this->name() + ".log", std::ios::app);
        logFile << "[ConvergenceTest] L2-norm(pressure) = " << l2Norm  << " hMax = " << hMax << std::endl;
    }

    Scalar analyticalSolution(const GlobalPosition& globalPos) const
    {
        using std::exp;
        const auto s = globalPos[dimWorld-1];
        return soilPressure_ + d_[0]*exp(s*C_) + d_[1]*exp(-s*C_);
    }

    const std::vector<Scalar>& analyticalPressure() const
    { return analyticalPressures_; }

    const std::vector<Scalar>& analyticalHead() const
    { return analyticalHead_; }

private:
    Scalar rootCollarPressure_, soilPressure_;
    Scalar radius_;
    Scalar density_;
    Scalar Kr_, Kx_;

    static constexpr Scalar eps_ = 1e-8;
    std::string name_;

    // analytical solution
    Dune::FieldVector<Scalar, 2> d_;
    Scalar C_;

    std::vector<Scalar> analyticalPressures_;
    std::vector<Scalar> analyticalHead_;
};

} // end namespace Dumux

#endif
