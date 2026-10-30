// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePTests
 * \brief Convergence test of the non-isothermal single-phase model with a manufactured solution
 */
#ifndef DUMUX_ONEPNI_CONVERGENCE_PROBLEM_HH
#define DUMUX_ONEPNI_CONVERGENCE_PROBLEM_HH

#include <cmath>
#include <vector>

#include <dune/common/fvector.hh>

#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/problemwithspatialparams.hh>
#include <dumux/common/constraintinfo.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/discretization/dirichletconstraints.hh>

namespace Dumux {

/*!
 * \ingroup OnePTests
 * \brief Convergence test of the non-isothermal single-phase model with a manufactured solution
 *
 * Stationary mass and energy balance on the unit square,
 * \f[
 *   \nabla \cdot (\varrho \mathbf{v}) = q_m, \quad
 *   \nabla \cdot (\varrho h \mathbf{v} - \lambda_\mathrm{eff} \nabla T) = q_e, \quad
 *   \mathbf{v} = -\frac{k}{\mu} \nabla p,
 * \f]
 * with the affine pressure \f$ p = p_0 + \mathbf{a} \cdot \mathbf{x} \f$, the temperature
 * \f$ T = T_0 + \Delta T \sin(\pi x) \sin(\pi y) \f$, the enthalpy \f$ h = c (T - T_\mathrm{ref}) + p/\varrho \f$
 * of an incompressible fluid and the effective conductivity
 * \f$ \lambda_\mathrm{eff} = \phi \lambda_f + (1 - \phi) \lambda_s \f$ of a porosity varying in space.
 * Density, viscosity and heat capacity are constant, so the velocity is uniform and
 * \f$ q_m = 0 \f$, \f$ q_e = \varrho c \mathbf{v} \cdot \nabla T + \mathbf{v} \cdot \nabla p - \nabla \cdot (\lambda_\mathrm{eff} \nabla T) \f$.
 * With \f$ \mathbf{a} = 0 \f$ the heat is transported by conduction only. The exact solution is
 * prescribed on the whole boundary.
 */
template<class TypeTag>
class OnePNIConvergenceProblem : public Experimental::ProblemWithSpatialParams<TypeTag>
{
    using ParentType = Experimental::ProblemWithSpatialParams<TypeTag>;

    using GridDiscretization = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridDiscretization::GridView;
    using GlobalPosition = typename GridView::template Codim<0>::Entity::Geometry::GlobalCoordinate;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename ModelTraits::Indices;
    using BoundaryTypes = Dumux::Experimental::BoundaryTypes<ModelTraits::numEq()>;

    using ConstraintInfo = Dumux::DirichletConstraintInfo<ModelTraits::numEq()>;
    using ConstraintValues = Dune::FieldVector<Scalar, ModelTraits::numEq()>;
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using DirichletConstraintData = Dumux::DirichletConstraintData<ConstraintInfo, ConstraintValues, GridIndexType>;

public:
    OnePNIConvergenceProblem(std::shared_ptr<const GridDiscretization> gridDiscretization)
    : ParentType(gridDiscretization)
    , pressureGradient_(getParam<GlobalPosition>("Problem.PressureGradient"))
    , referencePressure_(getParam<Scalar>("Problem.ReferencePressure"))
    , referenceTemperature_(getParam<Scalar>("Problem.ReferenceTemperature"))
    , temperatureAmplitude_(getParam<Scalar>("Problem.TemperatureAmplitude"))
    , density_(getParam<Scalar>("Component.LiquidDensity"))
    , viscosity_(getParam<Scalar>("Component.LiquidKinematicViscosity")*density_)
    , heatCapacity_(getParam<Scalar>("Component.LiquidHeatCapacity"))
    , fluidConductivity_(getParam<Scalar>("Component.LiquidThermalConductivity"))
    , solidConductivity_(getParam<Scalar>("Component.SolidThermalConductivity"))
    { appendDirichletConstraints_(); }

    //! The whole boundary is Dirichlet, given by the constraints
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    { return BoundaryTypes{}; }

    NumEqVector sourceAtPos(const GlobalPosition& globalPos) const
    {
        const auto x = globalPos[0];
        const auto y = globalPos[1];
        using std::sin; using std::cos;
        const Scalar pi = M_PI;

        // the temperature gradient and Laplacian
        GlobalPosition gradT(0.0);
        gradT[0] = temperatureAmplitude_*pi*cos(pi*x)*sin(pi*y);
        gradT[1] = temperatureAmplitude_*pi*sin(pi*x)*cos(pi*y);
        const Scalar laplaceT = -2.0*pi*pi*(exactTemperature(globalPos) - referenceTemperature_);

        // the effective conductivity and its gradient
        const auto& spatialParams = this->spatialParams();
        const Scalar porosity = spatialParams.porosityAtPos(globalPos);
        GlobalPosition gradPorosity(0.0);
        gradPorosity[0] = porosity/(1.0 + x*y)*y;
        gradPorosity[1] = porosity/(1.0 + x*y)*x;
        const Scalar lambda = porosity*fluidConductivity_ + (1.0 - porosity)*solidConductivity_;
        auto gradLambda = gradPorosity;
        gradLambda *= fluidConductivity_ - solidConductivity_;

        auto velocity = pressureGradient_;
        velocity *= -spatialParams.permeabilityAtPos(globalPos)/viscosity_;

        NumEqVector source(0.0);
        source[Indices::energyEqIdx] = density_*heatCapacity_*(velocity*gradT) + velocity*pressureGradient_
                                       - (gradLambda*gradT + lambda*laplaceT);
        return source;
    }

    Scalar exactPressure(const GlobalPosition& globalPos) const
    { return referencePressure_ + pressureGradient_*globalPos; }

    Scalar exactTemperature(const GlobalPosition& globalPos) const
    {
        using std::sin;
        return referenceTemperature_ + temperatureAmplitude_*sin(M_PI*globalPos[0])*sin(M_PI*globalPos[1]);
    }

    //! The Dirichlet values the assembler constrains the solution to
    const std::vector<DirichletConstraintData>& constraints() const
    { return constraints_; }

private:
    //! Every local dof on the boundary carries the exact pressure and temperature
    void appendDirichletConstraints_()
    {
        auto elemDisc = localView(this->gridDiscretization());
        for (const auto& element : elements(this->gridDiscretization().gridView()))
        {
            elemDisc.bind(element);

            for (const auto& boundaryFace : boundaryFaces(elemDisc))
                for (const auto& localDof : localDofs(elemDisc, boundaryFace))
                {
                    const auto globalPos = ipData(elemDisc, localDof).global();

                    ConstraintInfo info;
                    info.setAll();

                    ConstraintValues values(0.0);
                    values[Indices::pressureIdx] = exactPressure(globalPos);
                    values[Indices::temperatureIdx] = exactTemperature(globalPos);
                    constraints_.push_back(
                        DirichletConstraintData{std::move(info), std::move(values), localDof.dofIndex()}
                    );
                }
        }
    }

    GlobalPosition pressureGradient_;
    Scalar referencePressure_;
    Scalar referenceTemperature_;
    Scalar temperatureAmplitude_;
    Scalar density_;
    Scalar viscosity_;
    Scalar heatCapacity_;
    Scalar fluidConductivity_;
    Scalar solidConductivity_;
    std::vector<DirichletConstraintData> constraints_;
};

} // end namespace Dumux

#endif
