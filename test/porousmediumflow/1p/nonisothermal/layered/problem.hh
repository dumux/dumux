// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePTests
 * \brief Heat conduction through two layers of different thermal conductivity
 */
#ifndef DUMUX_ONEPNI_LAYERED_PROBLEM_HH
#define DUMUX_ONEPNI_LAYERED_PROBLEM_HH

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
 * \brief Heat conduction through two layers of different thermal conductivity
 *
 * The fluid is at rest and the heat flux \f$ q \f$ points in the x-direction. The temperature is
 * \f$ T_0 \f$ at \f$ x = 0 \f$ and linear in each layer with the gradient \f$ -q/\lambda_\mathrm{eff} \f$
 * of the layer, where \f$ \lambda_\mathrm{eff} = \phi \lambda_f + (1 - \phi) \lambda_s \f$.
 * The whole boundary carries the exact pressure and temperature.
 */
template<class TypeTag>
class OnePNILayeredProblem : public Experimental::ProblemWithSpatialParams<TypeTag>
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
    OnePNILayeredProblem(std::shared_ptr<const GridDiscretization> gridDiscretization)
    : ParentType(gridDiscretization)
    , pressure_(getParam<Scalar>("Problem.Pressure"))
    , temperatureLeft_(getParam<Scalar>("Problem.TemperatureLeft"))
    , heatFlux_(getParam<Scalar>("Problem.HeatFlux"))
    {
        const auto fluidConductivity = getParam<Scalar>("Component.LiquidThermalConductivity");
        const auto solidConductivity = getParam<Scalar>("Component.SolidThermalConductivity");
        const auto conductivity = [&] (Scalar porosity)
        { return porosity*fluidConductivity + (1.0 - porosity)*solidConductivity; };
        conductivityLeft_ = conductivity(this->spatialParams().porosityLeft());
        conductivityRight_ = conductivity(this->spatialParams().porosityRight());

        appendDirichletConstraints_();
    }

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        // the whole boundary is Dirichlet, which the assembler reads off the constraints;
        // a boundary type only marks the equations that receive a flux instead
        return BoundaryTypes{};
    }

    NumEqVector sourceAtPos(const GlobalPosition& globalPos) const
    { return NumEqVector(0.0); }

    Scalar pressure() const
    { return pressure_; }

    Scalar heatFlux() const
    { return heatFlux_; }

    Scalar exactTemperature(const GlobalPosition& globalPos) const
    {
        const Scalar x = globalPos[0];
        const Scalar xi = this->spatialParams().interfacePosition();
        if (x < xi)
            return temperatureLeft_ - heatFlux_*x/conductivityLeft_;
        else
            return temperatureLeft_ - heatFlux_*xi/conductivityLeft_ - heatFlux_*(x - xi)/conductivityRight_;
    }

    //! The Dirichlet values the assembler constrains the solution to
    const std::vector<DirichletConstraintData>& constraints() const
    { return constraints_; }

private:
    //! Every local dof on the boundary carries the exact pressure and temperature
    void appendDirichletConstraints_()
    {
        appendBoundaryDirichletConstraints(
            constraints_, this->gridDiscretization(),
            [] (const auto& elemDisc, const auto& boundaryFace, const auto& localDof)
            {
                ConstraintInfo info;
                info.setAll();
                return info;
            },
            [&] (const auto& elemDisc, const auto& localDof)
            {
                ConstraintValues values(0.0);
                values[Indices::pressureIdx] = pressure_;
                values[Indices::temperatureIdx] = exactTemperature(ipData(elemDisc, localDof).global());
                return values;
            }
        );
    }

    Scalar pressure_;
    Scalar temperatureLeft_;
    Scalar heatFlux_;
    Scalar conductivityLeft_;
    Scalar conductivityRight_;
    std::vector<DirichletConstraintData> constraints_;
};

} // end namespace Dumux

#endif
