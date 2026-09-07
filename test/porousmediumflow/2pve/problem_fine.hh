// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup VETest
 * \brief The immiscible 2p fine-level VE test problem.
 */

#ifndef DUMUX_TEST_TWOPVE_FINE_PROBLEM_HH
#define DUMUX_TEST_TWOPVE_FINE_PROBLEM_HH

#include <memory>

#include <dumux/common/numeqvector.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

#include "spatialparams_fine.hh"

namespace Dumux {

template<class TypeTag>
class TwoPVEFineProblem
{
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    enum {
        pressureH2OIdx = Indices::pressureIdx,
        saturationGasIdx = Indices::saturationIdx,
        contiGasEqIdx = Indices::conti0EqIdx + FluidSystem::comp1Idx,
    };
    enum {
        dim = GridView::dimension,
    };
    using WettingPhase = typename GetProp<TypeTag, Properties::FluidSystem>::WettingPhase;

public:
    using SpatialParamsFine = TwoPTestFineSpatialParams<GridGeometry, Scalar>;

    TwoPVEFineProblem(std::shared_ptr<const GridGeometry> gridGeometryFine,
                      const Scalar& fineCellHeight,
                      const Scalar& fineCellDepth)
    : spatialParams_(std::make_shared<SpatialParamsFine>(gridGeometryFine)),
      fineCellHeight_(fineCellHeight),
      fineCellDepth_(fineCellDepth)
    {
        injectionRate_ = getParam<double>("BoundaryConditions.InjectionRate");
    }

    /*!
     * \brief Return a reference to the fine-level spatial parameters
     */
    const SpatialParamsFine& spatialParams() const
    {
        return *spatialParams_;
    }

    /*!
     * \brief Return a pointer to the fine-level spatial parameters
     */
    std::shared_ptr<const SpatialParamsFine> spatialParamsPtr() const
    {
        return spatialParams_;
    }

    /*!
     * \brief Evaluates the initial values for a control volume of the fine-level grid. Initial conditions only need to be changed here.
     *
     * \param globalPos the global position of the scv belonging to the fine-level grid
     */
    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    {
        PrimaryVariables values(0.0);

        Scalar temperature = 326.0;
        Scalar pressureWInit = 1.0e7;
        Scalar densityW = WettingPhase::density(temperature,pressureWInit);

        const auto z = globalPos[dim-1] - spatialParams_->gridGeometry().bBoxMin()[dim-1];
        values[pressureH2OIdx] = pressureWInit - densityW * spatialParams_->gravity(globalPos).two_norm() * z;
        values[saturationGasIdx] = 0.0;

        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Neumann boundary segment on the fine-level grid. Neumann conditions only need to be changed here.
     *
     * \param globalPos the position of the integration point of the boundary segment belonging to the fine-level grid
     */
    NumEqVector neumannAtPos(const GlobalPosition& globalPos) const
    {
        if constexpr(dim==2)
        {
            NumEqVector values(0.0);
            if (onLeftBoundary_(globalPos))
                values[contiGasEqIdx] = injectionRate_; // kg / (m * s)
            return values;
        }
        else if constexpr(dim==3)
        {
            NumEqVector values(0.0);
            // put well at middle of y-direction boundary
            Scalar halfDepthBox = (spatialParams_->gridGeometry().bBoxMax()[1] - spatialParams_->gridGeometry().bBoxMin()[1])/2.0;

            if(onLeftBoundary_(globalPos) &&
               halfDepthBox <= (globalPos[1]+fineCellDepth_/2.0) &&
               halfDepthBox >  (globalPos[1]-fineCellDepth_/2.0) )
            {
                values[contiGasEqIdx] = injectionRate_/fineCellDepth_;
            }

            return values;
        }
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet boundary segment on the fine-level grid. Dirichlet conditions only need to be changed here.
     *
     * \param globalPos the global position of an element belonging to the fine level
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    {
        PrimaryVariables values;

        Scalar temperature = 326.0;
        Scalar pressureWInit = 1.0e7;
        Scalar densityW = WettingPhase::density(temperature,pressureWInit);

        const auto z = globalPos[dim-1] - spatialParams_->gridGeometry().bBoxMin()[dim-1];
        values[pressureH2OIdx] = pressureWInit - densityW * spatialParams_->gravity(globalPos).two_norm() * z;
        values[saturationGasIdx] = 0.0;

        return values;
    }


private:
    static constexpr Scalar eps_ = 1e-6;

    bool onLeftBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[0] < spatialParams_->gridGeometry().bBoxMin()[0] + eps_;
    }

    bool onRightBoundary_(const GlobalPosition &globalPos) const
    {
        return globalPos[0] > spatialParams_->gridGeometry().bBoxMax()[0] - eps_;
    }

    std::shared_ptr<const SpatialParamsFine> spatialParams_;
    Scalar fineCellHeight_, fineCellDepth_;
    Scalar injectionRate_;
};

} // end namespace Dumux

#endif
