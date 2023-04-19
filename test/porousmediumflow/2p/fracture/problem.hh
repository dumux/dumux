// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPTests
 * \brief A discrete fracture network embedded in an impermeable matrix.
 *
 * The fracture is a 2D network embedded in 3D.
 */
#ifndef DUMUX_TWOP_FRACTURE_TEST_PROBLEM_HH
#define DUMUX_TWOP_FRACTURE_TEST_PROBLEM_HH

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/common/numeqvector.hh>

#include <dumux/porousmediumflow/problem.hh>

namespace Dumux {

/*!
 * \ingroup TwoPTests
 * \brief Trichloroethene (DNAPL) transport through a fracture network (2d in 3d).
 */
template <class TypeTag>
class FractureProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;

    enum
    {
        // primary variable indices
        pressureIdx = Indices::pressureIdx,
        saturationIdx = Indices::saturationIdx,

        // equation indices
        contiTCEEqIdx = Indices::conti0EqIdx + FluidSystem::comp1Idx,

        // world dimension
        dimWorld = GridView::dimensionworld
    };
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    FractureProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry) {}


    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary segment.
     *
     * \param globalPos The global position where to set the BC types
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;

        values.setAllDirichlet();
        if (onInlet_(globalPos))
            values.setAllNeumann();
        if (globalPos[2] > 1.0 - eps_ || globalPos[2] < eps_)
            values.setAllNeumann();

        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Dirichlet boundary segment.
     *
     * \param globalPos The global position
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition &globalPos) const
    {
        const auto depth = this->gridGeometry().bBoxMax()[dimWorld-1] - globalPos[dimWorld-1];
        const auto g = this->spatialParams().gravity(globalPos)[dimWorld-1];

        PrimaryVariables values;
        values[pressureIdx] = 1e5 + 1000*g*depth;
        values[saturationIdx] = 0.0;
        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Neumann boundary segment.
     *
     * \param globalPos The position of the integration point of the boundary segment.
     *
     * For this method, the \a values parameter stores the mass flux
     * in normal direction of each phase. Negative values mean influx.
     */
    NumEqVector neumannAtPos(const GlobalPosition &globalPos) const
    {
        NumEqVector values(0.0);
        if (onInlet_(globalPos)) {
            values[contiTCEEqIdx] = -0.04; // kg / (m * s)
        }
        return values;
    }

    /*!
     * \brief Evaluates the initial values for a control volume.
     *
     * \param globalPos The global position
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    { return dirichletAtPos(globalPos);}


private:

    bool onInlet_(const GlobalPosition &globalPos) const
    { return globalPos[0] < eps_ && globalPos[1] > -0.5 - eps_; }

    static constexpr Scalar eps_ = 1.5e-7;
    std::string name_;
};

} // end namespace Dumux

#endif
