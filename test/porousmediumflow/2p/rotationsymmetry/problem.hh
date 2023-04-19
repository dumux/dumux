// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//

#ifndef DUMUX_TEST_TWOP_ROTATIONALSYMMETRY_PROBLEM_HH
#define DUMUX_TEST_TWOP_ROTATIONALSYMMETRY_PROBLEM_HH

#include <dumux/common/properties.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/porousmediumflow/problem.hh>

namespace Dumux {

/*!
 * \brief A rotational symmetric 2p problem: water and air separate in a density-driven process in a dome shaped domain
 */
template<class TypeTag>
class TwoPRotationalSymmetryProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using BoundaryTypes = Dumux::BoundaryTypes<PrimaryVariables::size()>;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    TwoPRotationalSymmetryProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        FluidSystem::init(/*tempMin=*/293.0, /*tempMax=*/294.0, /*numTemp=*/2,
                          /*pMin=*/1.0e4, /*pMax=*/1.0e6, /*numP=*/200);
    }

    /*!
     * \brief The boundary types at position globalPos
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;
        values.setAllNeumann();
        return values;
    }

    /*!
     * \brief The initial values at position globalPos
     */
    PrimaryVariables initialAtPos(const GlobalPosition &globalPos) const
    {
        PrimaryVariables values;
        const Scalar depth = this->gridGeometry().bBoxMax()[1] - globalPos[1];
        // hydrostatic pressure
        values[0] = 1e5 - 1000*this->spatialParams().gravity(globalPos)[1]*depth;
        // start with saturation 0.5 -> density driven demixing
        values[1] = 0.5;
        return values;
    }


};

} // end namespace Dumux

#endif
