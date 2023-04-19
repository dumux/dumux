// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PorousmediumflowModels
 * \brief Base class for all porous media problems.
 */
#ifndef DUMUX_POROUS_MEDIUM_FLOW_PROBLEM_HH
#define DUMUX_POROUS_MEDIUM_FLOW_PROBLEM_HH

#include <dumux/common/properties.hh>
#include <dumux/common/fvproblemwithspatialparams.hh>

namespace Dumux {

/*!
 * \ingroup PorousmediumflowModels
 * \brief Base class for all fully implicit porous media problems.
 *
 * TODO: derive from base problem property?
 */
template<class TypeTag>
class PorousMediumFlowProblem : public FVProblemWithSpatialParams<TypeTag>
{
    using ParentType = FVProblemWithSpatialParams<TypeTag>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:

    //! Use constructors of the base class
    using ParentType::ParentType;
};

} // end namespace Dumux

#endif
