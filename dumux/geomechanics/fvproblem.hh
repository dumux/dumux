// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup GeomechanicsModels
 * \brief Base class for all geomechanical problems
 */
#ifndef DUMUX_GEOMECHANICS_FV_PROBLEM_HH
#define DUMUX_GEOMECHANICS_FV_PROBLEM_HH

#include <dumux/common/typetraits/isvalid.hh>
#include <dumux/porousmediumflow/problem.hh>

namespace Dumux {

/*!
 * \ingroup GeomechanicsModels
 * \brief Base class for all geomechanical problems
 * \note We require only little additional functionality to the
 *       porous medium flow problem, which is why we inherit from that here.
 */
template<class TypeTag>
class GeomechanicsFVProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    static constexpr int numFP = GetPropType<TypeTag, Properties::ModelTraits>::numFluidPhases();

public:
    //! pull up the constructor of the parent class
    using ParentType::ParentType;
};

} // end namespace Dumux

#endif
