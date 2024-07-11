// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup ???
 * \brief The problem for the trace test
 */

#ifndef DUMUX_TRACE_PROBLEM_HH
#define DUMUX_TRACE_PROBLEM_HH

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>

#include <dumux/common/numeqvector.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/porousmediumflow/problem.hh>

namespace Dumux {

template<typename GlobalPosition>
double pressureField(const GlobalPosition& pos)
{ return pos[0]*pos[1]; }

template<typename GlobalPosition>
GlobalPosition pressureFieldGradient(const GlobalPosition& pos)
{
    GlobalPosition result;
    result[0] = pos[1];
    result[1] = pos[0];
    return result;
}

/*!
 * \ingroup ???
 * \brief The problem for the trace test
 */
template<class TypeTag>
class TraceTestProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;

    using GridView = typename GridGeometry::GridView;
    static constexpr int dimWorld = GridView::dimensionworld;
    using GlobalPosition = Dune::FieldVector<typename GridView::ctype, dimWorld>;

public:
    TraceTestProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {}

    //! Specifies which kind of boundary condition to use at the given position.
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        values.setAllDirichlet();
        if (globalPos[0] < 1e-6)
            values.setAllNeumann();
        return values;
    }

    //! Evaluates the Dirichlet boundary conditions.
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    { return pressureField(globalPos); }

    //! Evaluates the Neumann boundary conditions.
    template<typename ElementVolumeVariables, typename ElementFluxVariablesCache>
    NumEqVector<PrimaryVariables> neumann(const typename GridView::template Codim<0>::Entity& element,
                                          const typename GridGeometry::LocalView& fvGeometry,
                                          const ElementVolumeVariables& elemVolVars,
                                          const ElementFluxVariablesCache& elemFluxVarsCache,
                                          const typename GridGeometry::SubControlVolumeFace& scvf) const
    {
        return -1.0*this->spatialParams().permeabilityAtPos(scvf.ipGlobal())*(
            pressureFieldGradient(scvf.ipGlobal())
            *scvf.unitOuterNormal()
        );
    }
};

} // end namespace Dumux


#endif
