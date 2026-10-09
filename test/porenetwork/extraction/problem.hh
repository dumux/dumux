// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup PoreNetworkModels
 * \brief Single-phase flow through the pore network of a sphere packing between an inlet and an outlet
 */
#ifndef DUMUX_PNM_EXTRACTION_PERMEABILITY_PROBLEM_HH
#define DUMUX_PNM_EXTRACTION_PERMEABILITY_PROBLEM_HH

#include <dumux/common/boundarytypes.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/porousmediumflow/problem.hh>

namespace Dumux {

template <class TypeTag>
class SpherePackingPermeabilityProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;

public:
    template<class SpatialParams>
    SpherePackingPermeabilityProblem(std::shared_ptr<const GridGeometry> gridGeometry, std::shared_ptr<SpatialParams> spatialParams)
    : ParentType(gridGeometry, spatialParams)
    {
        inletLabel_ = getParam<int>("Problem.InletLabel");
        outletLabel_ = getParam<int>("Problem.OutletLabel");
        inletPressure_ = getParam<Scalar>("Problem.InletPressure");
        outletPressure_ = getParam<Scalar>("Problem.OutletPressure");
    }

    BoundaryTypes boundaryTypes(const Element& element, const SubControlVolume& scv) const
    {
        BoundaryTypes bcTypes;
        if (isInlet_(scv) || isOutlet_(scv))
            bcTypes.setAllDirichlet();
        else
            bcTypes.setAllNeumann();
        return bcTypes;
    }

    PrimaryVariables dirichlet(const Element& element, const SubControlVolume& scv) const
    { return PrimaryVariables(isInlet_(scv) ? inletPressure_ : outletPressure_); }

    PrimaryVariables initial(const Element&, const SubControlVolume&) const
    { return PrimaryVariables(outletPressure_); }

    int inletLabel() const { return inletLabel_; }
    int outletLabel() const { return outletLabel_; }

private:
    bool isInlet_(const SubControlVolume& scv) const
    { return this->gridGeometry().poreLabel(scv.dofIndex()) == inletLabel_; }

    bool isOutlet_(const SubControlVolume& scv) const
    { return this->gridGeometry().poreLabel(scv.dofIndex()) == outletLabel_; }

    int inletLabel_, outletLabel_;
    Scalar inletPressure_, outletPressure_;
};

} // end namespace Dumux

#endif
