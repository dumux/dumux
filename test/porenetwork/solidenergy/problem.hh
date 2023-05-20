// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef DUMUX_TEST_PORENETWORK_SOLID_ENERGY_PROBLEM_HH
#define DUMUX_TEST_PORENETWORK_SOLID_ENERGY_PROBLEM_HH

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/porousmediumflow/problem.hh>

namespace Dumux {

/*!
 * \brief Heat conduction problem with multiple solid spheres
 */
template <class TypeTag>
class SolidProblem : public PorousMediumFlowProblem<TypeTag>
{
    using ParentType = PorousMediumFlowProblem<TypeTag>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;

    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    template<class SpatialParams>
    SolidProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                    std::shared_ptr<SpatialParams> spatialParams)
    : ParentType(gridGeometry, spatialParams)
    {
        problemName_  =  getParam<std::string>("Problem.Name");
        initialTemperature_ = getParam<Scalar>("Problem.InitialTemperature");
        temperatureLeft_ = getParam<Scalar>("Problem.LeftTemperature");

        leftIndex_ = getParam<int>("Problem.LeftIndex");
        rightIndex_ = getParam<int>("Problem.RightIndex");
    }

    const std::string& name() const
    { return problemName_; }

    BoundaryTypes boundaryTypes(const Element &element, const SubControlVolume& scv) const
    {
        BoundaryTypes values;
        values.setAllNeumann();

        if (onLeftBoundary_(scv) || onRightBoundary_(scv))
            values.setAllDirichlet();

        return values;
    }

    template<class ElementVolumeVariables>
    NumEqVector source(const Element& element,
                       const FVElementGeometry& fvGeometry,
                       const ElementVolumeVariables& elemVolVars,
                       const SubControlVolume& scv) const
    {
        NumEqVector value = 0.0; //isolating boundary condition
        return value;
    }

    PrimaryVariables dirichlet(const Element& element, const SubControlVolume& scv) const
    {
        auto values = initialAtPos(scv.dofPosition()); //onRightBoundary_(scv)

        if (onLeftBoundary_(scv))
        {
            values = temperatureLeft_;
        }

        return values;
    }

    PrimaryVariables initialAtPos(const GlobalPosition& pos) const
    {
        PrimaryVariables values(initialTemperature_); //uniform initial temperature
        return values;
    }

private:

    bool onLeftBoundary_(const SubControlVolume& scv) const
    { return this->gridGeometry().poreLabel(scv.dofIndex()) == leftIndex_; }

    bool onRightBoundary_(const SubControlVolume& scv) const
    { return this->gridGeometry().poreLabel(scv.dofIndex()) == rightIndex_; }

    std::string problemName_;

    Scalar initialTemperature_;
    Scalar temperatureLeft_;

    int leftIndex_;
    int rightIndex_;
};

} // end namespace Dumux

#endif
