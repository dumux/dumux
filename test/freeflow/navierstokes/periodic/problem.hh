// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesTests
 * \brief Stationary test for the staggered grid Navier-Stokes model with periodic BC
 */

#ifndef DUMUX_TEST_FREEFLOW_NAVIERSTOKES_PERIODIC_PROBLEM_HH
#define DUMUX_TEST_FREEFLOW_NAVIERSTOKES_PERIODIC_PROBLEM_HH

#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesTests
 * \brief  Periodic test problem for the staggered grid
 *
 * A two-dimensional Navier-Stokes flow with a periodicity in one direction
 * is considered.
 */
template <class TypeTag, class BaseProblem>
class PeriodicTestProblem : public BaseProblem
{
    using ParentType = BaseProblem;
    using BoundaryTypes = typename ParentType::BoundaryTypes;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using InitialValues = typename ParentType::InitialValues;
    using Sources = typename ParentType::Sources;
    using DirichletValues = typename ParentType::DirichletValues;
    using BoundaryFluxes = typename ParentType::BoundaryFluxes;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    static constexpr auto dimWorld = GridGeometry::GridView::dimensionworld;
    using Element = typename FVElementGeometry::Element;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

public:
    PeriodicTestProblem(std::shared_ptr<const GridGeometry> gridGeometry, std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry, couplingManager)
    {
        usePressureDifference_ = getParam<bool>("Problem.UsePressureDifference", false);
    }

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on a given boundary control volume.
     *
     * \param globalPos The position of the center of the finite volume
     */
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition &globalPos) const
    {
        BoundaryTypes values;
        values.setAllDirichlet();
        return values;
    }

    /*!
     * \brief Returns Dirichlet boundary values at a given position.
     *
     * \param globalPos The global position
     */
    DirichletValues dirichletAtPos(const GlobalPosition & globalPos) const
    { return DirichletValues(0.0); }

    template<class ElementVolumeVariables>
    Sources source(const Element& element,
                   const FVElementGeometry& fvGeometry,
                   const ElementVolumeVariables& elemVolVars,
                   const SubControlVolume& scv) const
    {
        Sources source;

        if constexpr (ParentType::isMomentumProblem())
        {
            if (usePressureDifference_ && scv.dofPosition()[1] < this->gridGeometry().bBoxMin()[1] + eps_)
            {
                const auto& frontalScvf = (*scvfs(fvGeometry, scv).begin());
                source[Indices::momentumYBalanceIdx] = 100 * frontalScvf.area() / scv.volume();
            }
        }

        return source;
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    //! Enable internal Dirichlet constraints
    static constexpr bool enableInternalDirichletConstraints()
    { return !ParentType::isMomentumProblem(); }

    /*!
     * \brief Tag a degree of freedom to carry internal Dirichlet constraints.
     *        If true is returned for a dof, the equation for this dof is replaced
     *        by the constraint that its primary variable values must match the
     *        user-defined values obtained from the function internalDirichlet(),
     *        which must be defined in the problem.
     *
     * \param element The finite element
     * \param scv The sub-control volume
     */
    std::bitset<DirichletValues::dimension> hasInternalDirichletConstraint(const Element& element, const SubControlVolume& scv) const
    {
        std::bitset<DirichletValues::dimension> values;

        for (const auto& intersection : intersections(this->gridGeometry().gridView(), element))
        {
            const auto center = intersection.geometry().center();
            if (intersection.boundary() && intersection.neighbor() && center[1] > this->gridGeometry().bBoxMax()[1] - eps_)
                values.set(0);
        }

        return values;
    }

    /*!
     * \brief Define the values of internal Dirichlet constraints for a degree of freedom.
     * \param element The finite element
     * \param scv The sub-control volume
     */
    DirichletValues internalDirichlet(const Element& element, const SubControlVolume& scv) const
    {
        return DirichletValues(1.0);
    }

private:
    static constexpr Scalar eps_ = 1e-6;
    bool usePressureDifference_;
};

} // end namespace Dumux

#endif
