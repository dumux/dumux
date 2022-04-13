// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup NavierStokesNCTests
 * \brief Channel flow test for the multi-component staggered grid (Navier-)Stokes model.
 */

#ifndef DUMUX_COLUMN_1PNC_TEST_TEST_PROBLEM_HH
#define DUMUX_COLUMN_1PNC_TEST_TEST_PROBLEM_HH

#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

#include <dumux/freeflow/navierstokes/problem.hh>

#include <dumux/freeflow/navierstokes/momentum/fluxhelper.hh>
#include <dumux/freeflow/navierstokes/scalarfluxhelper.hh>

#include <dumux/io/gnuplotinterface.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesNCTests
 * \brief Test problem for the Maxwell-Stefan model
 */
template <class TypeTag>
class ColumnNCTestProblem : public NavierStokesProblem<TypeTag>
{
    using ParentType = NavierStokesProblem<TypeTag>;

    using BoundaryTypes = typename ParentType::BoundaryTypes;
    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using BoundaryFluxes = typename ParentType::BoundaryFluxes;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using InitialValues = typename ParentType::InitialValues;
    using DirichletValues = typename ParentType::DirichletValues;

    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

     using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

    static constexpr auto compOneIdx =  Indices::conti0EqIdx;
    static constexpr auto transportCompIdx = Indices::conti0EqIdx + 1;
    static constexpr auto transportEqIdx = Indices::conti0EqIdx + 1;

public:
    ColumnNCTestProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                        std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry, couplingManager)
    {
        FluidSystem::init();
    }

    /*!
     * \brief Returns a reference pressure at a given sub control volume face.
     *        This pressure is substracted from the actual pressure for the momentum balance
     *        which potentially helps to improve numerical accuracy by avoiding issues related do floating point arithmetic.
     */
    Scalar referencePressure(const Element& element,
                             const FVElementGeometry& fvGeometry,
                             const SubControlVolumeFace& scvf) const
    { return 1.0e5; }

    // \}
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

        if constexpr (ParentType::isMomentumProblem())
        {
            values.setDirichlet(Indices::velocityXIdx);
            values.setDirichlet(Indices::velocityYIdx);

           //if (onTopBoundary_(globalPos))
           //    values.setAllNeumann();
        }
        else
        {
            values.setAllNeumann();

        }

        return values;
    }

    /*!
     * \brief Evaluates the boundary conditions for a Neumann control volume.
     *
     * \param element The element for which the Neumann boundary condition is set
     * \param fvGeometry The fvGeometry
     * \param elemVolVars The element volume variables
     * \param elemFaceVars The element face variables
     * \param scvf The boundary sub control volume face
     */
    template<class ElementVolumeVariables, class ElementFluxVariablesCache>
    BoundaryFluxes neumann(const Element& element,
                           const FVElementGeometry& fvGeometry,
                           const ElementVolumeVariables& elemVolVars,
                           const ElementFluxVariablesCache& elemFluxVarsCache,
                           const SubControlVolumeFace& scvf) const
    {
        BoundaryFluxes values(0.0);
//         const auto& ipGlobal = scvf.ipGlobal();
//         if (onTopBoundary_(ipGlobal))
//         {
//             if constexpr (ParentType::isMomentumProblem())
//             {
//                 using FluxHelper = NavierStokesMomentumBoundaryFluxHelper;                 values = FluxHelper::fixedPressureMomentumFlux(*this, fvGeometry, scvf,
//                                                               elemVolVars, elemFluxVarsCache,
//                                                               referencePressure(element, fvGeometry, scvf), true /*zeroNormalVelocityGradient*/);
//             }
//             else
//             {
//                 const auto insideDensity = elemVolVars[scvf.insideScvIdx()].density();
//
//                 The resulting flux over the boundary is zero anyway (velocity is zero), but this will add some non-zero derivatives to the Jacobian and makes the BC more general.
//                 values[Indices::conti0EqIdx] = this->faceVelocity(element, fvGeometry, scvf) * insideDensity * scvf.unitOuterNormal();
//
//                 and now for the diffusion of co2
//                 static const Scalar dirichletMoleFraction = 2e-5;
//                 auto d = ipGlobal - element.geometry().center();
//                 d /= d.two_norm2();
//
//                 const auto& volVars = elemVolVars[scvf.insideScvIdx()];
//                 const auto tij = vtmv(scvf.unitOuterNormal(), volVars.diffusionCoefficient(0,0,1), d);
//                 auto molarDensity = volVars.molarDensity();
//                 const auto diffusiveFlux = -1.0*molarDensity*tij*(dirichletMoleFraction - volVars.moleFraction(0,1));
//                 values[Indices::conti0EqIdx+1] = diffusiveFlux;
//             }
//         }

        return values;
    }




   /*!
     * \brief Evaluates the boundary conditions for a Dirichlet control volume.
     *
     * \param globalPos The center of the finite volume which ought to be set.
     */
    DirichletValues dirichlet(const Element& element, const SubControlVolumeFace& scvf) const
    {
        DirichletValues values(0.0);
        return values;
    }

   /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluate the initial value for a control volume.
     *
     * \param globalPos The global position
     */
    InitialValues initialAtPos(const GlobalPosition &globalPos) const
    {
        InitialValues values(0.0);
        if constexpr (!ParentType::isMomentumProblem())
        {
            values[transportCompIdx] = 4e-7;
            values[Indices::pressureIdx] = 1.e+5+(this->gridGeometry().bBoxMax()[1] - globalPos[1])*999.7*9.81;
        }

        return values;
    }

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

        // the pure Neumann problem is only defined up to a constant
        // we create a well-posed problem by fixing the pressure at one dof

        if constexpr (!ParentType::isMomentumProblem())
        {
            const bool isLowerLeftCell = (scv.dofIndex() == 0);
            if (isLowerLeftCell)
                values.set(0); //only set true for 0 index (pressure)
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
        DirichletValues values(0.0);
        const auto globalPos = scv.center();
        values[0] = values[Indices::pressureIdx] = 1.e+5+(this->gridGeometry().bBoxMax()[1] - globalPos[1])*999.7*9.81;
        return values;
    }


private:

//     template<class SubControlVolume>
//     bool isUpperRow_(const SubControlVolume& scv) const
//     {
//         const auto globalPos = scv.center();
//
//         return (globalPos[1] > 5.995- eps_);
//     }

    bool onLeftBoundary_(const GlobalPosition& globalPos) const
    { return globalPos[0] < this->gridGeometry().bBoxMin()[0] + eps_; }

    bool onRightBoundary_(const GlobalPosition& globalPos) const
    { return globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps_; }

    bool onTopBoundary_(const GlobalPosition& globalPos) const
    { return globalPos[1] > this->gridGeometry().bBoxMax()[1] - eps_; }

    const Scalar eps_{1e-6};
};
} // end namespace Dumux

#endif
