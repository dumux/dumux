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
 * \ingroup NavierStokesTests
 * \brief Test for the staggered grid (Navier-)Stokes model with analytical solution.
 */
#ifndef DUMUX_TRIGONOMETRIC_3D_TEST_PROBLEM_HH
#define DUMUX_TRIGONOMETRIC_3D_TEST_PROBLEM_HH

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/freeflow/navierstokes/problem.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesTests
 * \brief  Test problem for the staggered grid.
 */
template <class TypeTag>
class TrignometricTestProblem : public NavierStokesProblem<TypeTag>
{
    using ParentType = NavierStokesProblem<TypeTag>;

    using BoundaryTypes = typename ParentType::BoundaryTypes;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using NumEqVector = typename ParentType::NumEqVector;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using PrimaryVariables = typename ParentType::PrimaryVariables;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;

    static constexpr auto dimWorld = GridGeometry::GridView::dimensionworld;
    using Element = typename FVElementGeometry::Element;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

public:
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    TrignometricTestProblem(std::shared_ptr<const GridGeometry> gridGeometry, std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry, couplingManager)
    {
        mu_ = getParam<Scalar>("Component.LiquidKinematicViscosity", 1.0);
    }

    TrignometricTestProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        mu_ = getParam<Scalar>("Component.LiquidKinematicViscosity", 1.0);
    }

    /*!
     * \name Problem parameters
     */
    // \{

    /*!
     * \brief Return the sources within the domain.
     *
     * \param globalPos The global position
     */
    NumEqVector sourceAtPos(const GlobalPosition &globalPos) const
    {
        if constexpr (ParentType::isMomentumProblem())
        {
            NumEqVector source;
            using std::sin;
            using std::cos;
            const Scalar x = globalPos[0];
            const Scalar y = globalPos[1];
            const Scalar z = globalPos[2];

            source[Indices::momentumXBalanceIdx] = sin(2*M_PI*x) * cos(2*M_PI*y) *cos(2*M_PI*z);
            source[Indices::momentumYBalanceIdx] = cos(2*M_PI*x) * sin(2*M_PI*y) *cos(2*M_PI*z);
            source[Indices::momentumZBalanceIdx] = -2*cos(2*M_PI*x) * cos(2*M_PI*y) *sin(2*M_PI*z);

            source *= 6*M_PI*M_PI;
            source += gradP_(globalPos);

            return source;
        }
        else
        {
            return NumEqVector(0.0);
        }
    }
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
    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;

        // set Dirichlet values for the velocity and pressure everywhere
        if constexpr (ParentType::isMomentumProblem())
        {
            values.setAllDirichlet();
        }
        else
            values.setNeumann(Indices::conti0EqIdx);

        return values;
    }

    /*!
     * \brief Return dirichlet boundary values at a given position
     *
     * \param globalPos The global position
     */
    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    {
        return analyticalSolution(globalPos);
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
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElementFluxVariablesCache& elemFluxVarsCache,
                        const SubControlVolumeFace& scvf) const
    {
        NumEqVector values(0.0);

        if constexpr (!ParentType::isMomentumProblem())
        {
            const auto insideDensity = elemVolVars[scvf.insideScvIdx()].density();
            values[Indices::conti0EqIdx] = this->faceVelocity(element, fvGeometry, scvf) * insideDensity * scvf.unitOuterNormal();
        }

        return values;
    }

    /*!
     * \brief Return the analytical solution of the problem at a given position
     *
     * \param globalPos The global position
     */
    PrimaryVariables analyticalSolution(const GlobalPosition& globalPos, Scalar time = 0.0) const
    {
        PrimaryVariables values(0.0);

        if constexpr (ParentType::isMomentumProblem())
        {
            using std::sin;
            using std::cos;
            const Scalar x = globalPos[0];
            const Scalar y = globalPos[1];
            const Scalar z = globalPos[2];
            values[Indices::velocityXIdx] = 0.5*sin(2*M_PI*x)*cos(2*M_PI*y)*cos(2*M_PI*z);
            values[Indices::velocityYIdx] = 0.5*cos(2*M_PI*x)*sin(2*M_PI*y)*cos(2*M_PI*z);
            values[Indices::velocityZIdx] = -cos(2*M_PI*x)*cos(2*M_PI*y)*sin(2*M_PI*z);;
        }
        else
            values[Indices::pressureIdx] = p_(globalPos);

        return values;
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    //! TODO should these be spatial params?
    Scalar pressureAtPos(const GlobalPosition& globalPos) const
    { return p_(globalPos); }

    Scalar densityAtPos(const GlobalPosition& globalPos) const
    { return 1.0; }

    Scalar effectiveViscosityAtPos(const GlobalPosition& globalPos) const
    { return 1.0; }

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
    std::bitset<PrimaryVariables::dimension> hasInternalDirichletConstraint(const Element& element, const SubControlVolume& scv) const
    {
        std::bitset<PrimaryVariables::dimension> values;

        auto fvGeometry = localView(this->gridGeometry());
        fvGeometry.bindElement(element);

        bool onBoundary = false;
        for (const auto& scvf : scvfs(fvGeometry))
            if(fvGeometry.scv(scvf.insideScvIdx()).dofIndex() == scv.dofIndex())
                onBoundary = std::max(onBoundary, scvf.boundary());

        if (onBoundary)
            values.set(0);

        return values;
    }

    /*!
     * \brief Define the values of internal Dirichlet constraints for a degree of freedom.
     * \param element The finite element
     * \param scv The sub-control volume
     */
    PrimaryVariables internalDirichlet(const Element& element, const SubControlVolume& scv) const
    { return PrimaryVariables(analyticalSolution(scv.dofPosition())[Indices::pressureIdx]); }

private:
    Scalar p_(const GlobalPosition& globalPos) const
    {
        using std::sin;
        const Scalar x = globalPos[0];
        const Scalar y = globalPos[1];
        const Scalar z = globalPos[2];
        return sin(2*M_PI*x)*sin(2*M_PI*y)*sin(2*M_PI*z);
    }

    GlobalPosition gradP_(const GlobalPosition& globalPos) const
    {
        using std::sin;
        using std::cos;
        const Scalar x = globalPos[0];
        const Scalar y = globalPos[1];
        const Scalar z = globalPos[2];
        return {2*M_PI*cos(2*M_PI*x)*sin(2*M_PI*y)*sin(2*M_PI*z),
                2*M_PI*sin(2*M_PI*x)*cos(2*M_PI*y)*sin(2*M_PI*z),
                2*M_PI*sin(2*M_PI*x)*sin(2*M_PI*y)*cos(2*M_PI*z)};
    }

    Scalar mu_;
};

} // end namespace Dumux

#endif
