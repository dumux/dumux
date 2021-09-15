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
 * \brief Test for the staggered grid (Navier-)Stokes model with analytical solution (Donea 2003, \cite Donea2003).
 */
#ifndef DUMUX_DONEA_TEST_PROBLEM_HH
#define DUMUX_DONEA_TEST_PROBLEM_HH

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/freeflow/navierstokes/problem.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesTests
 * \brief  Test problem for the staggered grid (Donea 2003, \cite Donea2003).
 *
 * A two-dimensional Stokes flow in a square domain is considered.
 * With the source terms as given in Donea 2003 \cite Donea2003, an analytical solution
 * is available and can be compared to the numerical solution.
 */
template <class TypeTag>
class DoneaTestProblem : public NavierStokesProblem<TypeTag>
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

    DoneaTestProblem(std::shared_ptr<const GridGeometry> gridGeometry, std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry, couplingManager)
    {
        useNeumann_ = getParam<bool>("Problem.UseNeumann", false);
        mu_ = getParam<Scalar>("Component.LiquidKinematicViscosity", 1.0);
    }

    DoneaTestProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        useNeumann_ = getParam<bool>("Problem.UseNeumann", false);
        mu_ = getParam<Scalar>("Component.LiquidKinematicViscosity", 1.0);
    }

    /*!
     * \name Problem parameters
     */
    // \{

    Scalar temperature() const { return 298.0; }

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
            const Scalar x = globalPos[0];
            const Scalar y = globalPos[1];

            source[Indices::momentumXBalanceIdx] = -2.0*mu_*dxxU_(x,y) - mu_*dyyU_(x,y) - mu_*dxyV_(x,y) + dxP_(x,y);
            source[Indices::momentumYBalanceIdx] = -2.0*mu_*dyyV_(x,y) - mu_*dxyU_(x,y) - mu_*dxxV_(x,y) + dyP_(x,y);
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
            if (useNeumann_)
            {
                static constexpr Scalar eps = 1e-8;
                if ((globalPos[0] > this->gridGeometry().bBoxMax()[0] - eps) || (globalPos[1] > this->gridGeometry().bBoxMax()[1] - eps))
                    values.setAllNeumann();
                else
                    values.setAllDirichlet();
            }
            else
            {
                values.setDirichlet(Indices::velocityXIdx);
                values.setDirichlet(Indices::velocityYIdx);
            }
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

        if constexpr (ParentType::isMomentumProblem())
        {
            const auto x = scvf.ipGlobal()[0];
            const auto y = scvf.ipGlobal()[1];

            Dune::FieldMatrix<Scalar, dimWorld, dimWorld> momentumFlux(0.0);
            momentumFlux[0][0] = -2.0*mu_*dxU_(x,y) + p_(x);
            momentumFlux[0][1] = -mu_*dyU_(x,y) - mu_*dxV_(x,y);
            momentumFlux[1][0] = momentumFlux[0][1];
            momentumFlux[1][1] = -2.0*mu_*dyV_(x,y) + p_(x);

            const auto normal = scvf.unitOuterNormal();
            momentumFlux.mv(normal, values);
        }
        else
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
        PrimaryVariables values;

        if constexpr (ParentType::isMomentumProblem())
        {
            const Scalar x = globalPos[0];
            const Scalar y = globalPos[1];
            values[Indices::velocityXIdx] = f2_(x)*df2_(y);
            values[Indices::velocityYIdx] = -f2_(y)*df2_(x);
        }
        else
            values[Indices::pressureIdx] = p_(globalPos[0]);

        return values;
    }

    // \}

    /*!
     * \name Volume terms
     */
    // \{

    //! TODO should these be spatial params?
    Scalar pressureAtPos(const GlobalPosition& globalPos) const
    { return p_(globalPos[0]); }

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

        if (!useNeumann_)
        {
            auto fvGeometry = localView(this->gridGeometry());
            fvGeometry.bindElement(element);

            bool onBoundary = false;
            for (const auto& scvf : scvfs(fvGeometry))
                onBoundary = std::max(onBoundary, scvf.boundary());

            if (onBoundary)
                values.set(0);

            // TODO: only use one cell or pass fvGeometry to hasInternalDirichletConstraint

            // if (scv.dofIndex() == 0)
            //     values.set(0);
            // the pure Neumann problem is only defined up to a constant
            // we create a well-posed problem by fixing the pressure at one dof
        }

        return values;
    }

    /*!
     * \brief Define the values of internal Dirichlet constraints for a degree of freedom.
     * \param element The finite element
     * \param scv The sub-control volume
     */
    PrimaryVariables internalDirichlet(const Element& element, const SubControlVolume& scv) const
    { return PrimaryVariables(analyticalSolution(scv.center())[Indices::pressureIdx]); }

private:
    Scalar p_(Scalar x) const
    { return x*(1.0-x); }

    Scalar dP_(Scalar x) const
    { return 1.0 - 2.0*x; }

    Scalar f2_(Scalar x) const
    { return p_(x)*p_(x); /*=x^2*(1-2x+x^2)=x^2-2x^3+x^4*/ }

    Scalar df2_(Scalar x) const
    { return 2.0*x - 6.0*x*x + 4.0*x*x*x; }

    Scalar ddf2_(Scalar x) const
    { return 2.0 - 12.0*x + 12.0*x*x; }

    Scalar dddf2_(Scalar x) const
    { return - 12.0 + 24.0*x; }

    Scalar dxP_ (Scalar x, Scalar y) const
    { return dP_(x); }

    Scalar dyP_ (Scalar x, Scalar y) const
    { return 0.0; }

    Scalar dxU_ (Scalar x, Scalar y) const
    { return df2_(x)*df2_(y); }

    Scalar dyU_ (Scalar x, Scalar y) const
    { return f2_(x)*ddf2_(y); }

    Scalar dxV_ (Scalar x, Scalar y) const
    { return -f2_(y)*ddf2_(x); }

     Scalar dyV_ (Scalar x, Scalar y) const
    { return -df2_(y)*df2_(x); }

    Scalar dxxU_ (Scalar x, Scalar y) const
    { return ddf2_(x)*df2_(y); }

    Scalar dxyU_ (Scalar x, Scalar y) const
    { return df2_(x)*ddf2_(y); }

    Scalar dyyU_ (Scalar x, Scalar y) const
    { return f2_(x)*dddf2_(y); }

    Scalar dyyV_ (Scalar x, Scalar y) const
    { return -ddf2_(y)*df2_(x); }

    Scalar dxyV_ (Scalar x, Scalar y) const
    { return -df2_(y)*ddf2_(x); }

    Scalar dxxV_ (Scalar x, Scalar y) const
    { return -f2_(y)*dddf2_(x); }

    bool useNeumann_;
    Scalar mu_;
};

} // end namespace Dumux

#endif
