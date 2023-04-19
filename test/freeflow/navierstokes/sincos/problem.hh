// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesTests
 * \brief Test for the staggered grid Navier-Stokes model with analytical solution.
 */
#ifndef DUMUX_SINCOS_TEST_PROBLEM_HH
#define DUMUX_SINCOS_TEST_PROBLEM_HH

#include <dune/common/fmatrix.hh>

#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>

#include <dumux/freeflow/navierstokes/boundarytypes.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesTests
 * \brief  Test problem for the staggered grid.
 *
 * The 2D, incompressible Navier-Stokes equations for zero gravity and a Newtonian
 * flow is solved and compared to an analytical solution (sums/products of trigonometric functions).
 * For the instationary case, the velocities and pressures are periodical in time. The Dirichlet boundary conditions are
 * consistent with the analytical solution and in the instationary case time-dependent.
 */
template <class TypeTag, class BaseProblem>
class SincosTestProblem : public BaseProblem
{
    using ParentType = BaseProblem;

    using BoundaryTypes = typename ParentType::BoundaryTypes;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using InitialValues = typename ParentType::InitialValues;
    using Sources = typename ParentType::Sources;
    using DirichletValues = typename ParentType::DirichletValues;
    using BoundaryFluxes = typename ParentType::BoundaryFluxes;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using FVElementGeometry = typename GridGeometry::LocalView;

    static constexpr auto dimWorld = GridGeometry::GridView::dimensionworld;
    using Element = typename FVElementGeometry::Element;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

public:
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    SincosTestProblem(std::shared_ptr<const GridGeometry> gridGeometry, std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry, couplingManager), time_(0.0)
    {
        isStationary_ = getParam<bool>("Problem.IsStationary");
        rho_ = getParam<Scalar>("Component.LiquidDensity");
        Scalar nu = getParam<Scalar>("Component.LiquidKinematicViscosity", 1.0);
        mu_ = rho_*nu; // dynamic viscosity
        useNeumann_ = getParam<bool>("Problem.UseNeumann", false);
    }

    /*!
     * \brief Return the sources within the domain.
     *
     * \param globalPos The global position
     */
    Sources sourceAtPos(const GlobalPosition &globalPos) const
    {
        Sources source(0.0);
        if constexpr (ParentType::isMomentumProblem())
        {
            const Scalar x = globalPos[0];
            const Scalar y = globalPos[1];
            const Scalar t = time_;

            source[Indices::momentumXBalanceIdx] = rho_*dtU_(x,y,t) - 2.0*mu_*dxxU_(x,y,t) - mu_*dyyU_(x,y,t) - mu_*dxyV_(x,y,t) + dxP_(x,y,t);
            source[Indices::momentumYBalanceIdx] = rho_*dtV_(x,y,t) - 2.0*mu_*dyyV_(x,y,t) - mu_*dxyU_(x,y,t) - mu_*dxxV_(x,y,t) + dyP_(x,y,t);

            if (this->enableInertiaTerms())
            {
                source[Indices::momentumXBalanceIdx] += rho_*dxUU_(x,y,t) + rho_*dyUV_(x,y,t);
                source[Indices::momentumYBalanceIdx] += rho_*dxUV_(x,y,t) + rho_*dyVV_(x,y,t);
            }
        }

        return source;
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

        if constexpr (ParentType::isMomentumProblem())
        {
            if (useNeumann_)
            {
                if (globalPos[1] > this->gridGeometry().bBoxMax()[1] - 1e-6
                || globalPos[0] > this->gridGeometry().bBoxMax()[0] - 1e-6
                || globalPos[1] < this->gridGeometry().bBoxMin()[1] + 1e-6)
                {
                    values.setNeumann(Indices::velocityXIdx);
                    values.setNeumann(Indices::velocityYIdx);
                }
                else
                {
                    values.setDirichlet(Indices::velocityXIdx);
                    values.setDirichlet(Indices::velocityYIdx);
                }
            }
            else
            {
                // set Dirichlet values for the velocity everywhere
                values.setAllDirichlet();
            }
        }
        else
            values.setAllNeumann();

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
        // set fixed pressure in one cell
        std::bitset<DirichletValues::dimension> values;

        if (!useNeumann_ && scv.dofIndex() == 0)
            values.set(Indices::pressureIdx);

        return values;
    }

    /*!
     * \brief Define the values of internal Dirichlet constraints for a degree of freedom.
     * \param element The finite element
     * \param scv The sub-control volume
     */
    DirichletValues internalDirichlet(const Element& element, const SubControlVolume& scv) const
    { return DirichletValues(analyticalSolution(scv.center())[Indices::pressureIdx]); }

    /*!
     * \brief Returns Dirichlet boundary values at a given position.
     *
     * \param globalPos The global position
     */
    DirichletValues dirichletAtPos(const GlobalPosition& globalPos) const
    {
        // use the values of the analytical solution
        return analyticalSolution(globalPos, time_);
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

        if constexpr (ParentType::isMomentumProblem())
        {
            const auto flux = [&](Scalar x, Scalar y)
            {
                Dune::FieldMatrix<Scalar, dimWorld, dimWorld> momentumFlux(0.0);
                const Scalar t = time_;

                momentumFlux[0][0] = -2*mu_*dxU_(x,y,t) + p_(x,y,t);
                momentumFlux[0][1] = -mu_*(dyU_(x,y,t) + dxV_(x,y,t));
                momentumFlux[1][0] = -mu_*(dyU_(x,y,t) + dxV_(x,y,t));
                momentumFlux[1][1] = -2*mu_*dyV_(x,y,t) + p_(x,y,t);

                if (this->enableInertiaTerms())
                {
                    momentumFlux[0][0] += rho_*u_(x,y,t)*u_(x,y,t);
                    momentumFlux[0][1] += rho_*u_(x,y,t)*v_(x,y,t);
                    momentumFlux[1][0] += rho_*v_(x,y,t)*u_(x,y,t);
                    momentumFlux[1][1] += rho_*v_(x,y,t)*v_(x,y,t);
                }

                return momentumFlux;
            };

            flux(scvf.ipGlobal()[0], scvf.ipGlobal()[1]).mv(scvf.unitOuterNormal(), values);
        }
        else
        {
            const auto insideDensity = elemVolVars[scvf.insideScvIdx()].density();
            values[Indices::conti0EqIdx] = insideDensity * (this->faceVelocity(element, fvGeometry, scvf) * scvf.unitOuterNormal());
        }

        return values;
    }

    /*!
     * \brief Returns the analytical solution of the problem at a given time and position.
     *
     * \param globalPos The global position
     * \param time The current simulation time
     */
    DirichletValues analyticalSolution(const GlobalPosition& globalPos, const Scalar time) const
    {
        const Scalar x = globalPos[0];
        const Scalar y = globalPos[1];
        const Scalar t = time;
        DirichletValues values;

        if constexpr (ParentType::isMomentumProblem())
        {
            values[Indices::velocityXIdx] = u_(x,y,t);
            values[Indices::velocityYIdx] = v_(x,y,t);
        }
        else
            values[Indices::pressureIdx] = p_(x,y,t);

        return values;
    }

    // \}

   /*!
     * \name Volume terms
     */
    // \{

    /*!
     * \brief Evaluates the initial value for a control volume.
     *
     * \param globalPos The global position
     */
    InitialValues initialAtPos(const GlobalPosition &globalPos) const
    {
        if (isStationary_)
            return InitialValues(0.0);
        else
            return analyticalSolution(globalPos, 0.0);
    }


    /*!
     * \brief Returns the analytical solution of the problem at a given position.
     *
     * \param globalPos The global position
     */
    DirichletValues analyticalSolution(const GlobalPosition& globalPos) const
    {
        return analyticalSolution(globalPos, time_);
    }

    /*!
     * \brief Updates the time
     */
    void updateTime(const Scalar time)
    {
        time_ = time;
    }

private:
    Scalar f_(Scalar t) const
    {
        using std::sin;
        if (isStationary_)
            return 1.0;
        else
            return sin(2.0 * t);
    }

    Scalar df_(Scalar t) const
    {
        using std::cos;
        if (isStationary_)
            return 0.0;
        else
            return 2.0 * cos(2.0 * t);
    }

    Scalar f1_(Scalar x) const
    { using std::cos; return -0.25 * cos(2.0 * x); }

    Scalar df1_(Scalar x) const
    { using std::sin; return 0.5 * sin(2.0 * x); }

    Scalar f2_(Scalar x) const
    { using std::cos; return -cos(x); }

    Scalar df2_(Scalar x) const
    { using std::sin; return sin(x); }

    Scalar ddf2_(Scalar x) const
    { using std::cos; return cos(x); }

    Scalar dddf2_(Scalar x) const
    { using std::sin; return -sin(x); }

    Scalar p_(Scalar x, Scalar y, Scalar t) const
    { return (f1_(x) + f1_(y)) * f_(t) * f_(t); }

    Scalar dxP_ (Scalar x, Scalar y, Scalar t) const
    { return df1_(x) * f_(t) * f_(t); }

    Scalar dyP_ (Scalar x, Scalar y, Scalar t) const
    { return df1_(y) * f_(t) * f_(t); }

    Scalar u_(Scalar x, Scalar y, Scalar t) const
    { return f2_(x)*df2_(y) * f_(t); }

    Scalar dtU_ (Scalar x, Scalar y, Scalar t) const
    { return f2_(x)*df2_(y) * df_(t); }

    Scalar dxU_ (Scalar x, Scalar y, Scalar t) const
    { return df2_(x)*df2_(y) * f_(t); }

    Scalar dyU_ (Scalar x, Scalar y, Scalar t) const
    { return f2_(x)*ddf2_(y) * f_(t); }

    Scalar dxxU_ (Scalar x, Scalar y, Scalar t) const
    { return ddf2_(x)*df2_(y) * f_(t); }

    Scalar dxyU_ (Scalar x, Scalar y, Scalar t) const
    { return df2_(x)*ddf2_(y) * f_(t); }

    Scalar dyyU_ (Scalar x, Scalar y, Scalar t) const
    { return f2_(x)*dddf2_(y) * f_(t); }

    Scalar v_(Scalar x, Scalar y, Scalar t) const
    { return -f2_(y)*df2_(x) * f_(t); }

    Scalar dtV_ (Scalar x, Scalar y, Scalar t) const
    { return -f2_(y)*df2_(x) * df_(t); }

    Scalar dxV_ (Scalar x, Scalar y, Scalar t) const
    { return -f2_(y)*ddf2_(x) * f_(t); }

    Scalar dyV_ (Scalar x, Scalar y, Scalar t) const
    { return -df2_(y)*df2_(x) * f_(t); }

    Scalar dyyV_ (Scalar x, Scalar y, Scalar t) const
    { return -ddf2_(y)*df2_(x) * f_(t); }

    Scalar dxyV_ (Scalar x, Scalar y, Scalar t) const
    { return -df2_(y)*ddf2_(x) * f_(t); }

    Scalar dxxV_ (Scalar x, Scalar y, Scalar t) const
    { return -f2_(y)*dddf2_(x) * f_(t); }

    Scalar dxUU_ (Scalar x, Scalar y, Scalar t) const
    { return 2.*u_(x,y,t)*dxU_(x,y,t); }

    Scalar dyVV_ (Scalar x, Scalar y, Scalar t) const
    { return 2.*v_(x,y,t)*dyV_(x,y,t); }

    Scalar dxUV_ (Scalar x, Scalar y, Scalar t) const
    { return v_(x,y,t)*dxU_(x,y,t) + u_(x,y,t)*dxV_(x,y,t); }

    Scalar dyUV_ (Scalar x, Scalar y, Scalar t) const
    { return v_(x,y,t)*dyU_(x,y,t) + u_(x,y,t)*dyV_(x,y,t); }

    Scalar rho_;
    Scalar mu_;
    Scalar time_;

    bool isStationary_;
    bool useNeumann_;
};

} // end namespace Dumux

#endif
