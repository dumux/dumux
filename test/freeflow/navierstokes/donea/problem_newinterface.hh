// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup NavierStokesTests
 * \brief Test for the (hybrid) CVFE (Navier-)Stokes models with analytical solution (Donea 2003, \cite Donea2003).
 */
#ifndef DUMUX_DONEA_TEST_PROBLEM_NEWINTERFACE_HH
#define DUMUX_DONEA_TEST_PROBLEM_NEWINTERFACE_HH

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dumux/common/math.hh>
#include <dumux/geometry/diameter.hh>

#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/discretization/extrusion.hh>
#include <dumux/discretization/method.hh>
#include <dumux/discretization/cvfe/localdof.hh>
#include <dumux/discretization/dirichletconstraints.hh>
#include <dumux/common/constraintinfo.hh>

namespace Dumux {

/*!
 * \ingroup NavierStokesTests
 * \brief  Test problem for the (hybrid) CVFE schemes (Donea 2003, \cite Donea2003).
 *
 * A two-dimensional Stokes flow in a square domain is considered.
 * With the source terms as given in Donea 2003 \cite Donea2003, an analytical solution
 * is available and can be compared to the numerical solution.
 */
template <class TypeTag, class BaseProblem>
class DoneaTestProblemNewInterface : public BaseProblem
{
    using ParentType = BaseProblem;

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using ModelTraits = GetPropType<TypeTag, Properties::ModelTraits>;
    using Sources = typename ParentType::Sources;
    using DirichletValues = typename ParentType::DirichletValues;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using ConstraintInfo = Dumux::DirichletConstraintInfo<ModelTraits::numEq()>;
    using ConstraintValues = Dune::FieldVector<Scalar, ModelTraits::numEq()>;
    using GridIndexType = typename IndexTraits<typename GridGeometry::GridView>::GridIndex;
    using DirichletConstraintData = Dumux::DirichletConstraintData<ConstraintInfo, ConstraintValues, GridIndexType>;
    using BoundaryTypes = typename ParentType::BoundaryTypes;
    using BoundaryFluxes = typename ParentType::BoundaryFluxes;

    static constexpr auto dimWorld = GridGeometry::GridView::dimensionworld;
    using Element = typename FVElementGeometry::Element;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

public:
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    DoneaTestProblemNewInterface(std::shared_ptr<const GridGeometry> gridGeometry, std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry, couplingManager)
    {
        useNeumann_ = getParam<bool>("Problem.UseNeumann", false);
        mu_ = getParam<Scalar>("Component.LiquidKinematicViscosity", 1.0);

        if constexpr (ParentType::isMomentumProblem())
        {
            CVFE::appendDirichletConstraints(*this,
                [&, this](const auto& fvGeometry, const auto&, const auto& localDof){
                    return this->dirichletAtPos(ipData(fvGeometry, localDof).global());
                },
                constraints_
            );
        }
        else if(!useNeumann_)
            appendInternalConstraints_();
    }

    DoneaTestProblemNewInterface(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        useNeumann_ = getParam<bool>("Problem.UseNeumann", false);
        mu_ = getParam<Scalar>("Component.LiquidKinematicViscosity", 1.0);

        CVFE::appendDirichletConstraints(*this,
            [&, this](const auto& fvGeometry, const auto&, const auto& localDof){
                return this->dirichletAtPos(ipData(fvGeometry, localDof).global());
            },
            constraints_
        );
    }

    /*!
     * \brief Return the sources within the domain.
     *
     * \param globalPos The global position
     */
    Sources sourceAtPos(const GlobalPosition &globalPos) const
    {
        if constexpr (ParentType::isMomentumProblem())
        {
            Sources source(0.0);
            const Scalar x = globalPos[0];
            const Scalar y = globalPos[1];

            source[Indices::momentumXBalanceIdx] = -2.0*mu_*dxxU_(x,y) - mu_*dyyU_(x,y) - mu_*dxyV_(x,y) + dxP_(x,y);
            source[Indices::momentumYBalanceIdx] = -2.0*mu_*dyyV_(x,y) - mu_*dxyU_(x,y) - mu_*dxxV_(x,y) + dyP_(x,y);
            return source;
        }
        else
        {
            return Sources(0.0);
        }
    }

    /*!
     * \name Boundary conditions
     */
    // \{

    /*!
     * \brief Specifies which kind of boundary condition should be
     *        used for which equation on the boundary.
     *
     * \param globalPos The position at the boundary
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
    DirichletValues dirichletAtPos(const GlobalPosition& globalPos) const
    { return analyticalSolution(globalPos); }

    /*!
     * \brief Return  Dirichlet boundary constraints and internal constraints.
     */
    const auto& constraints() const
    { return constraints_; }

    /*!
     * \brief Evaluates the boundary flux related to a localDof at a given interpolation point.
     *
     * \param fvGeometry The finite-volume geometry
     * \param elemVars All variables related to the element
     * \param elemFluxVarsCache The element flux variables cache
     * \param faceIpData Interpolation point data
     */
    template<class ElementVariables, class ElementFluxVariablesCache, class FaceIpData>
    BoundaryFluxes boundaryFlux(const FVElementGeometry& fvGeometry,
                                const ElementVariables& elemVars,
                                const ElementFluxVariablesCache& elemFluxVarsCache,
                                const FaceIpData& faceIpData) const
    {
        BoundaryFluxes values(0.0);

        if constexpr (ParentType::isMomentumProblem())
        {
            const auto x = faceIpData.global()[0];
            const auto y = faceIpData.global()[1];

            Dune::FieldMatrix<Scalar, dimWorld, dimWorld> momentumFlux(0.0);
            momentumFlux[0][0] = -2.0*mu_*dxU_(x,y) + p_(x);
            momentumFlux[0][1] = -mu_*dyU_(x,y) - mu_*dxV_(x,y);
            momentumFlux[1][0] = momentumFlux[0][1];
            momentumFlux[1][1] = -2.0*mu_*dyV_(x,y) + p_(x);

            const auto normal = faceIpData.unitOuterNormal();
            momentumFlux.mv(normal, values);
        }
        else
        {
            const auto& scvf = fvGeometry.scvf(faceIpData.scvfIndex());
            const auto insideDensity = elemVars[scvf.insideScvIdx()].density();
            values[Indices::conti0EqIdx] = this->velocity(fvGeometry, faceIpData) * insideDensity * scvf.unitOuterNormal();
        }

        return values;
    }

    // \}

    /*!
     * \brief Return the analytical solution of the problem at a given position
     *
     * \param globalPos The global position
     */
    DirichletValues analyticalSolution(const GlobalPosition& globalPos, Scalar time = 0.0) const
    {
        DirichletValues values(0.0);

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

    /*!
     * \brief Return the gradient of the analytical solution of the problem at a given position
     *
     * \param globalPos The global position
     */
    Dune::FieldVector<GlobalPosition, DirichletValues::dimension> gradAnalyticalSolution(const GlobalPosition& globalPos, Scalar time = 0.0) const
    {
        const Scalar x = globalPos[0];
        const Scalar y = globalPos[1];
        Dune::FieldVector<GlobalPosition, DirichletValues::dimension> values;

        if constexpr (ParentType::isMomentumProblem())
        {
            values[Indices::velocityXIdx][0] = dxU_(x,y);
            values[Indices::velocityXIdx][1] = dyU_(x,y);
            values[Indices::velocityYIdx][0] = dxV_(x,y);
            values[Indices::velocityYIdx][1] = dyV_(x,y);
        }
        else
        {
            values[Indices::pressureIdx][0] = dxP_(x,y);
            values[Indices::pressureIdx][1] = dyP_(x,y);
        }

        return values;
    }

    //! TODO should these be spatial params?
    Scalar pressureAtPos(const GlobalPosition& globalPos) const
    { return p_(globalPos[0]); }

    Scalar densityAtPos(const GlobalPosition& globalPos) const
    { return 1.0; }

    Scalar effectiveViscosityAtPos(const GlobalPosition& globalPos) const
    { return 1.0; }

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

    void appendInternalConstraints_()
    {
        static_assert(GridGeometry::discMethod == DiscretizationMethods::box, "Internal Dirichlet constraints only implemented for Box mass discretization scheme.");

        static constexpr Scalar eps = 1e-8;
        auto fvGeometry = localView(this->gridGeometry());
        for (const auto& element : elements(this->gridGeometry().gridView()))
        {
            fvGeometry.bind(element);
            for (const auto& scv : scvs(fvGeometry))
            {
                if  ((scv.dofPosition() - this->gridGeometry().bBoxMin()).two_norm() < eps)
                {
                    ConstraintInfo info;
                    info.setAll();

                    DirichletValues dirichletValues(analyticalSolution(scv.dofPosition())[Indices::pressureIdx]);
                    constraints_.push_back(DirichletConstraintData{std::move(info), std::move(dirichletValues), scv.dofIndex()});
                }
            }
        }
    }

    bool useNeumann_;
    Scalar mu_;
    std::vector<DirichletConstraintData> constraints_;
};

} // end namespace Dumux

#endif
