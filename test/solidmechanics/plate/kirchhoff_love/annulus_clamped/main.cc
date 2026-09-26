// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Annular Kirchhoff-Love plate clamped at both edges: the gauge of the shear
 *        gradient potential on a fully clamped hole boundary is an unknown constant
 */
#include <config.h>
#include <iostream>
#include <cmath>
#include <bitset>
#include <type_traits>

#include <dune/common/fmatrix.hh>
#include <dune/common/fvector.hh>
#include <dune/foamgrid/foamgrid.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/common/fvproblem.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/math.hh>
#include <dumux/discretization/evalsolution.hh>
#include <dumux/discretization/elementsolution.hh>
#include <dumux/discretization/pq1bubble.hh>
#include <dumux/discretization/box.hh>
#include <dumux/geometry/diameter.hh>

#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/istlsolvers.hh>

#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/newtonsolver.hh>
#include <dumux/multidomain/fvassembler.hh>

#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/grid/gridmanager_foam.hh>

#include <dumux/solidmechanics/plate/kirchhoff_love/model.hh>
#include <dumux/solidmechanics/plate/kirchhoff_love/couplingmanager.hh>

namespace Dumux {

/*!
 * \brief Axisymmetric deflection of an annular plate clamped at both edges under a uniform load
 *
 * w = -F r^4/(64 D) + c1 + c2 r^2 + c3 ln r + c4 r^2 ln r with w = w' = 0 at r = a and r = b.
 */
template<class Scalar>
struct ClampedClampedAnnulus
{
    ClampedClampedAnnulus(Scalar a, Scalar b, Scalar F, Scalar D)
    : a_(a), b_(b), F_(F), D_(D)
    {
        using std::log;
        Dune::FieldMatrix<Scalar, 4, 4> A;
        Dune::FieldVector<Scalar, 4> rhs;
        for (int i = 0; i < 2; ++i)
        {
            const auto r = (i == 0) ? a : b;
            A[2*i] = {1.0, r*r, log(r), r*r*log(r)};
            A[2*i + 1] = {0.0, 2.0*r, 1.0/r, 2.0*r*log(r) + r};
            rhs[2*i] = F*r*r*r*r/(64.0*D);
            rhs[2*i + 1] = F*r*r*r/(16.0*D);
        }
        A.solve(c_, rhs);
    }

    Scalar deflection(Scalar r) const
    {
        using std::log;
        return -F_*r*r*r*r/(64.0*D_) + c_[0] + c_[1]*r*r + c_[2]*log(r) + c_[3]*r*r*log(r);
    }

    Scalar laplacian(Scalar r) const
    {
        using std::log;
        return -F_*r*r/(4.0*D_) + 4.0*c_[1] + c_[3]*(4.0*log(r) + 4.0);
    }

    //! phi = -D (Delta w - Delta w(b)) is the gauge with phi = 0 at the outer edge
    Scalar shearGradPotential(Scalar r) const
    { return -D_*(laplacian(r) - laplacian(b_)); }

    Scalar innerPotential() const
    { return shearGradPotential(a_); }

private:
    Scalar a_, b_, F_, D_;
    Dune::FieldVector<Scalar, 4> c_;
};

template<class TypeTag>
class RotationProblem : public FVProblem<TypeTag>
{
    using ParentType = FVProblem<TypeTag>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

public:
    RotationProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                    std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry)
    , couplingManager_(couplingManager)
    , poissonRatio_(getParam<Scalar>("Problem.PoissonRatio"))
    {
        const auto E = getParam<Scalar>("Problem.E");
        const auto t = getParam<Scalar>("Problem.Thickness");
        stiffness_ = E*t*t*t/(12.0*(1.0 - poissonRatio_*poissonRatio_));
    }

    Scalar D(const GlobalPosition& globalPos) const { return stiffness_; }
    Scalar poissonRatio(const GlobalPosition& globalPos) const { return poissonRatio_; }

    BoundaryTypes boundaryTypes(const Element& element, const SubControlVolume& scv) const
    {
        BoundaryTypes values;
        values.setAllDirichlet();
        return values;
    }

    template<class ElementVolumeVariables, class ElementFluxVariablesCache>
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElementFluxVariablesCache& elemFluxVarsCache,
                        const SubControlVolumeFace& scvf) const
    { return NumEqVector(0.0); }

    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(0.0); }

    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(0.0); }

    const CouplingManager& couplingManager() const
    { return *couplingManager_; }

private:
    std::shared_ptr<CouplingManager> couplingManager_;
    Scalar poissonRatio_, stiffness_;
};

/*!
 * \brief Deformation sub-problem: w = 0 on both edges, phi = 0 on the outer edge and
 *        phi = innerPotential on the inner edge, psi natural with one interior pin.
 *
 * With naturalInner set, the inner edge keeps w = 0 in place of the transverse balance and
 * assembles the compatibility equation naturally, so that the residual of that equation at
 * the inner nodes is the box flux of grad w - theta out of the ring of inner boxes.
 */
template<class TypeTag>
class DeformationProblem : public FVProblem<TypeTag>
{
    using ParentType = FVProblem<TypeTag>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using CouplingManager = GetPropType<TypeTag, Properties::CouplingManager>;

public:
    DeformationProblem(std::shared_ptr<const GridGeometry> gridGeometry,
                       std::shared_ptr<CouplingManager> couplingManager)
    : ParentType(gridGeometry)
    , couplingManager_(couplingManager)
    , force_(getParam<Scalar>("Problem.Force"))
    , innerRadius_(getParam<Scalar>("Problem.InnerRadius"))
    , outerRadius_(getParam<Scalar>("Problem.OuterRadius"))
    {
        const auto E = getParam<Scalar>("Problem.E");
        const auto t = getParam<Scalar>("Problem.Thickness");
        const auto nu = getParam<Scalar>("Problem.PoissonRatio");
        stiffness_ = E*t*t*t/(12.0*(1.0 - nu*nu));
    }

    Scalar stiffness() const { return stiffness_; }
    void setInnerPotential(Scalar value) { innerPotential_ = value; }
    void setNaturalInner(bool value) { naturalInner_ = value; }
    bool onInnerEdge(const GlobalPosition& globalPos) const
    { return globalPos.two_norm() < 0.5*(innerRadius_ + outerRadius_); }

    BoundaryTypes boundaryTypes(const Element& element, const SubControlVolume& scv) const
    {
        BoundaryTypes values;
        if (naturalInner_ && onInnerEdge(scv.dofPosition()))
        {
            values.setDirichlet(Indices::verticalDeformationIdx, Indices::shearGradPotentialEqIdx);
            values.setNeumann(Indices::deformationEqIdx);
            values.setNeumann(Indices::shearCurlPotentialEqIdx);
        }
        else
        {
            values.setDirichlet(Indices::shearGradPotentialIdx);
            values.setDirichlet(Indices::verticalDeformationIdx);
            values.setNeumann(Indices::shearCurlPotentialEqIdx);
        }
        return values;
    }

    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    {
        PrimaryVariables values(0.0);
        if (onInnerEdge(globalPos))
            values[Indices::shearGradPotentialIdx] = innerPotential_;
        return values;
    }

    template<class ElementVolumeVariables, class ElementFluxVariablesCache>
    NumEqVector neumann(const Element& element,
                        const FVElementGeometry& fvGeometry,
                        const ElementVolumeVariables& elemVolVars,
                        const ElementFluxVariablesCache& elemFluxVarsCache,
                        const SubControlVolumeFace& scvf) const
    {
        NumEqVector values(0.0);
        const auto rotation = this->couplingManager().rotation(fvGeometry, scvf);
        const auto tangent = [&](){
            auto tangent = scvf.unitOuterNormal();
            std::swap(tangent[0], tangent[1]);
            tangent[1] = -tangent[1];
            return tangent;
        }();
        values[Indices::shearCurlPotentialEqIdx] = -vtmv(tangent, 1.0, rotation);
        return values;
    }

    NumEqVector sourceAtPos(const GlobalPosition& globalPos) const
    {
        NumEqVector source(0.0);
        source[Indices::shearGradPotentialEqIdx] = force_;
        return source;
    }

    static constexpr bool enableInternalDirichletConstraints()
    { return true; }

    std::bitset<NumEqVector::dimension>
    hasInternalDirichletConstraint(const Element& element, const SubControlVolume& scv) const
    {
        std::bitset<NumEqVector::dimension> values;
        if (scv.dofIndex() == pinDofIndex_)
            values.set(Indices::shearCurlPotentialIdx);
        return values;
    }

    PrimaryVariables internalDirichlet(const Element& element, const SubControlVolume& scv) const
    { return PrimaryVariables(0.0); }

    void setPinDofIndex(std::size_t idx) { pinDofIndex_ = idx; }

    PrimaryVariables initialAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(0.0); }

    const CouplingManager& couplingManager() const
    { return *couplingManager_; }

private:
    std::shared_ptr<CouplingManager> couplingManager_;
    Scalar force_, innerRadius_, outerRadius_, stiffness_;
    Scalar innerPotential_ = 0.0;
    bool naturalInner_ = false;
    std::size_t pinDofIndex_ = 0;
};

} // end namespace Dumux

namespace Dumux::Properties {

namespace TTag {
struct AnnulusClampedCommon
{
    using Grid = Dune::FoamGrid<2, 2>;
    using EnableGridVolumeVariablesCache = std::true_type;
    using EnableGridFluxVariablesCache = std::true_type;
    using EnableGridGeometryCache = std::true_type;
};
struct AnnulusClampedRotation
{ using InheritsFrom = std::tuple<AnnulusClampedCommon, KirchhoffLovePlateRotation, PQ1BubbleModel>; };
struct AnnulusClampedDeformation
{ using InheritsFrom = std::tuple<AnnulusClampedCommon, KirchhoffLovePlateDeformation, BoxModel>; };
} // end namespace TTag

template<class TypeTag>
struct Problem<TypeTag, TTag::AnnulusClampedRotation>
{ using type = RotationProblem<TypeTag>; };
template<class TypeTag>
struct Problem<TypeTag, TTag::AnnulusClampedDeformation>
{ using type = DeformationProblem<TypeTag>; };

template<class TypeTag>
struct CouplingManager<TypeTag, TTag::AnnulusClampedRotation>
{
    using MDTraits = MultiDomainTraits<TTag::AnnulusClampedRotation, TTag::AnnulusClampedDeformation>;
    using type = KirchhoffLovePlateCouplingManager<MDTraits>;
};
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::AnnulusClampedDeformation>
{
    using MDTraits = MultiDomainTraits<TTag::AnnulusClampedRotation, TTag::AnnulusClampedDeformation>;
    using type = KirchhoffLovePlateCouplingManager<MDTraits>;
};

} // end namespace Dumux::Properties

int main(int argc, char** argv)
{
    using namespace Dumux;

    Dumux::initialize(argc, argv);
    Parameters::init(argc, argv);

    using RotationTypeTag = Properties::TTag::AnnulusClampedRotation;
    using DeformationTypeTag = Properties::TTag::AnnulusClampedDeformation;
    using CommonTypeTag = Properties::TTag::AnnulusClampedCommon;

    using Grid = GetPropType<CommonTypeTag, Properties::Grid>;
    GridManager<Grid> gridManager;
    gridManager.init();
    const auto& leafGridView = gridManager.grid().leafGridView();

    using RotationGridGeometry = GetPropType<RotationTypeTag, Properties::GridGeometry>;
    using DeformationGridGeometry = GetPropType<DeformationTypeTag, Properties::GridGeometry>;
    auto rotationGridGeometry = std::make_shared<RotationGridGeometry>(leafGridView);
    auto deformationGridGeometry = std::make_shared<DeformationGridGeometry>(leafGridView);

    using CouplingManager = GetPropType<RotationTypeTag, Properties::CouplingManager>;
    auto couplingManager = std::make_shared<CouplingManager>(rotationGridGeometry, deformationGridGeometry);

    using RotationProblem = GetPropType<RotationTypeTag, Properties::Problem>;
    auto rotationProblem = std::make_shared<RotationProblem>(rotationGridGeometry, couplingManager);
    using DeformationProblem = GetPropType<DeformationTypeTag, Properties::Problem>;
    auto deformationProblem = std::make_shared<DeformationProblem>(deformationGridGeometry, couplingManager);

    const auto& gg = *deformationGridGeometry;
    for (const auto& vertex : vertices(gg.gridView()))
    {
        const auto r = vertex.geometry().center().two_norm();
        if (r > 0.7 && r < 0.8)
        { deformationProblem->setPinDofIndex(gg.vertexMapper().index(vertex)); break; }
    }

    constexpr auto rotationIdx = CouplingManager::rotationIdx;
    constexpr auto deformationIdx = CouplingManager::deformationIdx;
    using Traits = MultiDomainTraits<RotationTypeTag, DeformationTypeTag>;
    using SolutionVector = typename Traits::SolutionVector;
    SolutionVector x;
    rotationProblem->applyInitialSolution(x[rotationIdx]);
    deformationProblem->applyInitialSolution(x[deformationIdx]);

    using RotationGridVariables = GetPropType<RotationTypeTag, Properties::GridVariables>;
    auto rotationGridVariables = std::make_shared<RotationGridVariables>(rotationProblem, rotationGridGeometry);
    using DeformationGridVariables = GetPropType<DeformationTypeTag, Properties::GridVariables>;
    auto deformationGridVariables = std::make_shared<DeformationGridVariables>(deformationProblem, deformationGridGeometry);

    couplingManager->init(rotationProblem, deformationProblem, x);
    rotationGridVariables->init(x[rotationIdx]);
    deformationGridVariables->init(x[deformationIdx]);

    using Assembler = MultiDomainFVAssembler<Traits, CouplingManager, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(
        std::make_tuple(rotationProblem, deformationProblem),
        std::make_tuple(rotationGridGeometry, deformationGridGeometry),
        std::make_tuple(rotationGridVariables, deformationGridVariables),
        couplingManager
    );
    using LAT = LinearAlgebraTraitsFromAssembler<Assembler>;
    using LinearSolver = UMFPackIstlSolver<SeqLinearSolverTraits, LAT>;
    auto linearSolver = std::make_shared<LinearSolver>();
    using NewtonSolver = Dumux::MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
    auto nonLinearSolver = std::make_shared<NewtonSolver>(assembler, linearSolver, couplingManager);

    const auto a = getParam<double>("Problem.InnerRadius");
    const auto b = getParam<double>("Problem.OuterRadius");
    const ClampedClampedAnnulus<double> exact(a, b, getParam<double>("Problem.Force"), deformationProblem->stiffness());

    double hMax = 0.0;
    for (const auto& element : elements(gg.gridView()))
        hMax = std::max(hMax, Dumux::diameter(element.geometry()));

    const auto errors = [&](const SolutionVector& sol)
    {
        double l2(0.0), l2Ref(0.0), l2Phi(0.0), l2PhiRef(0.0);
        for (const auto& element : elements(gg.gridView()))
        {
            const auto geometry = element.geometry();
            const auto elemSol = elementSolution(element, sol[deformationIdx], gg);
            const auto& quad = Dune::QuadratureRules<double, Grid::dimension>::rule(geometry.type(), 3);
            for (const auto& qp : quad)
            {
                const auto globalPos = geometry.global(qp.position());
                const auto values = evalSolution(element, geometry, gg, elemSol, globalPos);
                const auto dx = qp.weight()*geometry.integrationElement(qp.position());
                const auto wExact = exact.deflection(globalPos.two_norm());
                l2 += (values[1] - wExact)*(values[1] - wExact)*dx;
                l2Ref += wExact*wExact*dx;
                const auto phiExact = exact.shearGradPotential(globalPos.two_norm());
                l2Phi += (values[0] - phiExact)*(values[0] - phiExact)*dx;
                l2PhiRef += phiExact*phiExact*dx;
            }
        }
        return std::make_pair(std::sqrt(l2/l2Ref), std::sqrt(l2Phi/l2PhiRef));
    };

    const auto innerFlux = [&](const SolutionVector& sol)
    {
        deformationProblem->setNaturalInner(true);
        assembler->assembleResidual(sol);
        deformationProblem->setNaturalInner(false);
        double flux = 0.0;
        const auto& res = assembler->residual()[deformationIdx];
        for (const auto& vertex : vertices(gg.gridView()))
            if (deformationProblem->onInnerEdge(vertex.geometry().center()))
                flux += res[gg.vertexMapper().index(vertex)][1];
        return flux;
    };

    const auto solveWith = [&](double innerPotential)
    {
        deformationProblem->setInnerPotential(innerPotential);
        deformationProblem->applyInitialSolution(x[deformationIdx]);
        rotationProblem->applyInitialSolution(x[rotationIdx]);
        nonLinearSolver->solve(x);
        return x;
    };

    std::cout << "Max element diameter: " << hMax << std::endl;

    const auto x0 = solveWith(0.0);
    const auto [ew0, ephi0] = errors(x0);
    const auto c0 = innerFlux(x0);
    std::cout << "phi = 0 on the hole: rel L2 error w " << ew0 << ", phi " << ephi0
              << ", inner box flux of grad w - theta " << c0 << std::endl;

    const auto x1 = solveWith(1.0);
    const auto c1 = innerFlux(x1);
    const auto lambda = -c0/(c1 - c0);
    std::cout << "floating gauge: lambda = " << lambda << " (exact " << exact.innerPotential() << ")" << std::endl;

    const auto xs = solveWith(lambda);
    const auto [ews, ephis] = errors(xs);
    std::cout << "phi = lambda on the hole: inner box flux of grad w - theta " << innerFlux(xs) << std::endl;
    std::cout << "Relative L2-error deformation w: " << ews << std::endl;
    std::cout << "Relative L2-error shear gradient potential phi: " << ephis << std::endl;

    VtkOutputModule<DeformationGridVariables, std::tuple_element_t<1, SolutionVector>>
        vtkWriter(*deformationGridVariables, xs[deformationIdx], deformationProblem->name());
    vtkWriter.addVolumeVariable([](const auto& v){ return v.verticalDeformation(); }, "w");
    vtkWriter.addVolumeVariable([](const auto& v){ return v.shearCurlPotential(); }, "psi");
    vtkWriter.addVolumeVariable([](const auto& v){ return v.shearGradPotential(); }, "phi");
    vtkWriter.write(0.0);

    deformationProblem->setNaturalInner(true);
    deformationProblem->applyInitialSolution(x[deformationIdx]);
    rotationProblem->applyInitialSolution(x[rotationIdx]);
    nonLinearSolver->solve(x);
    const auto [ewn, ephin] = errors(x);
    std::cout << "compatibility natural on the hole: rel L2 error w " << ewn << ", phi " << ephin << std::endl;

    return 0;
}
