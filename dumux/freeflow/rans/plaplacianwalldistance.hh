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
 * \ingroup RANSModel
 * \copydoc Dumux::PLaplaceWallDistance
 */
#ifndef DUMUX_RANS_P_LAPLACE_WALL_DISTANCE
#define DUMUX_RANS_P_LAPLACE_WALL_DISTANCE

#include <vector>
#include <functional>
#include <dumux/assembly/fvassembler.hh>
#include <dumux/common/indextraits.hh>
#include <dumux/common/fvproblem.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/discretization/evalgradients.hh>
#include <dumux/discretization/evalsolution.hh>
#include <dumux/discretization/fvgridvariables.hh>
#include <dumux/discretization/box.hh>
#include <dumux/discretization/box/fluxvariablescache.hh>
#include <dumux/flux/fluxvariablescaching.hh>
#include <dumux/nonlinear/newtonsolver.hh>
#include <dumux/linear/seqsolverbackend.hh>

namespace Dumux::WallPropertiesDetails {

template<class Scalar>
class VolumeVariables
{
public:
    using PrimaryVariables = Dune::FieldVector<Scalar, 1>;

    template<class ElemSol, class Problem, class Element, class Scv>
    void update(const ElemSol& elemSol,
                const Problem& problem,
                const Element& element,
                const Scv& scv)
    { priVars_ = elemSol[scv.localDofIndex()]; }

    const PrimaryVariables &priVars() const
    { return priVars_; }

    Scalar priVar(const int pvIdx) const
    { return priVars_[pvIdx]; }

    Scalar distance() const
    { return priVars_[0]; }

    Scalar extrusionFactor() const
    { return 1.0; }

private:
    PrimaryVariables priVars_;
};

template<class TypeTag>
class Problem : public FVProblem<TypeTag>
{
    using ParentType =  FVProblem<TypeTag>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<GetPropType<TypeTag, Properties::PrimaryVariables>>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using BoundaryTypes = Dumux::BoundaryTypes<PrimaryVariables::size()>;

public:

    template<class CondsiderWallPos>
    Problem(std::shared_ptr<const GridGeometry> gridGeometry,
            const CondsiderWallPos& condsiderWallPos,
            const std::string& paramGroup = "")
    : ParentType(gridGeometry, paramGroup)
    {
        considerPos_ = condsiderWallPos;
        restrictonSet_ = true;
    }

    Problem(std::shared_ptr<const GridGeometry> gridGeometry,
            const std::string& paramGroup = "")
    : ParentType(gridGeometry, paramGroup)
    {
        considerPos_ = [](const GlobalPosition& pos) { return true; };
    }

    bool solveNonLinear() const
    { return solveNonLinear_; }

    void switchToNonLinearSolve(bool b)
    { solveNonLinear_ = b; }

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes bcTypes;
        bcTypes.setAllDirichlet();

        if (restrictonSet_ && !considerPos_(globalPos))
            bcTypes.setAllNeumann();

        return bcTypes;
    }

    PrimaryVariables dirichletAtPos(const GlobalPosition& p) const
    { return PrimaryVariables(0.0); }

    PrimaryVariables sourceAtPos(const GlobalPosition& p) const
    { return PrimaryVariables(-1.0); }

private:
    bool solveNonLinear_ = false;
    bool restrictonSet_ = false;
    std::function<bool (const GlobalPosition& pos)> considerPos_;
};

template<class TypeTag>
class LocalResidual : public GetPropType<TypeTag, Properties::BaseLocalResidual>
{
    using ParentType = GetPropType<TypeTag, Properties::BaseLocalResidual>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using VolumeVariables = GetPropType<TypeTag, Properties::VolumeVariables>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using NumEqVector = Dumux::NumEqVector<GetPropType<TypeTag, Properties::PrimaryVariables>>;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using ElementVolumeVariables = typename GetPropType<TypeTag, Properties::GridVolumeVariables>::LocalView;
    using ElementFluxVariablesCache = typename GetPropType<TypeTag, Properties::GridFluxVariablesCache>::LocalView;

    static constexpr auto dimWorld = GridGeometry::GridView::dimensionworld;

public:
    using ParentType::ParentType;

    NumEqVector computeStorage(const Problem& problem,
                               const SubControlVolume& scv,
                               const VolumeVariables& volVars) const
    { return NumEqVector(0.0); }

    NumEqVector computeFlux(const Problem& problem,
                            const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolumeFace& scvf,
                            const ElementFluxVariablesCache& elemFluxVarsCache) const
    {
        // evaluate gradX at integration point and interpolate density
        const auto& fluxVarsCache = elemFluxVarsCache[scvf];

        Dune::FieldVector<Scalar, dimWorld> gradP(0.0);
        for (auto&& scv : scvs(fvGeometry))
        {
            const auto& volVars = elemVolVars[scv];

            // the global shape function gradient
            gradP.axpy(volVars.distance(), fluxVarsCache.gradN(scv.indexInElement()));
        }

        Scalar factor = 1.0;

        if (problem.solveNonLinear())
        {
            const auto absGrad = gradP.two_norm();
            static const auto order = getParamFromGroup<int>(problem.paramGroup(), "Problem.Order", 2);
            factor *= Dune::power(absGrad, order-2);
        }

        return vtmv(scvf.unitOuterNormal(), factor, gradP)*scvf.area();
    }
};

}

namespace Dumux::Properties {

// Create new type tags
namespace TTag {

template<class G>
struct PoissonWallProblem
{
    using InheritsFrom = std::tuple<BoxModel>;
    using Grid = G;
};
} // end namespace TTag

template<class TypeTag, class T>
struct Grid<TypeTag, TTag::PoissonWallProblem<T>> { using type = typename TypeTag::Grid; };

template<class TypeTag, class T>
struct Scalar<TypeTag, TTag::PoissonWallProblem<T>> { using type = double; };

template<class TypeTag, class T>
struct VolumeVariables<TypeTag, TTag::PoissonWallProblem<T>> { using type = Dumux::WallPropertiesDetails::VolumeVariables<double>; };

template<class TypeTag, class T>
struct Problem<TypeTag, TTag::PoissonWallProblem<T>> { using type = Dumux::WallPropertiesDetails::Problem<TypeTag>; };

template<class TypeTag, class T>
struct LocalResidual<TypeTag, TTag::PoissonWallProblem<T>> { using type = Dumux::WallPropertiesDetails::LocalResidual<TypeTag>; };

template<class TypeTag, class T>
struct PrimaryVariables<TypeTag, TTag::PoissonWallProblem<T>> { using type = Dune::FieldVector<double, 1>; };

//! The flux variables cache class for models involving flow in porous media
template<class TypeTag, class T>
struct FluxVariablesCache<TypeTag, TTag::PoissonWallProblem<T>> { using type = Dumux::BoxFluxVariablesCache<GetPropType<TypeTag, Properties::Scalar>, GetPropType<TypeTag, Properties::GridGeometry>>; };

//! The flux variables cache filler (FluxVariablesCache is the data type,
//! the filler knows how to build up the caches for the stencil efficiently)
template<class TypeTag, class T>
struct FluxVariablesCacheFiller<TypeTag, TTag::PoissonWallProblem<T>> { using type = Dumux::FluxVariablesCaching::EmptyCacheFiller; };

template<class TypeTag, class T>
struct ModelTraits<TypeTag, TTag::PoissonWallProblem<T>>
{
    struct type
    {
        static constexpr auto numEq() { return 1; }
    };
};

} // end namespace Properties;

namespace Dumux {

/*!
 * \ingroup RANSModel
 * \brief Class to calculate the wall distance based on a p-Laplace distance map
 *         See http://home.eps.hw.ac.uk/~ab226/papers/pde4dist.pdf
 */
template<class GridGeometry>
class PLaplaceWallDistance
{
    using GridView = typename GridGeometry::GridView;
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using ScvfGeometry = std::decay_t<decltype(std::declval<SubControlVolumeFace>().geometry())>;
    using Scalar = typename GridView::Grid::ctype;

public:
    PLaplaceWallDistance(const GridGeometry& gridGeometry, const std::string& paramGroup = "")
    : gridGeometry_(gridGeometry)
    , paramGroup_(paramGroup) {}

    void updateWallDistance()
    {
        updateWallDistanceImpl_([](const auto& p) { return true; }, false);
    }

    template<class ConsiderPosFunction>
    void updateWallDistance(const ConsiderPosFunction& considerPos)
    {
        updateWallDistanceImpl_(considerPos, true);
    }

    template<class ConsiderPosFunction>
    void updateWallDistanceImpl_(const ConsiderPosFunction& considerPos, bool restriction)
    {
        distance_.clear();

        using TypeTag = Properties::TTag::PoissonWallProblem<typename GridGeometry::Grid>;

        using GG = GetPropType<TypeTag, Properties::GridGeometry>;
        using Problem = GetPropType<TypeTag, Properties::Problem>;

        auto boxGG = std::make_shared<GG>(gridGeometry_.gridView());
        boxGG->update();
        auto problem = restriction ? std::make_shared<Problem>(boxGG, considerPos, paramGroup_)
                                   : std::make_shared<Problem>(boxGG, paramGroup_);

        using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
        auto gridVariables = std::make_shared<GridVariables>(problem, boxGG);

        // the solution vector
        using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
        SolutionVector x(boxGG->numDofs());
        gridVariables->init(x);

        // the assembler with time loop for instationary problem
        using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
        auto assembler = std::make_shared<Assembler>(problem, boxGG, gridVariables);

        // the linear solver
        using LinearSolver = UMFPackBackend;
        auto linearSolver = std::make_shared<LinearSolver>();

        // the non-linear solver
        using NewtonSolver = NewtonSolver<Assembler, LinearSolver>;
        NewtonSolver nonLinearSolver(assembler, linearSolver);

        nonLinearSolver.solve(x);

        problem->switchToNonLinearSolve(true);

        nonLinearSolver.solve(x);

        std::cout << "done " << std::endl;

        auto fvGeometry = localView(*boxGG);

        auto dist = x;

        static const auto p = getParamFromGroup<int>(paramGroup_, "Problem.Order", 2);

        for (const auto& element : elements(boxGG->gridView()))
        {
            fvGeometry.bindElement(element);

            const auto& elemSol = elementSolution(element, x, *boxGG);

            for (const auto& scv : scvs(fvGeometry))
            {
                const auto grad = evalGradients(element,
                                                element.geometry(),
                                                *boxGG,
                                                elemSol,
                                                scv.dofPosition())[0];

                auto result = -1.0 * Dune::power(grad.two_norm(), p-1);
                const auto f = (static_cast<double>(p)/static_cast<double>(p-1) * x[scv.dofIndex()] + Dune::power(grad.two_norm(), p));
                result += std::pow(f, static_cast<double>(p-1)/static_cast<double>(p));
                dist[scv.dofIndex()] = result;
            }
        }

        if constexpr (GridGeometry::discMethod == DiscretizationMethod::box)
        {
            distance_.resize(dist.size());
            for (int i = 0; i < distance_.size(); ++i)
                distance_[i] = dist[i][0];
        }
        else
        {
            distance_.resize(gridGeometry_.gridView().size(0));
            for (const auto& element : elements(boxGG->gridView()))
            {
                const auto eIdx = boxGG->elementMapper().index(element);
                const auto& elementGeometry = element.geometry();
                const auto elemSol = elementSolution(element, dist, *boxGG);
                distance_[eIdx] = evalSolution(element,
                                               elementGeometry,
                                               *boxGG,
                                               elemSol,
                                               elementGeometry.center());
            }
        }
    }

    const std::vector<Scalar>& wallDinstance() const
    { return distance_; }

private:
    std::vector<Scalar> distance_;
    const GridGeometry& gridGeometry_;
    const std::string paramGroup_;
};

} // end namespace Dumux

#endif
