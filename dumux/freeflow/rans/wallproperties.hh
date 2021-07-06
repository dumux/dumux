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
 * \copydoc Dumux::RANSIOFields
 */
#ifndef DUMUX_RANS_WALL_PROPERTIES_HH
#define DUMUX_RANS_WALL_PROPERTIES_HH

#include <vector>
#include <dumux/common/indextraits.hh>
#include <dumux/common/fvproblem.hh>
#include <dumux/geometry/distance.hh>
#include <dumux/discretization/fvgridvariables.hh>
#include <dumux/discretization/box.hh>
#include <dumux/discretization/evalgradients.hh>
#include <dumux/flux/fluxvariablescaching.hh>
#include <dumux/discretization/box/fluxvariablescache.hh>
#include <dumux/assembly/fvassembler.hh>
#include <dumux/nonlinear/newtonsolver.hh>
#include <dumux/linear/seqsolverbackend.hh>

namespace Dumux {

/*!
 * \ingroup RANSModel
 * \brief Class to calculate the wall distance based on a simple brute force search
 */
template<class GridGeometry>
class BoundarySearchWallProperties
{
    using GridView = typename GridGeometry::GridView;
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using ScvfGeometry = std::decay_t<decltype(std::declval<SubControlVolumeFace>().geometry())>;
    using Scalar = typename GridView::Grid::ctype;

    struct WallScvfData
    {
        GridIndexType eIdx;
        GridIndexType scvfIdx;
    };

public:

    BoundarySearchWallProperties(const GridGeometry& gridGeometry)
    : gridGeometry_(gridGeometry)  {}

    void updateWallDistance()
    {
        updateWallDistance([](const auto& scvf) { return true; });
    }

    template<class ConsiderFaceFunction>
    void updateWallDistance(const ConsiderFaceFunction& considerFace)
    {
        struct TempScvfInfo
        {
            ScvfGeometry geometry;
            GridIndexType eIdx;
            GridIndexType scvfIdx;
        };

        std::vector<TempScvfInfo> tempScvfInfo;
        tempScvfInfo.reserve(gridGeometry_.numBoundaryScvf());

        wallScvfData_.clear();
        wallScvfData_.resize(gridGeometry_.numDofs());
        distance_.clear();
        distance_.resize(gridGeometry_.numDofs(), std::numeric_limits<Scalar>::max());

        auto fvGeometry = localView(gridGeometry_);

        // first loop over all elements: find all wall scvfs
        for (const auto& element : elements(gridGeometry_.gridView()))
        {
            fvGeometry.bindElement(element);
            if (!fvGeometry.hasBoundaryScvf())
                continue;

            const auto eIdx = gridGeometry_.elementMapper().index(element);

            for (const auto& scvf : scvfs(fvGeometry))
            {
                if (scvf.boundary() && considerFace(scvf))
                    tempScvfInfo.push_back({scvf.geometry(), eIdx, scvf.index()});
            }
        }

        // second loop over all elements: find distances of all dofs to wall
        for (const auto& element : elements(gridGeometry_.gridView()))
        {
            fvGeometry.bindElement(element);

            for (const auto& scv : scvs(fvGeometry))
            {
                for (const auto& wallScvfInfo : tempScvfInfo)
                {
                    // TODO 3D
                    const auto d = distancePointSegment(scv.dofPosition(), wallScvfInfo.geometry);
                    if (d < distance_[scv.dofIndex()])
                    {
                        distance_[scv.dofIndex()] = d;
                        wallScvfData_[scv.dofIndex()] = WallScvfData{wallScvfInfo.eIdx, wallScvfInfo.scvfIdx};
                    }
                }
            }
        }
    }

    const std::vector<Scalar>& wallDinstance() const
    { return distance_; }

private:
    std::vector<Scalar> distance_;
    std::vector<WallScvfData> wallScvfData_;
    const GridGeometry& gridGeometry_;
};

}


namespace Dumux::WallPropertiesDetails {

template<class Scalar>
class VolumeVariables
{

public:
    //! Export the type used for the primary variables
    using PrimaryVariables = Dune::FieldVector<Scalar, 1>;

    /*!
    * \brief Updates all quantities for a given control volume.
    *
    * \param elemSol A vector containing all primary variables connected to the element
    * \param problem The object specifying the problem which ought to
    *                be simulated
    * \param element An element which contains part of the control volume
    * \param scv The sub-control volume
    */
    template<class ElemSol, class Problem, class Element, class Scv>
    void update(const ElemSol& elemSol,
                const Problem& problem,
                const Element& element,
                const Scv& scv)
    {
        priVars_ = elemSol[scv.localDofIndex()];
    }

    /*!
    * \brief Returns the vector of primary variables.
    */
    const PrimaryVariables &priVars() const
    { return priVars_; }

    /*!
    * \brief Returns a component of primary variable vector.
    *
    * \param pvIdx The index of the primary variable of interest
    */
    Scalar priVar(const int pvIdx) const
    { return priVars_[pvIdx]; }

    Scalar distance() const
    { return priVars_[0]; }

    /*!
    * \brief Returns how much the sub-control volume is extruded.
    *
    * This means the factor by which a lower-dimensional (1D or 2D)
    * entity needs to be expanded to get a full dimensional cell. The
    * default is 1.0 which means that 1D problems are actually
    * thought as pipes with a cross section of 1 m^2 and 2D problems
    * are assumed to extend 1 m to the back.
    */
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
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    using ParentType::ParentType;

    bool solveNonLinear() const
    { return solveNonLinear_; }

    void switchToNonLinearSolve(bool b)
    { solveNonLinear_ = b; }

    PrimaryVariables dirichletAtPos(const GlobalPosition& p) const
    {
        return PrimaryVariables(0.0);
    }

    PrimaryVariables sourceAtPos(const GlobalPosition& p) const
    {
        return PrimaryVariables(-1.0);
    }


private:
    bool solveNonLinear_ = false;
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
    {
        return NumEqVector(0.0);
    }

    NumEqVector computeFlux(const Problem& problem,
                            const Element& element,
                            const FVElementGeometry& fvGeometry,
                            const ElementVolumeVariables& elemVolVars,
                            const SubControlVolumeFace& scvf,
                            const ElementFluxVariablesCache& elemFluxVarsCache) const
    {
        const auto& insideVolVars = elemVolVars[scvf.insideScvIdx()];
        const auto& outsideVolVars = elemVolVars[scvf.outsideScvIdx()];

        // evaluate gradX at integration point and interpolate density
        const auto& fluxVarsCache = elemFluxVarsCache[scvf];
        const auto& shapeValues = fluxVarsCache.shapeValues();

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
            static const auto order = getParam<int>("Problem.Order", 2);
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
 * \brief Class to calculate the wall distance based on a Poisson distance map
 */
template<class GridGeometry>
class PoissonWallProperties
{
    using GridView = typename GridGeometry::GridView;
    using GridIndexType = typename IndexTraits<GridView>::GridIndex;
    using SubControlVolumeFace = typename GridGeometry::SubControlVolumeFace;
    using ScvfGeometry = std::decay_t<decltype(std::declval<SubControlVolumeFace>().geometry())>;
    using Scalar = typename GridView::Grid::ctype;


    struct WallScvfData
    {
        GridIndexType eIdx;
        GridIndexType scvfIdx;
    };

public:
    PoissonWallProperties(const GridGeometry& gridGeometry)
    : gridGeometry_(gridGeometry) {}

    void solve()
    {
        using TypeTag = Properties::TTag::PoissonWallProblem<typename GridGeometry::Grid>;

        using GG = GetPropType<TypeTag, Properties::GridGeometry>;
        using GridVolVars = GetPropType<TypeTag, Properties::GridVolumeVariables>;
        using Problem = GetPropType<TypeTag, Properties::Problem>;

        auto boxGG = std::make_shared<GG>(gridGeometry_.gridView());
        boxGG->update();
        auto problem = std::make_shared<Problem>(boxGG);

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

        Dune::VTKWriter<typename GG::GridView> writer(boxGG->gridView());


        auto fvGeometry = localView(*boxGG);

        auto dist = x;

        static const auto p = getParam<int>("Problem.Order", 2);


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

        writer.addVertexData(dist, "distance");

        writer.write("result_poisson");


    }

private:
    std::vector<Scalar> distance_;
    std::vector<WallScvfData> wallScvfData_;
    const GridGeometry& gridGeometry_;
};

} // end namespace Dumux

#endif
