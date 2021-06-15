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
 * \brief Test for the staggered grid Navier-Stokes model (Donea 2003, \cite Donea2003).
 */

#include <config.h>

#include <ctime>
#include <iostream>

#include <dune/common/parallel/mpihelper.hh>
#include <dune/common/timer.hh>
#include <dune/grid/io/file/dgfparser/dgfexception.hh>
#include <dune/grid/io/file/vtk.hh>
#include <dune/istl/io.hh>

#include <dumux/assembly/fvassembler.hh>
#include <dumux/assembly/diffmethod.hh>
#include <dumux/common/dumuxmessage.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/io/grid/gridmanager.hh>
#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/io/vtkfunction.hh>
#include <dumux/linear/seqsolverbackend.hh>
#include <dumux/nonlinear/newtonsolver.hh>


#include <dumux/multidomain/fvassembler.hh>
#include <dumux/multidomain/traits.hh>
#include <dumux/multidomain/staggeredfreeflow/couplingmanager.hh>
#include <dumux/multidomain/newtonsolver.hh>

#include <dumux/io/vtk/intersectionwriter.hh>
#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/freeflow/navierstokes/velocityoutput.hh>

#include "../l2error.hh"
#include "../analyticalsolution.hh"
#include "problem_new.hh"

namespace Dumux::Properties{

// Set the problem property
template<class TypeTag>
struct CouplingManager<TypeTag, TTag::DoneaTestNew>
{
private:
    using Traits = MultiDomainTraits<TTag::DoneaTestNewMomentum, TTag::DoneaTestNewMass>;
public:
    using type = StaggeredFreeFlowCouplingManager<Traits>;
};

}

template<class Problem, class GridVolumeVariables, class SolutionVector>
auto getVelocityGradient(const Problem& problem, const GridVolumeVariables& gridVolVars, const SolutionVector& sol)
{
    const auto& gridGeometry = problem.gridGeometry();
    struct Result
    {
        std::array<std::vector<double>, 4> analytical;
        std::array<std::vector<double>, 4> numerical;
    } result;

    for (auto& r : result.analytical)
        r.resize(gridGeometry.gridView().size(0));

    for (auto& r : result.numerical)
        r.resize(gridGeometry.gridView().size(0));


    auto fvGeometry = localView(gridGeometry);
    auto elemVolVars = localView(gridVolVars);
    for (const auto& element : elements(gridGeometry.gridView()))
    {
        fvGeometry.bind(element);
        elemVolVars.bind(element, fvGeometry, sol);
        const auto eIdx = gridGeometry.elementMapper().index(element);

        const auto gradAnalytical = problem.velocityGradient(element.geometry().center());
        result.analytical[0][eIdx] = gradAnalytical[0][0];
        result.analytical[1][eIdx] = gradAnalytical[0][1];
        result.analytical[2][eIdx] = gradAnalytical[1][0];
        result.analytical[3][eIdx] = gradAnalytical[1][1];

        for (const auto& scvf : scvfs(fvGeometry))
        {
            if (scvf.isFrontal() && !scvf.boundary())
            {
                const auto gradNumerical = Dumux::StaggeredVelocityGradients::velocityGradient(fvGeometry, scvf, elemVolVars, true);
                result.numerical[0][eIdx] = gradNumerical[0][0];
                result.numerical[1][eIdx] = gradNumerical[0][1];
                result.numerical[2][eIdx] = gradNumerical[1][0];
                result.numerical[3][eIdx] = gradNumerical[1][1];
            }
        }
    }

    return result;
}

int main(int argc, char** argv)
{
    using namespace Dumux;

    // define the type tag for this problem
    using MomentumTypeTag = Properties::TTag::DoneaTestNewMomentum;
    using MassTypeTag = Properties::TTag::DoneaTestNewMass;

    // initialize MPI, finalize is done automatically on exit
    const auto& mpiHelper = Dune::MPIHelper::instance(argc, argv);

    // print dumux start message
    if (mpiHelper.rank() == 0)
        DumuxMessage::print(/*firstCall=*/true);

    // parse command line arguments and input file
    Parameters::init(argc, argv);

    // try to create a grid (from the given grid file or the input file)
    using GridManager = Dumux::GridManager<GetPropType<MomentumTypeTag, Properties::Grid>>;
    GridManager gridManager;
    gridManager.init();

    ////////////////////////////////////////////////////////////
    // run instationary non-linear problem on this grid
    ////////////////////////////////////////////////////////////

    // we compute on the leaf grid view
    const auto& leafGridView = gridManager.grid().leafGridView();

    // create the finite volume grid geometry
    Dune::Timer timer;
    using MomentumGridGeometry = GetPropType<MomentumTypeTag, Properties::GridGeometry>;
    auto momentumGridGeometry = std::make_shared<MomentumGridGeometry>(leafGridView);
    momentumGridGeometry->update();

    using MassGridGeometry = GetPropType<MassTypeTag, Properties::GridGeometry>;
    auto massGridGeometry = std::make_shared<MassGridGeometry>(leafGridView);
    massGridGeometry->update();

    // the coupling manager
    using Traits = MultiDomainTraits<MomentumTypeTag, MassTypeTag>;
    using CouplingManager = StaggeredFreeFlowCouplingManager<Traits>;

    auto couplingManager = std::make_shared<CouplingManager>();

    // the problem (boundary conditions)
    using MomentumProblem = GetPropType<MomentumTypeTag, Properties::Problem>;
    auto momentumProblem = std::make_shared<MomentumProblem>(momentumGridGeometry, couplingManager);

    using MassProblem = GetPropType<MassTypeTag, Properties::Problem>;
    auto massProblem = std::make_shared<MassProblem>(massGridGeometry, couplingManager);

    // the solution vector
    constexpr auto momentumIdx = Dune::index_constant<0>();
    constexpr auto massIdx = Dune::index_constant<1>();
    using SolutionVector = typename Traits::SolutionVector;
    SolutionVector x;
    x[momentumIdx].resize(momentumGridGeometry->numDofs());
    x[massIdx].resize(massGridGeometry->numDofs());

    // the grid variables
    using MomentumGridVariables = GetPropType<MomentumTypeTag, Properties::GridVariables>;
    auto momentumGridVariables = std::make_shared<MomentumGridVariables>(momentumProblem, momentumGridGeometry);

    using MassGridVariables = GetPropType<MassTypeTag, Properties::GridVariables>;
    auto massGridVariables = std::make_shared<MassGridVariables>(massProblem, massGridGeometry);

    couplingManager->init(momentumProblem, massProblem, std::make_tuple(momentumGridVariables, massGridVariables), x);
    massGridVariables->init(x[massIdx]);
    momentumGridVariables->init(x[momentumIdx]);

    using Assembler = MultiDomainFVAssembler<Traits, CouplingManager, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(std::make_tuple(momentumProblem, massProblem),
                                                 std::make_tuple(momentumGridGeometry, massGridGeometry),
                                                 std::make_tuple(momentumGridVariables, massGridVariables),
                                                 couplingManager);

    // intialize the vtk output module
    using IOFields = GetPropType<MassTypeTag, Properties::IOFields>;
    VtkOutputModule vtkWriter(*massGridVariables, x[massIdx], massProblem->name());
    IOFields::initOutputModule(vtkWriter); // Add model specific output fields
    vtkWriter.addVelocityOutput(std::make_shared<NavierStokesVelocityOutput<MassGridVariables>>());
    const auto exactPressure = getScalarAnalyticalSolution(*massProblem)[GetPropType<MassTypeTag, Properties::ModelTraits>::Indices::pressureIdx];
    const auto exactVelocity = getVelocityAnalyticalSolution(*momentumProblem);
    vtkWriter.addField(exactPressure, "pressureExact");
    vtkWriter.addField(exactVelocity, "velocityExact");


    // the linear solver
    using LinearSolver = Dumux::UMFPackBackend;
    auto linearSolver = std::make_shared<LinearSolver>();

    // the non-linear solver
    using NewtonSolver = MultiDomainNewtonSolver<Assembler, LinearSolver, CouplingManager>;
    NewtonSolver nonLinearSolver(assembler, linearSolver, couplingManager);

    // linearize & solve
    nonLinearSolver.solve(x);

    static const bool addGradients = getParam<bool>("Vtk.AddGradients", false);
    if (addGradients)
    {
        const auto gradV = getVelocityGradient(*momentumProblem, momentumGridVariables->curGridVolVars(), x[momentumIdx]);
        vtkWriter.addField(gradV.analytical[0], "dudxExact");
        vtkWriter.addField(gradV.analytical[1], "dudyExact");
        vtkWriter.addField(gradV.analytical[2], "dvdxExact");
        vtkWriter.addField(gradV.analytical[3], "dvdyExact");
        vtkWriter.addField(gradV.numerical[0], "dudx");
        vtkWriter.addField(gradV.numerical[1], "dudy");
        vtkWriter.addField(gradV.numerical[2], "dvdx");
        vtkWriter.addField(gradV.numerical[3], "dvdy");
        vtkWriter.write(1.0);
    }
    else
        vtkWriter.write(1.0);

    using GridView = std::decay_t<decltype(momentumGridGeometry->gridView())>;
    using SolVec = std::decay_t<decltype(x[momentumIdx])>;

    class ScalarFunction
    {
        public:
        typedef Dune::VTK::SkeletonFunctionTraits<std::decay_t<decltype(momentumGridGeometry->gridView())>, typename std::decay_t<decltype(momentumGridGeometry->gridView())>::ctype>
        Traits;

        ScalarFunction(const GridView& gv, const SolVec& s) : gv_(gv), s_(s) {}

        //! return number of components
        unsigned dimRange() const { return 1; }

        void evaluate(const typename Traits::Cell& c,
                        const typename Traits::Domain& xl,
                        typename Traits::Range& result) const
        {
            assert(c.conforming());


            const auto globalIdx = gv_.indexSet().subIndex(c.inside(), c.indexInInside(), 1);
            result.resize(1, s_[globalIdx]);
        }

        private:
        const GridView gv_;
        const SolVec& s_;
    };

    auto faceData = std::make_shared<ScalarFunction>(momentumGridGeometry->gridView(), x[momentumIdx]);
    ConformingIntersectionWriter vtk(momentumGridGeometry->gridView());
    vtk.addCellData(faceData, "velcocityScalar");
    vtk.write("facedata", Dune::VTK::ascii);

    if (getParam<bool>("Problem.PrintL2Error"))
    {
        const auto pressureL2error = calculateL2Error(*massProblem, x[massIdx]);
        const auto velocityL2error = calculateL2Error(*momentumProblem, x[momentumIdx]);

        std::cout << std::setprecision(8) << "** L2 error (abs/rel) for "
                        << std::setw(6) << massGridGeometry->numDofs() << " cc dofs and " << momentumGridGeometry->numDofs()
                        << " face dofs (total: " << massGridGeometry->numDofs() + momentumGridGeometry->numDofs() << "): "
                        << std::scientific
                        << "L2(p) = " << pressureL2error.absolute[0] << " / " << pressureL2error.relative[0]
                        << " , L2(vx) = " << velocityL2error.absolute[0] << " / " << velocityL2error.relative[0]
                        << " , L2(vy) = " << velocityL2error.absolute[1] << " / " << velocityL2error.relative[1]
                        << std::endl;
    }

    timer.stop();

    const auto& comm = Dune::MPIHelper::getCollectiveCommunication();
    std::cout << "Simulation took " << timer.elapsed() << " seconds on "
              << comm.size() << " processes.\n"
              << "The cumulative CPU time was " << timer.elapsed()*comm.size() << " seconds.\n";

    ////////////////////////////////////////////////////////////
    // finalize, print dumux message to say goodbye
    ////////////////////////////////////////////////////////////

    // print dumux end message
    if (mpiHelper.rank() == 0)
    {
        Parameters::print();
        DumuxMessage::print(/*firstCall=*/false);
    }

    return 0;
}
