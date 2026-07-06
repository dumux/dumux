// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//

#include <config.h>

#include <type_traits>
#include <complex>

#include <dune/grid/yaspgrid.hh>

#include <dumux/common/initialize.hh>
#include <dumux/common/properties.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/numeqvector.hh>
#include <dumux/common/boundarytypes.hh>
#include <dumux/common/fvproblem.hh>

#include <dumux/io/gridwriter.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>

#include <dumux/discretization/box.hh>

#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/pdesolver.hh>
#include <dumux/assembly/fvassembler.hh>

#include "model.hh"

namespace Dumux {
template<class TypeTag>
class ComplexHelmholtzTestProblem : public FVProblem<TypeTag>
{
    using ParentType = FVProblem<TypeTag>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using GlobalPosition = typename GridGeometry::LocalView::Element::Geometry::GlobalCoordinate;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;
    using NumEqVector = Dumux::NumEqVector<PrimaryVariables>;
    using BoundaryTypes = Dumux::BoundaryTypes<GetPropType<TypeTag, Properties::ModelTraits>::numEq()>;
    using PointSource = GetPropType<TypeTag, Properties::PointSource>;

public:
    ComplexHelmholtzTestProblem(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {
        waveNumberSquared_ = {M_PI*M_PI*2.0, 0.0};
    }

    BoundaryTypes boundaryTypesAtPos(const GlobalPosition& globalPos) const
    {
        BoundaryTypes values;
        values.setAllDirichlet();
        return values;
    }

    PrimaryVariables dirichletAtPos(const GlobalPosition& globalPos) const
    { return PrimaryVariables(0.0); }

    NumEqVector sourceAtPos(const GlobalPosition& globalPos) const
    {
        NumEqVector source(0.0);

        const double R = 0.1;
        if (std::hypot(globalPos[0]-0.37, globalPos[1]-0.43) < R)
            source[0] = 1.0/(M_PI*R*R);

        return source;
    }

    void setWaveNumberSquared(const std::complex<double>& waveNumberSquared)
    { waveNumberSquared_ = waveNumberSquared; }

    std::complex<double> waveNumberSquared() const
    { return waveNumberSquared_; }

private:
    std::complex<double> waveNumberSquared_;
};
} // end namespace Dumux

namespace Dumux::Properties::TTag {

struct ComplexHelmholtzTest
{
    using InheritsFrom = std::tuple<ComplexHelmholtzModel, BoxModel>;

    using Scalar = double;
    using Grid = Dune::YaspGrid<2>;

    template<class TypeTag>
    using Problem = ComplexHelmholtzTestProblem<TypeTag>;

    using EnableGridVolumeVariablesCache = std::true_type;
    using EnableGridFluxVariablesCache = std::true_type;
    using EnableGridGeometryCache = std::true_type;
};

} // end namespace Dumux::Properties::TTag

int main(int argc, char** argv)
{
    using namespace Dumux;

    Dumux::initialize(argc, argv);

    Parameters::init(argc, argv);

    using TypeTag = Properties::TTag::ComplexHelmholtzTest;

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Grid = GetPropType<TypeTag, Properties::Grid>;
    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;

    GridManager<Grid> gridManager;
    gridManager.init();

    auto gridGeometry = std::make_shared<GridGeometry>(gridManager.grid().leafGridView());
    auto problem = std::make_shared<Problem>(gridGeometry);
    SolutionVector sol(gridGeometry->dofMapper().size());
    sol = std::complex<double>(0.0);
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(sol);

    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    using LinearSolver = ILUBiCGSTABIstlSolver<LinearSolverTraits<GridGeometry>, LinearAlgebraTraitsFromAssembler<Assembler>>;
    using Solver = Dumux::LinearPDESolver<Assembler, LinearSolver>;

    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables);
    auto linearSolver = std::make_shared<LinearSolver>(gridGeometry->gridView(), gridGeometry->dofMapper());
    Solver solver(assembler, linearSolver);

    // do a frequency sweep over main resonance frequencies of the domain
    for (int m = 1; m < 10; ++m)
    {
        for (int n = 1; n <= m; ++n)
        {
            const auto c = 1.0; // wave speed
            const auto k = M_PI*std::sqrt(m*m + n*n);
            const auto omega = c*k; // angular frequency
            const auto gamma = 0.1; // damping
            const auto kSq = std::complex<double>(omega*omega/(c*c), -gamma*omega); // wave number squared
            problem->setWaveNumberSquared(kSq);

            solver.solve(sol);

            // for the output, we sample the solution over one period of the oscillation
            std::vector<double> solutionAtSamplingPoints(sol.size(), 0.0);
            IO::GridWriter<typename Grid::LeafGridView> vtkWriter(
                IO::Format::pvd_with(IO::Format::vtu.with({
                    .encoder = IO::Encoding::ascii,
                    .compressor = IO::Compression::none,
                    .data_format = IO::VTK::DataFormat::inlined
                })),
                gridGeometry->gridView(),
                problem->name() + "_" + std::to_string(m) + "_" + std::to_string(n)
            );

            vtkWriter.setPointField("u", [&](const auto& v){ return solutionAtSamplingPoints[gridGeometry->dofMapper().index(v)]; });

            const int samplingPointsPerPeriod = 50;
            for (int j = 0; j < samplingPointsPerPeriod; ++j)
            {
                const auto t = j*2.0*M_PI/(omega*samplingPointsPerPeriod);
                for (size_t idx = 0; idx < sol.size(); ++idx)
                {
                    using namespace std::complex_literals;
                    solutionAtSamplingPoints[idx] = std::real(sol[idx][0]*std::exp(1i*omega*t));
                }
                vtkWriter.write(t*omega); // phase as time value
            }
        }
    }

    return 0;
}
