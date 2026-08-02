// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \brief Checks that restricting the assembly element loop to a subdomain reproduces the corresponding
 *        rows of the monolithic Jacobian and residual exactly.
 *
 * Looping \f$ N_i \f$ must give the square block, and extending the loop to \f$ N_i \cup \Gamma_i \f$
 * must additionally give the columns outside \f$ N_i \f$. Both are compared entry-wise against an
 * ordinary FVAssembler at the same, deliberately non-converged solution.
 */
#include <config.h>

#include <cmath>
#include <iostream>
#include <memory>
#include <vector>

#include <dune/common/parallel/mpihelper.hh>

#include <dumux/assembly/diffmethod.hh>
#include <dumux/assembly/fvassembler.hh>
#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/io/grid/gridmanager_yasp.hh>
#include <dumux/linear/dunevectors.hh>

#include <dumux/nonlinear/preconditioning/assembler.hh>
#include <dumux/nonlinear/preconditioning/partition.hh>

#include "2p/properties.hh"

namespace {

int failures = 0;

void check(bool condition, const std::string& what)
{
    if (!condition)
    {
        std::cerr << "FAIL: " << what << std::endl;
        ++failures;
    }
}

} // end anonymous namespace

int main(int argc, char** argv)
{
    using namespace Dumux;
    using TypeTag = Properties::TTag::TwoPIncompressibleTpfa;

    initialize(argc, argv);
    Parameters::init(argc, argv);

    using Grid = GetPropType<TypeTag, Properties::Grid>;
    GridManager<Grid> gridManager;
    gridManager.init();

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(gridManager.grid().leafGridView());

    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);

    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x(gridGeometry->numDofs());
    problem->applyInitialSolution(x);

    // a converged state would make the comparison degenerate, so perturb into a genuinely nonlinear one
    for (std::size_t i = 0; i < x.size(); ++i)
    {
        x[i][0] += 1e4*std::sin(3.0*double(i));
        x[i][1] += 0.1*std::sin(7.0*double(i) + 1.0);
    }

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    using Assembler = FVAssembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables);
    assembler->setLinearSystem();
    assembler->assembleJacobianAndResidual(x);
    const auto& globalJacobian = assembler->jacobian();
    const auto& globalResidual = assembler->residual();

    using Index = NonlinearPreconditioning::Partition::Index;
    const auto numElements = gridGeometry->gridView().size(0);

    // the influence map the assembler itself uses: which rows a visited element writes into
    std::vector<std::vector<Index>> influence(numElements);
    for (Index i = 0; i < static_cast<Index>(numElements); ++i)
        for (const auto& dataJ : gridGeometry->connectivityMap()[i])
            influence[i].push_back(dataJ.globalJ);

    // split along x into two subdomains
    std::vector<NonlinearPreconditioning::SubdomainIndex> cores(numElements);
    for (const auto& element : elements(gridGeometry->gridView()))
    {
        const auto eIdx = gridGeometry->elementMapper().index(element);
        cores[eIdx] = element.geometry().center()[0] < 3.0 ? 0 : 1;
    }

    for (std::size_t overlap : {0ul, 1ul, 2ul})
    {
        auto partition = std::make_shared<NonlinearPreconditioning::Partition>(influence, cores, overlap);
        const auto report = NonlinearPreconditioning::PartitionValidator::validate(*partition, influence);
        check(report.valid(), "partition is structurally valid at overlap " + std::to_string(overlap));

        for (NonlinearPreconditioning::SubdomainIndex s = 0; s < partition->numSubdomains(); ++s)
        {
            const std::string tag = " (overlap " + std::to_string(overlap)
                                    + ", subdomain " + std::to_string(s) + ")";

            // the square block: loop N_i only
            {
                NonlinearPreconditioning::SubdomainAssembler<TypeTag, DiffMethod::numeric> sub(
                    problem, gridGeometry, gridVariables, partition, influence, s, false);
                sub.assembleJacobianAndResidual(x);

                check(sub.jacobian().numColumnSinkWrites() == 0, "no column-sink writes" + tag);

                std::size_t compared = 0;
                for (const auto row : partition->dofs(s))
                {
                    const auto lr = partition->localIndex(s, row);
                    check((sub.residual()[row] - globalResidual[row]).infinity_norm() < 1e-20,
                          "residual row matches the monolithic assembly" + tag);

                    for (const auto col : partition->dofs(s))
                    {
                        if (!globalJacobian.exists(row, col))
                            continue;
                        const auto lc = partition->localIndex(s, col);
                        check(sub.jacobian().matrix().exists(lr, lc), "square entry present" + tag);
                        if (sub.jacobian().matrix().exists(lr, lc))
                        {
                            auto diff = sub.jacobian().matrix()[lr][lc];
                            diff -= globalJacobian[row][col];
                            check(diff.infinity_norm() < 1e-20, "square entry matches exactly" + tag);
                            ++compared;
                        }
                    }
                }
                check(compared > 0, "square comparison was not vacuous" + tag);
            }

            // the rectangular block: loop N_i union fringe, opening the global columns
            {
                const auto columns = partition->dofsWithFringe(s);
                std::vector<Index> colOf(partition->numDofs(), NonlinearPreconditioning::Partition::noIndex);
                for (Index k = 0; k < columns.size(); ++k)
                    colOf[columns[k]] = k;

                NonlinearPreconditioning::SubdomainAssembler<TypeTag, DiffMethod::numeric> sub(
                    problem, gridGeometry, gridVariables, partition, influence, s, true);
                sub.assembleJacobianAndResidual(x);

                check(sub.jacobian().numColumnSinkWrites() == 0, "no column-sink writes, fringe" + tag);

                std::size_t comparedFringe = 0;
                for (const auto row : partition->dofs(s))
                {
                    const auto lr = partition->localIndex(s, row);
                    for (const auto col : columns)
                    {
                        if (!globalJacobian.exists(row, col))
                            continue;
                        const auto lc = colOf[col];
                        check(sub.jacobian().matrix().exists(lr, lc), "rectangular entry present" + tag);
                        if (sub.jacobian().matrix().exists(lr, lc))
                        {
                            auto diff = sub.jacobian().matrix()[lr][lc];
                            diff -= globalJacobian[row][col];
                            check(diff.infinity_norm() < 1e-20, "rectangular entry matches exactly" + tag);
                            if (partition->localIndex(s, col) == NonlinearPreconditioning::Partition::noIndex)
                                ++comparedFringe;
                        }
                    }
                }
                if (!partition->fringe(s).empty())
                    check(comparedFringe > 0, "fringe columns were actually exercised" + tag);
            }
        }
    }

    if (failures > 0)
    {
        std::cerr << failures << " check(s) failed" << std::endl;
        return 1;
    }

    std::cout << "Restricted subdomain assembly reproduces the monolithic assembly exactly" << std::endl;
    return 0;
}
