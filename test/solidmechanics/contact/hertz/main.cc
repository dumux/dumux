// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Hertz contact of two elastic cylinders (plane strain) by penalty contact of two bodies
 *        in one mesh
 *
 * For every approach of the platens, the force per unit length \f$ F \f$ is the contact force
 * on the upper body. Hertz's solution for two equal cylinders of radius \f$ R \f$ gives the
 * contact half-width \f$ a = \sqrt{4 F R (1 - \nu^2)/(\pi E)} \f$ and the pressure
 * \f$ p(x) = p_0 \sqrt{1 - x^2/a^2} \f$ with \f$ p_0 = 2F/(\pi a) \f$ (Johnson, Contact
 * Mechanics, 1985, section 4.3). Reported per approach: \f$ F \f$, \f$ a/R \f$, the peak nodal
 * pressure against \f$ p_0 \f$, and the relative L2 error of the nodal pressures against
 * \f$ p(x) \f$ over the surface nodes within \f$ 3a \f$ of the axis.
 */
#include <config.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numbers>
#include <vector>

#include <dumux/common/exceptions.hh>
#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/assembly/assembler.hh>
#include <dumux/io/grid/gridmanager_alu.hh>
#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/nonlinear/newtonsolver.hh>
#include <dumux/solidmechanics/contact/contactassembler.hh>
#include <dumux/solidmechanics/contact/penaltycontact.hh>
#include <dumux/solidmechanics/contact/surface.hh>

#include "properties.hh"

int main(int argc, char** argv)
{
    using namespace Dumux;
    using TypeTag = Properties::TTag::Hertz;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    initialize(argc, argv);
    Parameters::init(argc, argv);

    using Grid = GetPropType<TypeTag, Properties::Grid>;
    GridManager<Grid> gridManager;
    gridManager.init();
    const auto& leafGridView = gridManager.grid().leafGridView();

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);

    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);
    const Scalar radius = problem->radius();

    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x(gridGeometry->numDofs());
    x = 0.0;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    using Assembler = Experimental::Assembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables);

    // the curved parts of the boundary near the contact
    using Surface = Contact::ContactSurface<typename GridGeometry::GridView>;
    auto surface = std::make_shared<Surface>(leafGridView, gridGeometry->dofMapper(),
                                             [&](const auto& is){ return std::abs(is.geometry().center()[1]) < 0.5*radius; });
    using ContactModel = Contact::PenaltyContact<Surface>;
    auto contact = std::make_shared<ContactModel>(surface);
    using SystemAssembler = Contact::ContactAssembler<Assembler, ContactModel>;
    auto systemAssembler = std::make_shared<SystemAssembler>(assembler, contact, getParam<bool>("Contact.CouplePairs", true));

    using LinearSolver = UMFPackIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<SystemAssembler>>;
    auto linearSolver = std::make_shared<LinearSolver>();
    NewtonSolver<SystemAssembler, LinearSolver> nonlinearSolver(systemAssembler, linearSolver);

    // surface nodes of the upper body
    std::vector<bool> isUpper(surface->numNodes(), false);
    for (std::size_t f = 0; f < surface->numFacets(); ++f)
    {
        Scalar centerY = 0.0;
        for (const auto s : surface->facet(f))
            centerY += surface->referencePosition(s)[1]/2.0;
        if (centerY > 0.0)
            for (const auto s : surface->facet(f))
                isUpper[s] = true;
    }

    std::vector<Scalar> contactPressure(gridGeometry->numDofs(), 0.0);
    VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, x, problem->name());
    vtkWriter.addField(x, "d");
    vtkWriter.addField(contactPressure, "contactPressure");

    const Scalar E = problem->youngsModulus();
    const Scalar nu = problem->poissonRatio();
    std::cout << "approach force a/R p0 pmax/p0-1 profileError activePairs" << std::endl;
    // the platens advance in increments of at most ApproachIncrement, halved when Newton's method fails
    const Scalar maxIncrement = getParam<Scalar>("Problem.ApproachIncrement");
    Scalar current = 0.0;
    auto xOld = x;
    for (const auto approach : getParam<std::vector<Scalar>>("Problem.Approaches"))
    {
        Scalar increment = maxIncrement;
        while (current < approach - 1e-14)
        {
            const Scalar next = std::min(current + increment, approach);
            problem->setApproach(next);
            try
            {
                systemAssembler->updateCandidates(x);
                nonlinearSolver.solve(x);
                while (systemAssembler->addMissingCandidates(x))
                    nonlinearSolver.solve(x);
                if (!systemAssembler->candidatesCover(xOld, x))
                    DUNE_THROW(NumericalProblem, "Surface moved further than the contact search covers");
            }
            catch (const NumericalProblem&)
            {
                x = xOld;
                increment *= 0.5;
                if (increment < 1e-6*maxIncrement)
                    DUNE_THROW(NumericalProblem, "Approach increment too small");
                continue;
            }
            current = next;
            xOld = x;
            increment = std::min(2.0*increment, maxIncrement);
        }

        SolutionVector r(x.size());
        r = 0.0;
        contact->addResidual(r, x);
        Scalar force = 0.0;
        for (std::size_t s = 0; s < surface->numNodes(); ++s)
            if (isUpper[s])
                force -= r[surface->dofIndex(s)][1];

        const Scalar a = std::sqrt(4.0*force*radius*(1.0 - nu*nu)/(std::numbers::pi*E));
        const Scalar p0 = 2.0*force/(std::numbers::pi*a);
        const auto hertzPressure = [&](Scalar xPos){ return p0*std::sqrt(std::max(0.0, 1.0 - xPos*xPos/(a*a))); };

        std::fill(contactPressure.begin(), contactPressure.end(), 0.0);
        std::vector<Scalar> nodalPressure(surface->numNodes(), 0.0);
        const auto pairs = contact->activePairs(x);
        for (const auto& p : pairs)
        {
            nodalPressure[p.node] = p.pressure;
            contactPressure[surface->dofIndex(p.node)] = nodalPressure[p.node];
        }

        Scalar pMax = 0.0, errorSquared = 0.0, normSquared = 0.0;
        const auto positions = surface->positions(x);
        for (std::size_t s = 0; s < surface->numNodes(); ++s)
        {
            pMax = std::max(pMax, nodalPressure[s]);
            if (std::abs(positions[s][0]) < 3.0*a)
            {
                const Scalar area = surface->referenceArea(s);
                const Scalar reference = hertzPressure(positions[s][0]);
                errorSquared += area*(nodalPressure[s] - reference)*(nodalPressure[s] - reference);
                normSquared += area*reference*reference;
            }
        }

        std::cout << std::setprecision(6) << approach << " " << force << " " << a/radius << " " << p0 << " "
                  << pMax/p0 - 1.0 << " " << std::sqrt(errorSquared/normSquared) << " " << pairs.size() << std::endl;
        vtkWriter.write(approach);
    }

    return 0;
}
