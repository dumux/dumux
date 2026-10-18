// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Collapse of a thin elastic ring under external pressure with self-contact of its inner
 *        surface
 *
 * The ring buckles into the elliptic mode at \f$ p = 3D/R^3 \f$ (bending stiffness
 * \f$ D = E t^3/(12(1-\nu^2)) \f$), its opposite walls touch at a point and then along a growing
 * line. The pressure is raised in increments, halved when Newton's method fails. Reported for
 * each pressure \f$ p R^3/D \f$ listed in Problem.Pressures: the semi-minor and semi-major axes of
 * the midline, the area enclosed by the midline over \f$ \pi R^2 \f$, the contact force between
 * the two halves \f$ Q R^2/D \f$, the position along the wall of the largest contact pressure
 * over \f$ R \f$, and the number of active contact pairs.
 */
#include <config.h>

#include <algorithm>
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numbers>
#include <vector>

#include <dune/grid/common/gridfactory.hh>

#include <dumux/common/exceptions.hh>
#include <dumux/common/initialize.hh>
#include <dumux/common/parameters.hh>
#include <dumux/common/properties.hh>
#include <dumux/assembly/assembler.hh>
#include <dumux/io/vtkoutputmodule.hh>
#include <dumux/linear/istlsolvers.hh>
#include <dumux/linear/linearalgebratraits.hh>
#include <dumux/linear/linearsolvertraits.hh>
#include <dumux/nonlinear/newtonsolver.hh>
#include <dumux/solidmechanics/contact/contactassembler.hh>
#include <dumux/solidmechanics/contact/penaltycontact.hh>
#include <dumux/solidmechanics/contact/surface.hh>

#include "properties.hh"

namespace Dumux {

//! A ring of quadrilaterals, cells elements around and layers through the wall
template<class Grid, class Scalar>
std::shared_ptr<Grid> makeRingGrid(Scalar radius, Scalar thickness, unsigned int cells, unsigned int layers)
{
    Dune::GridFactory<Grid> factory;
    for (unsigned int i = 0; i < cells; ++i)
    {
        const Scalar phi = 2.0*std::numbers::pi*i/cells;
        for (unsigned int j = 0; j <= layers; ++j)
        {
            const Scalar r = radius - 0.5*thickness + thickness*j/layers;
            factory.insertVertex({r*std::cos(phi), r*std::sin(phi)});
        }
    }
    const auto vertex = [&](unsigned int i, unsigned int j){ return (i % cells)*(layers + 1) + j; };
    for (unsigned int i = 0; i < cells; ++i)
        for (unsigned int j = 0; j < layers; ++j)
            factory.insertElement(Dune::GeometryTypes::quadrilateral,
                                  {vertex(i, j), vertex(i, j + 1), vertex(i + 1, j), vertex(i + 1, j + 1)});
    return std::shared_ptr<Grid>(factory.createGrid());
}

} // end namespace Dumux

int main(int argc, char** argv)
{
    using namespace Dumux;
    using TypeTag = Properties::TTag::RingCollapse;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    initialize(argc, argv);
    Parameters::init(argc, argv);

    const auto radius = getParam<Scalar>("Problem.Radius");
    const auto thickness = getParam<Scalar>("Problem.Thickness");
    const auto cells = getParam<unsigned int>("Grid.Cells");
    const auto layers = getParam<unsigned int>("Grid.Layers");
    if (cells % 4 != 0)
        DUNE_THROW(ParameterException, "Grid.Cells has to be a multiple of four to put nodes on the axes");

    using Grid = GetPropType<TypeTag, Properties::Grid>;
    const auto grid = makeRingGrid<Grid>(radius, thickness, cells, layers);
    const auto& leafGridView = grid->leafGridView();

    using GridGeometry = GetPropType<TypeTag, Properties::GridGeometry>;
    auto gridGeometry = std::make_shared<GridGeometry>(leafGridView);

    using Problem = GetPropType<TypeTag, Properties::Problem>;
    auto problem = std::make_shared<Problem>(gridGeometry);
    const Scalar stiffness = problem->bendingStiffness();

    using SolutionVector = GetPropType<TypeTag, Properties::SolutionVector>;
    SolutionVector x(gridGeometry->numDofs());
    x = 0.0;

    using GridVariables = GetPropType<TypeTag, Properties::GridVariables>;
    auto gridVariables = std::make_shared<GridVariables>(problem, gridGeometry);
    gridVariables->init(x);

    using Assembler = Experimental::Assembler<TypeTag, DiffMethod::numeric>;
    auto assembler = std::make_shared<Assembler>(problem, gridGeometry, gridVariables);

    using Surface = Contact::ContactSurface<typename GridGeometry::GridView>;
    auto surface = std::make_shared<Surface>(leafGridView, gridGeometry->dofMapper(),
                                             [&](const auto& is){ return is.geometry().center().two_norm() < radius; });
    using ContactModel = Contact::PenaltyContact<Surface>;
    auto contact = std::make_shared<ContactModel>(surface);
    using SystemAssembler = Contact::ContactAssembler<Assembler, ContactModel>;
    auto systemAssembler = std::make_shared<SystemAssembler>(assembler, contact, getParam<bool>("Contact.CouplePairs", true));

    using LinearSolver = UMFPackIstlSolver<SeqLinearSolverTraits, LinearAlgebraTraitsFromAssembler<SystemAssembler>>;
    auto linearSolver = std::make_shared<LinearSolver>();
    NewtonSolver<SystemAssembler, LinearSolver> nonlinearSolver(systemAssembler, linearSolver);

    // inner and outer node of the wall at each angle, for the midline
    std::vector<std::size_t> innerDof(cells), outerDof(cells);
    std::vector<Dune::FieldVector<Scalar, 2>> referencePosition(gridGeometry->numDofs());
    {
        auto elemDisc = localView(*gridGeometry);
        for (const auto& element : elements(leafGridView))
        {
            elemDisc.bind(element);
            for (const auto& localDof : localDofs(elemDisc))
            {
                const auto pos = ipData(elemDisc, localDof).global();
                referencePosition[localDof.dofIndex()] = pos;
                const Scalar phi = std::atan2(pos[1], pos[0]);
                const auto i = static_cast<std::size_t>(std::lround(phi/(2.0*std::numbers::pi)*cells) + cells) % cells;
                if (pos.two_norm() < radius - 0.49*thickness)
                    innerDof[i] = localDof.dofIndex();
                else if (pos.two_norm() > radius + 0.49*thickness)
                    outerDof[i] = localDof.dofIndex();
            }
        }
    }
    const auto midline = [&](std::size_t i)
    {
        auto p = referencePosition[innerDof[i]];
        p += x[innerDof[i]];
        p += referencePosition[outerDof[i]];
        p += x[outerDof[i]];
        p *= 0.5;
        return p;
    };

    std::vector<Scalar> contactPressure(gridGeometry->numDofs(), 0.0);
    VtkOutputModule<GridVariables, SolutionVector> vtkWriter(*gridVariables, x, problem->name());
    vtkWriter.addField(x, "d");
    vtkWriter.addField(contactPressure, "contactPressure");
    vtkWriter.write(0.0);

    std::cout << "pressure semiMinorAxis semiMajorAxis area contactForce peakPosition activePairs" << std::endl;
    const Scalar maxIncrement = getParam<Scalar>("Problem.PressureIncrement");
    Scalar current = 0.0;
    auto xOld = x;
    for (const auto target : getParam<std::vector<Scalar>>("Problem.Pressures"))
    {
        Scalar increment = maxIncrement;
        while (current < target - 1e-12)
        {
            const Scalar next = std::min(current + increment, target);
            problem->setPressure(next*stiffness/(radius*radius*radius));
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
                    DUNE_THROW(NumericalProblem, "Pressure increment too small");
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
            if (surface->referencePosition(s)[1] > 0.0)
                force -= r[surface->dofIndex(s)][1];

        std::fill(contactPressure.begin(), contactPressure.end(), 0.0);
        const auto pairs = contact->activePairs(x);
        const auto positions = surface->positions(x);
        Scalar peakPressure = 0.0, peakPosition = 0.0;
        for (const auto& p : pairs)
        {
            contactPressure[surface->dofIndex(p.node)] = p.pressure;
            const auto& pos = positions[p.node];
            if (pos[0] >= 0.0 && pos[1] > 0.0 && p.pressure > peakPressure)
            {
                peakPressure = p.pressure;
                peakPosition = pos[0];
            }
        }

        Scalar area = 0.0;
        for (std::size_t i = 0; i < cells; ++i)
        {
            const auto a = midline(i), b = midline((i + 1) % cells);
            area += 0.5*(a[0]*b[1] - a[1]*b[0]);
        }

        std::cout << std::setprecision(8) << target << " " << midline(cells/4)[1]/radius << " " << midline(0)[0]/radius << " "
                  << area/(std::numbers::pi*radius*radius) << " " << force*radius*radius/stiffness << " "
                  << peakPosition/radius << " " << pairs.size() << std::endl;
        vtkWriter.write(target);
    }

    return 0;
}
