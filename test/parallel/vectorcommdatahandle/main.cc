// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//

#include <config.h>

#include <array>
#include <bitset>
#include <cmath>
#include <iostream>

#include <dune/common/exceptions.hh>
#include <dune/common/parallel/mpihelper.hh>

#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/yaspgrid.hh>

#include <dumux/common/initialize.hh>
#include <dumux/parallel/vectorcommdatahandle.hh>

int main(int argc, char** argv)
{
    Dumux::initialize(argc, argv);
    const auto& mpiHelper = Dune::MPIHelper::instance();

    constexpr int dim = 2;
    using Grid = Dune::YaspGrid<dim>;
    using GridView = Grid::LeafGridView;

    const Dune::FieldVector<double, dim> upperRight(1.0);
    const std::array<int, dim> cells = {8, 8};
    auto grid = std::make_shared<Grid>(upperRight, cells);
    grid->globalRefine(1);

    const auto gridView = grid->leafGridView();

    const Dune::MCMGLayout miniLikeLayout = [](Dune::GeometryType gt, int dimGrid)
    {
        const int codim = dimGrid - gt.dim();
        return (codim == 0 || codim == dimGrid) ? std::size_t(1) : std::size_t(0);
    };

    using Mapper = Dune::MultipleCodimMultipleGeomTypeMapper<GridView>;
    Mapper mapper(gridView, miniLikeLayout);

    using Vector = std::vector<double>;
    Vector vSequential(mapper.size());
    Vector vMulti(mapper.size());

    for (std::size_t i = 0; i < mapper.size(); ++i)
    {
        const double value = static_cast<double>((i % 17) + 1);
        vSequential[i] = value;
        vMulti[i] = value;
    }

    if (mpiHelper.size() > 1)
    {
        Dumux::VectorCommDataHandleSum<Mapper, Vector, 0, double> codim0Handle(mapper, vSequential);
        Dumux::VectorCommDataHandleSum<Mapper, Vector, dim, double> codimDimHandle(mapper, vSequential);

        gridView.communicate(codim0Handle, Dune::InteriorBorder_InteriorBorder_Interface, Dune::ForwardCommunication);
        gridView.communicate(codimDimHandle, Dune::InteriorBorder_InteriorBorder_Interface, Dune::ForwardCommunication);

        std::bitset<dim + 1> activeCodims;
        activeCodims.set(0);
        activeCodims.set(dim);
        Dumux::MultiCodimVectorCommDataHandleSum<Mapper, Vector, dim + 1, double> multiHandle(mapper, vMulti, activeCodims);
        gridView.communicate(multiHandle, Dune::InteriorBorder_InteriorBorder_Interface, Dune::ForwardCommunication);
    }

    double maxDiff = 0.0;
    for (std::size_t i = 0; i < mapper.size(); ++i)
        maxDiff = std::max(maxDiff, std::abs(vSequential[i] - vMulti[i]));

    const auto globalMaxDiff = mpiHelper.getCommunication().max(maxDiff);
    if (globalMaxDiff > 1e-14)
        DUNE_THROW(Dune::Exception, "Multi-codim vector communication mismatch: maxDiff=" << globalMaxDiff);

    if (mpiHelper.rank() == 0)
        std::cout << "[test_parallel_vectorcommdatahandle_multicodim] passed with maxDiff=" << globalMaxDiff << std::endl;

    return 0;
}
