// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \brief Test for finite volume element geometry, sub control volume, and sub
          control volume faces
 */
#include <config.h>

#include <iostream>
#include <utility>

#include <dune/common/test/iteratortest.hh>
#include <dune/common/fvector.hh>
#include <dune/grid/utility/structuredgridfactory.hh>
#include <dune/grid/yaspgrid.hh>

#include <dumux/common/initialize.hh>
#include <dumux/discretization/pq2/fvgridgeometry.hh>

#ifndef DOXYGEN
namespace Dumux::Detail {
template<class T>
class NoopFunctor {
public:
  NoopFunctor() {}
  void operator()(const T& t){}
};
} // end namespace Dumux::Detail
#endif

int main (int argc, char *argv[])
{
    using namespace Dumux;

    // maybe initialize MPI and/or multithreading backend
    initialize(argc, argv);

    std::cout << "Checking the FVGeometries, SCVs and SCV faces" << std::endl;

    using Grid = Dune::YaspGrid<3>;

    constexpr int dim = Grid::dimension;
    using Scalar = typename Grid::LeafGridView::ctype;

    using GridGeometry = PQ2FVGridGeometry<Scalar, typename Grid::LeafGridView, ENABLE_CACHING>;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using SubControlVolumeFace = typename FVElementGeometry::SubControlVolumeFace;
    using GlobalPosition = typename SubControlVolume::Traits::GlobalPosition;

    // make a grid
    GlobalPosition lower(0.0);
    GlobalPosition upper(1.0);
    std::array<unsigned int, dim> els{{1, 1, 1}};
    std::shared_ptr<Grid> grid = Dune::StructuredGridFactory<Grid>::createCubeGrid(lower, upper, els);
    auto leafGridView = grid->leafGridView();

    GridGeometry gridGeometry(leafGridView);

    // iterate over elements. For every element get fv geometry and loop over scvs and scvfaces
    for (const auto& element : elements(leafGridView))
    {
        auto eIdx = gridGeometry.elementMapper().index(element);
        std::cout << std::endl << "Checking fvGeometry of element " << eIdx << " with center: " << element.geometry().center() <<std::endl;
        auto fvGeometry = localView(gridGeometry);

        // bind the local view to the element
        if (fvGeometry.isBound()) DUNE_THROW(Dune::Exception, "Local view should not be bound at this point");
        fvGeometry.bind(element);
        if (!fvGeometry.isBound()) DUNE_THROW(Dune::Exception, "Local view should be bound at this point");

        // make sure the bound element fits
        auto eIdxBound = gridGeometry.elementMapper().index(fvGeometry.element());
        if (eIdx != eIdxBound)
            DUNE_THROW(Dune::Exception, "Bound element index does not match");

        auto range = scvs(fvGeometry);
        Detail::NoopFunctor<SubControlVolume> op;
        if(0 != testForwardIterator(range.begin(), range.end(), op))
            DUNE_THROW(Dune::Exception, "Iterator does not fulfill the forward iterator concept");

        const auto& localBasis = fvGeometry.feLocalBasis();
        const auto& localCoefficients = fvGeometry.feLocalCoefficients();
        using ShapeValue = typename Dune::FieldVector<double, 1>;
        std::vector<ShapeValue> shapeValues;
        std::size_t numCvLocalDofs = 0;

        for (const auto& localDof : cvLocalDofs(fvGeometry))
        {
            ++numCvLocalDofs;
            const auto& localKey = localCoefficients.localKey(localDof.index());
            if (localKey.codim() != dim)
                DUNE_THROW(Dune::Exception, "Expected cvLocalDof to have codim == dim");

            const auto scvIdx = localKey.subEntity();
            const auto& scv = fvGeometry.scv(scvIdx);
            if (scv.localDofIndex() != localDof.index())
                DUNE_THROW(Dune::Exception, "scv local dof index does not match localDof index");

            std::cout << "-- scv " << scv.localDofIndex() << " center at: " << scv.center() << " , volume: " << scv.volume()
                      << " , dof Position at: " << scv.dofPosition() << " , dof index: " << scv.dofIndex() << " , isOverlapping: " << scv.isOverlapping() << std::endl;

            // evaluate shape functions and gradients at the integration point
            const auto ipLocal = element.geometry().local(scv.dofPosition());
            localBasis.evaluateFunction(ipLocal, shapeValues);

            std::cout << "ipLocal: " << ipLocal << std::endl;
            std::cout << "shapeValues: ";
            for(int i=0; i< shapeValues.size(); i++)
                std::cout << shapeValues[i];
            std::cout << std::endl;

            std::cout << "localKeys: ";
            for(int i=0; i< shapeValues.size(); i++)
                std::cout << localCoefficients.localKey(i);
            std::cout << std::endl;
        }

        if (numCvLocalDofs != fvGeometry.numScv())
            DUNE_THROW(Dune::Exception, "Number of cvLocalDofs does not match number of SCVs");

        std::size_t numNonCvLocalDofs = 0;
        for (const auto& localDof : nonCVLocalDofs(fvGeometry))
        {
            ++numNonCvLocalDofs;
            const auto& localKey = localCoefficients.localKey(localDof.index());
            if (localKey.codim() == dim)
                DUNE_THROW(Dune::Exception, "Expected nonCVLocalDof to have codim != dim");

            const auto data = ipData(fvGeometry, localDof);
            std::cout << "-- nonCV local dof " << localDof.index()
                      << " global dof index: " << localDof.dofIndex()
                      << " local position: " << data.local()
                      << " global position: " << data.global() << std::endl;
        }

        if (numCvLocalDofs + numNonCvLocalDofs != fvGeometry.numLocalDofs())
            DUNE_THROW(Dune::Exception, "Number of cv and non-cv local dofs does not match total number of local dofs");

        auto range2 = scvfs(fvGeometry);
        Detail::NoopFunctor<SubControlVolumeFace> op2;
        if(0 != testForwardIterator(range2.begin(), range2.end(), op2))
            DUNE_THROW(Dune::Exception, "Iterator does not fulfill the forward iterator concept");

        std::size_t boundaryCount = 0;
        for (auto&& scvf : scvfs(fvGeometry))
        {
            const auto geo = fvGeometry.geometry(scvf);
            const auto faceIpData = ipData(fvGeometry, scvf);
            std::cout << "-- scvf " << scvf.index() << " ip at: " << scvf.ipGlobal() << " normal: " << scvf.unitOuterNormal()
                      << " area: " << scvf.area() << " isOverlapping: " << (scvf.isOverlapping() ? "true" : "false")
                      << " insideScvIdx: " << scvf.insideScvIdx()
                      << " geo type: " << geo.type()
                      << " geo corners: " << geo.corners();

            if ((faceIpData.global()-scvf.ipGlobal()).two_norm() > 1e-14)
                DUNE_THROW(Dune::Exception, "ipData(fvGeometry, scvf).global() does not match scvf.ipGlobal()");

            if ((geo.center()-scvf.center()).two_norm() > 1e-14)
                DUNE_THROW(Dune::Exception, "Center of scvf-geometry and scvf do not match! " << geo.center() << ", " << scvf.center());

            if(!scvf.boundary())
                std::cout << " outsideScvIdx: " << scvf.outsideScvIdx();
            if (scvf.boundary())
            {
                ++boundaryCount;
                std::cout << " (on boundary).";
            }

            // verify that boundary faces have no neighbor
            if (scvf.boundary() && scvf.numOutsideScvs() != 0)
                DUNE_THROW(Dune::Exception, "Boundary face states that it has a neighbor");

            // verify that non-boundary faces have a single neighbor
            if (!scvf.boundary() && scvf.numOutsideScvs() != 1)
                DUNE_THROW(Dune::Exception, "Expected non-boundary face to have a single neighbor");

            std::cout << std::endl;
        }

        for (const auto& intersection : intersections(leafGridView, element))
        {
            // Local dof intersection range only properly works for boundary intersection
            if(!intersection.boundary())
                continue;

            std::size_t numLocalDofsOnIntersection = 0;
            for (const auto& localDof : localDofs(fvGeometry, intersection))
            {
                ++numLocalDofsOnIntersection;
                const auto intIpData = ipData(fvGeometry, localDof);
                const auto iCenter = intersection.geometry().center();
                if ( std::abs( (intIpData.global()-iCenter) * intersection.centerUnitOuterNormal() ) > 1e-12)
                    DUNE_THROW(Dune::Exception, "ipData(fvGeometry, localDof).global() is not on intersection.");
            }

            if (numLocalDofsOnIntersection == 0)
                DUNE_THROW(Dune::Exception, "Expected at least one local dof on each intersection");
        }

        if ((boundaryCount>0) != fvGeometry.hasBoundaryScvf())
            DUNE_THROW(Dune::InvalidStateException, "fvGeometry.hasBoundaryScvf() reports " << fvGeometry.hasBoundaryScvf()
                            << " but the number of boundary scvfs is " << boundaryCount);
    }
}
