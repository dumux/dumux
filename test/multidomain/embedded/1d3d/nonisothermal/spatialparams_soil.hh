// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef BENCHMARKS_SOIL_SPATIALPARAMS_HH
#define BENCHMARKS_SOIL_SPATIALPARAMS_HH

#include <dumux/porousmediumflow/fvspatialparams1p.hh>

namespace Dumux {

template<class GridGeometry, class Scalar>
class SoilSpatialParams
: public FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar, SoilSpatialParams<GridGeometry, Scalar>>
{
    using ThisType = SoilSpatialParams<GridGeometry, Scalar>;
    using ParentType = FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar, ThisType>;
    using GridView = typename GridGeometry::GridView;
    using Grid = typename GridView::Grid;
    using SecondaryGridView = Dune::UGGrid<3>::LeafGridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using CoordScalar = typename GridView::ctype;
    static constexpr int dimWorld = GridView::dimensionworld;
    using Tensor = Dune::FieldMatrix<CoordScalar, dimWorld, dimWorld>;

public:
    // export permeability type
    using PermeabilityType = Scalar;

    SoilSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
    : ParentType(gridGeometry)
    {}

    /*!
    * \brief Defines the intrinsic permeability \f$\mathrm{[m^2]}\f$.
    *
    * \param element The element
    * \param scv The sub control volume
    * \param elemSol The element solution vector
    */
    template<class ElementSolution>
    PermeabilityType permeability(const Element& element,
                                const SubControlVolume& scv,
                                const ElementSolution& elemSol) const
    {
        return 1.0e-20;
    }

    /*!
    * \brief Defines the porosity \f$\mathrm{[-]}\f$.
    *
    * \param element The current finite element
    * \param scv The sub control volume
    * \param elemSol The current element solution vector
    */
    template<class ElementSolution>
    Scalar porosity(const Element& element,
                    const SubControlVolume& scv,
                    const ElementSolution& elemSol) const
    {
        return 0.0;
    }
};

} // end namespace Dumux

#endif
