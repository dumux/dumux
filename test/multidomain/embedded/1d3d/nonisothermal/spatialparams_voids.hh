// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
#ifndef BENCHMARKS_VOIDS_SPATIALPARAMS_HH
#define BENCHMARKS_VOIDS_SPATIALPARAMS_HH

#include <dumux/porousmediumflow/fvspatialparams1p.hh>
#include <dumux/common/parameters.hh>
#include <dumux/io/grid/griddata.hh>

namespace Dumux {

template<class GridGeometry, class Scalar>
class VoidsSpatialParams
: public FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar, VoidsSpatialParams<GridGeometry, Scalar>>
{
    using ThisType = VoidsSpatialParams<GridGeometry, Scalar>;
    using ParentType = FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar, ThisType>;
    using GlobalPosition = typename GridGeometry::GlobalCoordinate;
    using Element = typename GridGeometry::GridView::template Codim<0>::Entity;
    using SubControlVolume = typename GridGeometry::SubControlVolume;
    using Grid = typename GridGeometry::Grid;

public:
    // export permeability type
    using PermeabilityType = Scalar;

    VoidsSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry,
                      std::shared_ptr<const GridData<Grid>> gridData)
    : ParentType(gridGeometry), gridData_(gridData)

    {
        const auto& gv = gridGeometry->gridView();
        radii_.resize(gv.size(0));
        for (const auto& element : elements(gv))
        {
            const auto eIdx = gv.indexSet().index(element);
            auto level0element = element;
            for (auto levelIdx = element.level(); levelIdx != 0; levelIdx--)
                level0element = level0element.father();

            const auto& params = gridData_->parameters(level0element);
            radii_[eIdx] = params[0];
        }
    }

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
        const Scalar r = radius(this->gridGeometry().elementMapper().index(element));
        const auto priVars = elemSol[scv.localDofIndex()];
        const Scalar density = Components::SimpleH2O<Scalar>::liquidDensity(priVars[1]/*temperature*/, priVars[0]/*pressure*/);
        return r * r * density;
    }

    /*!
     * \brief Returns the radius of the circular pipe for the current sub-control volume in [m].
     *
     * \param eIdxGlobal the index of the element
     */
    Scalar radius(unsigned int eIdxGlobal) const
    {
        return radii_[eIdxGlobal];
    }

    /*!
     * \brief Returns all radii of the circular pipe in [m].
     */
    const std::vector<Scalar>& getRadii() const
    {
        return radii_;
    }

    /*!
     * \brief Returns the porosity \f$[-]\f$.
     *
     * \param globalPos the scv center
     */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    {
        return 1;
    }

    /*!
     * \brief Returns how much the domain is extruded at a given sub-control volume.
     *
     * The extrusion factor here extrudes the 1d line to a circular tube with
     * cross-section area pi*r^2.
     */
    template<class ElementSolution>
    Scalar extrusionFactor(const Element &element,
                           const SubControlVolume &scv,
                           const ElementSolution& elemSol) const
    {
        const auto eIdx = this->gridGeometry().elementMapper().index(element);
        const auto r = radius(eIdx);
        return M_PI*r*r;
    }

private:
    std::shared_ptr<const GridData<Grid>> gridData_;
    std::vector<Scalar> radii_;
};

} // end namespace Dumux

#endif
