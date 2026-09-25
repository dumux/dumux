// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightText: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup TwoPVEModel
 * \brief The base class for the coarse-level spatial parameters of the two-phase VE model.
 */

#ifndef DUMUX_TWOPVE_SPATIAL_PARAMS_HH
#define DUMUX_TWOPVE_SPATIAL_PARAMS_HH

#include <memory>
#include <vector>

#include <dune/grid/common/rangegenerators.hh>

#include <dumux/porousmediumflow/fvspatialparamsmp.hh>
#include <dumux/porousmediumflow/2pve/columnmapping.hh>

namespace Dumux {

/*!
 * \ingroup TwoPVEModel
 * \brief The base class for the coarse-level spatial parameters of the two-phase VE model
 *
 * The coarse-level permeability and porosity of a column are the vertical integrals of the
 * fine-level permeability and porosity over the column, in \f$\mathrm{[m^3]}\f$ and \f$\mathrm{[m]}\f$.
 * The fine-level spatial parameters have to provide `permeabilityAtElement(fineElement)` and
 * `porosityAtElement(fineElement)`. The implementation has to provide
 * `fluidMatrixInteractionAtPos(globalPos)` returning a Brooks-Corey material law and
 * `temperatureAtPos(globalPos)`.
 *
 * \tparam GridGeometry the coarse-level grid geometry
 * \tparam Scalar the scalar type
 * \tparam SpatialParamsFine the fine-level spatial parameters
 * \tparam Implementation the class deriving from this class
 */
template<class GridGeometry, class Scalar, class SpatialParamsFine, class Implementation>
class TwoPVESpatialParams
: public FVPorousMediumFlowSpatialParamsMP<GridGeometry, Scalar, Implementation>
{
    using ParentType = FVPorousMediumFlowSpatialParamsMP<GridGeometry, Scalar, Implementation>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using SubControlVolume = typename GridGeometry::LocalView::SubControlVolume;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using ColumnMapping = VEColumnMapping<GridGeometry, Scalar>;

public:
    using PermeabilityType = Scalar;

    /*!
     * \brief Computes the coarse-level permeability and porosity of all columns
     *
     * \param gridGeometry      coarse-level grid geometry
     * \param columnMapping     mapping between coarse-level columns and fine-level elements
     * \param spatialParamsFine fine-level spatial parameters
     * \param fineCellHeight    vertical height of the fine-level elements
     */
    TwoPVESpatialParams(std::shared_ptr<const GridGeometry> gridGeometry,
                        const ColumnMapping& columnMapping,
                        std::shared_ptr<const SpatialParamsFine> spatialParamsFine,
                        const Scalar fineCellHeight)
    : ParentType(gridGeometry)
    , spatialParamsFine_(spatialParamsFine)
    , permeability_(gridGeometry->elementMapper().size(), 0.0)
    , porosity_(gridGeometry->elementMapper().size(), 0.0)
    {
        for (const auto& element : elements(gridGeometry->gridView()))
        {
            const auto columnIdx = gridGeometry->elementMapper().index(element);
            for (const auto& fineElement : columnMapping.column(columnIdx))
            {
                permeability_[columnIdx] += spatialParamsFine_->permeabilityAtElement(fineElement)*fineCellHeight;
                porosity_[columnIdx] += spatialParamsFine_->porosityAtElement(fineElement)*fineCellHeight;
            }
        }
    }

    /*!
     * \brief Returns the fine-level spatial parameters
     */
    const SpatialParamsFine& spatialParamsFine() const
    { return *spatialParamsFine_; }

    /*!
     * \brief Returns the coarse-level permeability \f$\mathrm{[m^3]}\f$ of a column
     *
     * \param element coarse-level element
     * \param scv     sub-control volume of the element
     * \param elemSol element solution
     */
    template<class ElementSolution>
    PermeabilityType permeability(const Element& element,
                                  const SubControlVolume& scv,
                                  const ElementSolution& elemSol) const
    { return permeabilityAtElement(element); }

    /*!
     * \brief Returns the coarse-level permeability \f$\mathrm{[m^3]}\f$ of a column
     *
     * \param element coarse-level element
     */
    PermeabilityType permeabilityAtElement(const Element& element) const
    { return permeability_[this->gridGeometry().elementMapper().index(element)]; }

    /*!
     * \brief Returns the coarse-level porosity \f$\mathrm{[m]}\f$ of a column
     *
     * \param element coarse-level element
     * \param scv     sub-control volume of the element
     * \param elemSol element solution
     */
    template<class ElementSolution>
    Scalar porosity(const Element& element,
                    const SubControlVolume& scv,
                    const ElementSolution& elemSol) const
    { return porosityAtElement(element); }

    /*!
     * \brief Returns the coarse-level porosity \f$\mathrm{[m]}\f$ of a column
     *
     * \param element coarse-level element
     */
    Scalar porosityAtElement(const Element& element) const
    { return porosity_[this->gridGeometry().elementMapper().index(element)]; }

    /*!
     * \brief Returns the index of the wetting phase, which the VE model requires to be the first phase
     *
     * \param globalPos the global position
     */
    template<class FluidSystem>
    int wettingPhaseAtPos(const GlobalPosition& globalPos) const
    { return FluidSystem::phase0Idx; }

private:
    std::shared_ptr<const SpatialParamsFine> spatialParamsFine_;
    std::vector<PermeabilityType> permeability_;
    std::vector<Scalar> porosity_;
};

} // end namespace Dumux

#endif
