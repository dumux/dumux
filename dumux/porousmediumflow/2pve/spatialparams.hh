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

#include <cmath>
#include <memory>
#include <type_traits>
#include <utility>
#include <vector>

#include <dune/common/exceptions.hh>
#include <dune/grid/common/rangegenerators.hh>

#include <dumux/porousmediumflow/fvspatialparamsmp.hh>
#include <dumux/porousmediumflow/2pve/columnmapping.hh>

namespace Dumux {

/*!
 * \ingroup TwoPVEModel
 * \brief The base class for the coarse-level spatial parameters of the two-phase VE model
 *
 * The vertical integrals of fine-level permeability and porosity are stored per horizontal
 * coarse element. The volume variables use column averages and the column height as extrusion factor.
 * The fine-level spatial parameters have to provide `permeabilityAtElement(fineElement)` and
 * `porosityAtElement(fineElement)`, and the porosity has to be uniform within each column.
 * The implementation has to provide
 * `fluidMatrixInteractionAtPos(globalPos)` returning a Brooks-Corey material law and
 * `temperatureAtPos(globalPos)`.
 *
 * \tparam GridGeometry the coarse-level grid geometry
 * \tparam FineGridGeometry the full-dimensional reconstruction grid geometry
 * \tparam Scalar the scalar type
 * \tparam SpatialParamsFine the fine-level spatial parameters
 * \tparam Implementation the class deriving from this class
 */
template<class GridGeometry, class FineGridGeometry, class Scalar, class SpatialParamsFine, class Implementation>
class TwoPVESpatialParams
: public FVPorousMediumFlowSpatialParamsMP<GridGeometry, Scalar, Implementation>
{
    using ParentType = FVPorousMediumFlowSpatialParamsMP<GridGeometry, Scalar, Implementation>;
    using GridView = typename GridGeometry::GridView;
    using Element = typename GridView::template Codim<0>::Entity;
    using SubControlVolume = typename GridGeometry::LocalView::SubControlVolume;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;
    using ColumnMapping = TwoPVEColumnMapping<GridGeometry, FineGridGeometry, Scalar>;

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
    , permeabilityIntegral_(gridGeometry->elementMapper().size(), 0.0)
    , porosityIntegral_(gridGeometry->elementMapper().size(), 0.0)
    , columnHeight_(gridGeometry->elementMapper().size(), 0.0)
    {
        using FinePermeability = std::decay_t<decltype(spatialParamsFine_->permeabilityAtElement(std::declval<const typename FineGridGeometry::GridView::template Codim<0>::Entity&>()))>;
        static_assert(std::is_convertible_v<FinePermeability, Scalar>, "The VE upscaling requires a scalar fine-level permeability");

        for (const auto& element : elements(gridGeometry->gridView()))
        {
            const auto columnIdx = gridGeometry->elementMapper().index(element);
            const auto& column = columnMapping.column(columnIdx);
            const Scalar columnPorosity = spatialParamsFine_->porosityAtElement(column.front());
            for (const auto& fineElement : column)
            {
                using std::abs;
                const Scalar porosityFine = spatialParamsFine_->porosityAtElement(fineElement);
                // the gas plume distance follows from the water balance of the unweighted saturation profile
                if (abs(porosityFine - columnPorosity) > 1e-10*abs(columnPorosity))
                    DUNE_THROW(Dune::InvalidStateException, "The VE model requires a uniform porosity within each column, column " << columnIdx << " has porosities " << columnPorosity << " and " << porosityFine);

                permeabilityIntegral_[columnIdx] += spatialParamsFine_->permeabilityAtElement(fineElement)*fineCellHeight;
                porosityIntegral_[columnIdx] += porosityFine*fineCellHeight;
            }

            const Scalar columnHeight = column.size()*fineCellHeight;
            columnHeight_[columnIdx] = columnHeight;
        }
    }

    /*!
     * \brief Returns the fine-level spatial parameters
     */
    const SpatialParamsFine& spatialParamsFine() const
    { return *spatialParamsFine_; }

    /*!
     * \brief Returns the column-average intrinsic permeability \f$\mathrm{[m^2]}\f$
     *
     * \param element coarse-level element
     * \param scv     sub-control volume of the element
     * \param elemSol element solution
     */
    template<class ElementSolution>
    PermeabilityType permeability(const Element& element,
                                  const SubControlVolume& scv,
                                  const ElementSolution& elemSol) const
    { return permeabilityIntegralAtElement(element)/columnHeightAtElement(element); }

    /*!
     * \brief Returns the dimensionless column-average porosity
     *
     * \param element coarse-level element
     * \param scv     sub-control volume of the element
     * \param elemSol element solution
     */
    template<class ElementSolution>
    Scalar porosity(const Element& element,
                    const SubControlVolume& scv,
                    const ElementSolution& elemSol) const
    { return porosityIntegralAtElement(element)/columnHeightAtElement(element); }

    //! Column height [m], used to extrude the horizontal cell and its faces.
    template<class ElementSolution>
    Scalar extrusionFactor(const Element& element,
                           const SubControlVolume&,
                           const ElementSolution&) const
    { return columnHeightAtElement(element); }

    /*!
     * \brief Returns the index of the wetting phase, which the VE model requires to be the first phase
     *
     * \param globalPos the global position
     */
    template<class FluidSystem>
    int wettingPhaseAtPos(const GlobalPosition& globalPos) const
    { return FluidSystem::phase0Idx; }

    //! Vertical integral of the intrinsic permeability [m^3].
    Scalar permeabilityIntegralAtElement(const Element& element) const
    { return permeabilityIntegral_[this->gridGeometry().elementMapper().index(element)]; }

    //! Vertical integral of porosity [m], used for storage per horizontal measure.
    Scalar porosityIntegralAtElement(const Element& element) const
    { return porosityIntegral_[this->gridGeometry().elementMapper().index(element)]; }

    //! Height of the reconstructed column [m].
    Scalar columnHeightAtElement(const Element& element) const
    { return columnHeight_[this->gridGeometry().elementMapper().index(element)]; }

    //! Horizontal gravity is zero; vertical gravity enters through fine-level reconstruction.
    GlobalPosition gravity(const GlobalPosition&) const
    { return GlobalPosition(0.0); }

private:
    std::shared_ptr<const SpatialParamsFine> spatialParamsFine_;
    std::vector<PermeabilityType> permeabilityIntegral_;
    std::vector<Scalar> porosityIntegral_;
    std::vector<Scalar> columnHeight_;
};

} // end namespace Dumux

#endif
