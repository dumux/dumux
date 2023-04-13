// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
/*!
 * \file
 * \ingroup OnePTests
 * \brief The spatial parameters class for the test problem using the
 *        1p box model.
 */

#ifndef DUMUX_1P_TEST_SPATIALPARAMS_HH
#define DUMUX_1P_TEST_SPATIALPARAMS_HH

#include <dumux/porousmediumflow/properties.hh>
#include <dumux/porousmediumflow/fvspatialparams1p.hh>
#include <dumux/material/gstatrandomfield.hh>

namespace Dumux {

/*!
 * \ingroup OnePTests
 * \brief The spatial parameters class for the test problem using the
 *        1p box model.
 */
template<class GridGeometry, class Scalar>
class OnePTestSpatialParams
: public FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar,
                                         OnePTestSpatialParams<GridGeometry, Scalar>>
{
    using GridView = typename GridGeometry::GridView;
    using FVElementGeometry = typename GridGeometry::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using IndexSet = typename GridView::IndexSet;

    using ThisType = OnePTestSpatialParams<GridGeometry, Scalar>;
    using ParentType = FVPorousMediumFlowSpatialParamsOneP<GridGeometry, Scalar, ThisType>;

    enum {
        dim=GridView::dimension,
        dimWorld=GridView::dimensionworld
    };

    using Element = typename GridView::template Codim<0>::Entity;
    using GlobalPosition = typename Element::Geometry::GlobalCoordinate;

public:
    // export permeability type
    using PermeabilityType = Scalar;

    OnePTestSpatialParams(std::shared_ptr<const GridGeometry> gridGeometry)
        : ParentType(gridGeometry),
          randomPermeability_(gridGeometry->gridView().size(dim), 0.0),
          indexSet_(gridGeometry->gridView().indexSet())
    {
        randomField_ = getParam<bool>("SpatialParams.RandomField", false);
        permeability_ = getParam<Scalar>("SpatialParams.Permeability");
        if(!randomField_)
            permeabilityLens_ = getParam<Scalar>("SpatialParams.PermeabilityLens");
        else
            initRandomField(*gridGeometry);

        lensLowerLeft_ = getParam<GlobalPosition>("SpatialParams.LensLowerLeft");
        lensUpperRight_ = getParam<GlobalPosition>("SpatialParams.LensUpperRight");
    }

    /*!
     * \brief Function for defining the (intrinsic) permeability \f$[m^2]\f$.
     *
     * \param element The element
     * \param scv The sub-control volume
     * \param elemSol The element solution vector
     * \return The intrinsic permeability
     */
    template<class ElementSolution>
    PermeabilityType permeability(const Element& element,
                                  const SubControlVolume& scv,
                                  const ElementSolution& elemSol) const
    {
        if (isInLens_(scv.dofPosition()))
        {
            if(randomField_)
                return randomPermeability_[indexSet_.index(element)];
            else
                return permeabilityLens_;
        }
        else
            return permeability_;
    }

    /*! \brief Defines the porosity in [-].
   *
   * \param globalPos The global position where we evaluate
   */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return 0.4; }

    /*!
     * \brief This method allows the generation of a statistical field using gstat
     *
     * \param gg The finite-volume grid geometry used by the problem
     */
    void initRandomField(const GridGeometry& gg)
    {
        const auto& gridView = gg.gridView();
        const auto& elementMapper = gg.elementMapper();
        const auto gStatControlFile = getParam<std::string>("Gstat.ControlFile");
        const auto gStatInputFile = getParam<std::string>("Gstat.InputFile");
        const auto outputFilePrefix = getParam<std::string>("Gstat.OutputFilePrefix");

        // create random permeability object
        using RandomField = GstatRandomField<GridView, Scalar>;
        RandomField randomPermeabilityField(gridView, elementMapper);
        randomPermeabilityField.create(gStatControlFile,
                                       gStatInputFile,
                                       outputFilePrefix + ".dat",
                                       RandomField::FieldType::log10,
                                       true);
        randomPermeability_.resize(gridView.size(dim), 0.0);

        // copy vector from the temporary gstat object
        randomPermeability_ = randomPermeabilityField.data();
    }

    //! get the permeability field for output
    const std::vector<Scalar>& getPermField() const
    { return randomPermeability_; }

private:
    bool isInLens_(const GlobalPosition &globalPos) const
    {
        for (int i = 0; i < dimWorld; ++i) {
            if (globalPos[i] < lensLowerLeft_[i] + eps_ || globalPos[i] > lensUpperRight_[i] - eps_)
                return false;
        }
        return true;
    }

    bool randomField_;
    GlobalPosition lensLowerLeft_;
    GlobalPosition lensUpperRight_;

    Scalar permeability_, permeabilityLens_;
    std::vector<Scalar> randomPermeability_;

    const IndexSet& indexSet_;
    static constexpr Scalar eps_ = 1.5e-7;
};

} // end namespace

#endif
