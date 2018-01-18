// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup OnePTests
 * \brief The spatial parameters class for the test problem using the
 *        1p box model
 */
#ifndef DUMUX_1P_TEST_SPATIALPARAMS_HH
#define DUMUX_1P_TEST_SPATIALPARAMS_HH

#include <dumux/material/spatialparams/fv1p.hh>
#include <dumux/material/spatialparams/gstatrandomfield.hh>

namespace Dumux
{

//forward declaration
template<class TypeTag>
class OnePTestSpatialParams;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(OnePTestSpatialParams);
}

/*!
 * \ingroup OnePTests
 * \brief The spatial parameters class for the test problem using the
 *        1p box model
 */
template<class TypeTag>
class OnePTestSpatialParams : public FVSpatialParamsOneP<TypeTag>
{
    using ParentType = FVSpatialParamsOneP<TypeTag>;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVGridGeometry)::LocalView;
    using SubControlVolume = typename FVElementGeometry::SubControlVolume;
    using IndexSet = typename GridView::IndexSet;

    enum {
        dim=GridView::dimension,
        dimWorld=GridView::dimensionworld
    };

    using GlobalPosition = Dune::FieldVector<Scalar,dimWorld>;
    using Element = typename GridView::template Codim<0>::Entity;

public:
    // export permeability type
    using PermeabilityType = Scalar;

    OnePTestSpatialParams(const Problem& problem)
        : ParentType(problem),
          randomPermeability_(problem.fvGridGeometry().gridView().size(dim), 0.0),
          indexSet_(problem.fvGridGeometry().gridView().indexSet())
    {
        randomField_ = getParam<bool>("SpatialParams.RandomField", false);
        permeability_ = getParam<Scalar>("SpatialParams.Permeability");
        if(!randomField_)
            permeabilityLens_ = getParam<Scalar>("SpatialParams.PermeabilityLens");
        else
            initRandomField(problem.fvGridGeometry());

        lensLowerLeft_ = getParam<GlobalPosition>("SpatialParams.LensLowerLeft");
        lensUpperRight_ = getParam<GlobalPosition>("SpatialParams.LensUpperRight");
    }

    /*!
     * \brief Function for defining the (intrinsic) permeability \f$[m^2]\f$.
     *
     * \param element The element
     * \param scv The sub control volume
     * \param elemSol The element solution vector
     * \return the intrinsic permeability
     */
    PermeabilityType permeability(const Element& element,
                                  const SubControlVolume& scv,
                                  const ElementSolutionVector& elemSol) const
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

    /*! \brief Define the porosity in [-].
   *
   * \param globalPos The global position where we evaluate
   */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return 0.4; }

    /*!
     * \brief This method allows the generation of a statistical field using gstat
     *
     * \param gridView The GridView used by the problem
     */
    template<class FVGridGeometry>
    void initRandomField(const FVGridGeometry& gg)
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
