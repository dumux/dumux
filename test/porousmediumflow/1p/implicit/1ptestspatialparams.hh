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
 *
 * \brief The spatial parameters class for the test problem using the
 *        1p box model
 */
#ifndef DUMUX_1P_TEST_SPATIALPARAMS_HH
#define DUMUX_1P_TEST_SPATIALPARAMS_HH

#include <dumux/material/spatialparams/implicit1p.hh>
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

// Set properties of the porous medium
NEW_PROP_TAG(SpatialParamsRandomField);
SET_BOOL_PROP(OnePTestSpatialParams, SpatialParamsRandomField, false);
}

/*!
 * \ingroup OnePModel
 * \ingroup ImplicitTestProblems
 *
 * \brief The spatial parameters class for the test problem using the
 *        1p box model
 */
template<class TypeTag>
class OnePTestSpatialParams : public ImplicitSpatialParamsOneP<TypeTag>
{
    using ParentType = ImplicitSpatialParamsOneP<TypeTag>;
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using FVElementGeometry = typename GET_PROP_TYPE(TypeTag, FVElementGeometry);
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    using VolumeVariables = typename GET_PROP_TYPE(TypeTag, VolumeVariables);
    using IndexSet = typename GridView::IndexSet;
    using ScalarVector = std::vector<Scalar>;

    enum {
        dim=GridView::dimension,
        dimWorld=GridView::dimensionworld
    };

    using GlobalPosition = Dune::FieldVector<Scalar,dimWorld>;
    using Element = typename GridView::template Codim<0>::Entity;

public:
    OnePTestSpatialParams(const Problem& problem, const GridView& gridView)
        : ParentType(problem, gridView),
          randomPermeability_(gridView.size(dim), 0.0),
          indexSet_(gridView.indexSet())
    {
        randomField_ = GET_PARAM_FROM_GROUP(TypeTag, bool, SpatialParams, RandomField);
        permeability_ = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.Permeability);
        if(!randomField_)
            permeabilityLens_ = GET_RUNTIME_PARAM(TypeTag, Scalar, SpatialParams.PermeabilityLens);
        else
            initRandomField(gridView);

        lensLowerLeft_ = GET_RUNTIME_PARAM(TypeTag, GlobalPosition, SpatialParams.LensLowerLeft);
        lensUpperRight_ = GET_RUNTIME_PARAM(TypeTag, GlobalPosition, SpatialParams.LensUpperRight);
    }

    /*!
     * \brief Return the intrinsic permeability for the current sub-control volume in [m^2].
     *
     * \param element The current finite element
     * \param fvGeometry The current finite volume geometry of the element
     * \param scvIdx The index sub-control volume face where the
     *                      intrinsic velocity ought to be calculated.
     */
    Scalar permeability(const Element& element,
                                 const SubControlVolume& scv,
                                 const ElementSolutionVector& elemSol) const
    {
        if (isInLens_(scv.dofPosition()))
        {
            if(randomField_)
                return randomPermeability_[scv.dofIndex()];
            else
                return permeabilityLens_;
        }
        else
            return permeability_;
    }

    /*! \brief Define the porosity in [-].
   *
   * \param element The finite element
   * \param fvGeometry The finite volume geometry
   * \param scvIdx The local index of the sub-control volume where
   */
    Scalar porosityAtPos(const GlobalPosition& globalPos) const
    { return 0.4; }

    /*!
     * \brief This method allows the generation of a statistical field using gstat
     *
     * \param gridView The GridView used by the problem
     */
    void initRandomField(const GridView& gridView)
    {
        std::string gStatControlFile = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Gstat, ControlFile);
        std::string gStatInputFile = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Gstat, InputFile);
        std::string outputFilePrefix = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, std::string, Gstat, OutputFilePrefix);

        // create random permeability object
        GstatRandomField<GridView, Scalar> randomPermeabilityField(gridView);
        randomPermeabilityField.create(gStatControlFile,
                                       gStatInputFile,
                                       outputFilePrefix + ".dat",
                                       GstatRandomField<GridView, Scalar>::FieldType::log10,
                                       true);
        randomPermeability_.resize(gridView.size(dim), 0.0);

        // copy vector from the temporary gstat object
        for (const auto& element : elements(gridView))
        {
            auto index = indexSet_.index(element.template subEntity<dim> (/*scvIdx=*/0));
            randomPermeability_[index] = randomPermeabilityField.data(element);
        }

        // output the random field to vtk
        randomPermeabilityField.writeVtk(outputFilePrefix, "absolute permeability");
    }

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
    ScalarVector randomPermeability_;

    const IndexSet& indexSet_;
    static constexpr Scalar eps_ = 1.5e-7;
};

} // end namespace

#endif
