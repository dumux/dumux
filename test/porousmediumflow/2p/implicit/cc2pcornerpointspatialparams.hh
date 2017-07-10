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
#ifndef DUMUX_CC2P_CORNERPOINT_SPATIAL_PARAMS_HH
#define DUMUX_CC2P_CORNERPOINT_SPATIAL_PARAMS_HH

#include <dumux/material/spatialparams/implicit.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedvangenuchten.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>
#include <dumux/porousmediumflow/2p/implicit/model.hh>

namespace Dumux
{

//forward declaration
template<class TypeTag>
class CC2PCornerPointSpatialParams;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(CC2PCornerPointSpatialParams);

// Set the spatial parameters
SET_TYPE_PROP(CC2PCornerPointSpatialParams, SpatialParams, CC2PCornerPointSpatialParams<TypeTag>);

// Set the material Law
SET_PROP(CC2PCornerPointSpatialParams, MaterialLaw)
{
private:
    // define the material law which is parameterized by effective
    // saturations
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);
    using EffectiveLaw = RegularizedVanGenuchten<Scalar>;
public:
    // define the material law parameterized by absolute saturations
    using type = EffToAbsLaw<EffectiveLaw>;
};
}

template<class TypeTag>
class CC2PCornerPointSpatialParams : public ImplicitSpatialParams<TypeTag>
{
    using ParentType = ImplicitSpatialParams<TypeTag>;
    using Problem = typename GET_PROP_TYPE(TypeTag, Problem);
    using GridView = typename GET_PROP_TYPE(TypeTag, GridView);
    using ElementSolutionVector = typename GET_PROP_TYPE(TypeTag, ElementSolutionVector);
    using GridCreator = typename GET_PROP_TYPE(TypeTag, GridCreator);
    using Scalar = typename GET_PROP_TYPE(TypeTag, Scalar);

    enum {
        dim = GridView::dimension,
        dimWorld = GridView::dimensionworld
    };

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
    using Tensor = Dune::FieldMatrix<Scalar, dimWorld, dimWorld>;

    using Element = typename GridView::template Codim<0>::Entity;
    using SubControlVolume = typename GET_PROP_TYPE(TypeTag, SubControlVolume);
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;

public:
    // export permeability type
    using PermeabilityType = Tensor;

    //get the material law from the property system
    using MaterialLaw = typename GET_PROP_TYPE(TypeTag, MaterialLaw);
    using MaterialLawParams = typename MaterialLaw::Params;
    /*!
     * \brief The constructor
     *
     * \param gridView The grid view
     */
    CC2PCornerPointSpatialParams(const Problem& problem, const GridView& gridView)
    : ParentType(problem, gridView)
    {
        homogeneous_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, bool, Problem, Homogeneous);

        const std::vector<int>& globalCell = GridCreator::grid().globalCell();

        if (GridCreator::deck().hasKeyword("PORO")) {
            std::cout << "Found PORO..." << std::endl;
            std::vector<double> eclVector = GridCreator::deck().getKeyword("PORO").getRawDoubleData();
            porosity_.resize(globalCell.size());

            for (size_t i = 0; i < globalCell.size(); ++i) {
                if (homogeneous_)
                    porosity_[i] = 0.2;
                else
                    porosity_[i] = eclVector[globalCell[i]];
            }
        }

        if (GridCreator::deck().hasKeyword("PERMX")) {
            std::cout << "Found PERMX..." << std::endl;
            std::vector<double> eclVector = GridCreator::deck().getKeyword("PERMX").getRawDoubleData();
            permX_.resize(globalCell.size());

            for (size_t i = 0; i < globalCell.size(); ++i) {
                // assume that values are given in mD = 9.86923e-16 m^2
                if (homogeneous_)
                    permX_[i] = 9.86923e-16;
                else
                    permX_[i] = 9.86923e-16*eclVector[globalCell[i]];
            }
        }

        if (GridCreator::deck().hasKeyword("PERMZ")) {
            std::cout << "Found PERMZ..." << std::endl;
            std::vector<double> eclVector = GridCreator::deck().getKeyword("PERMZ").getRawDoubleData();
            permZ_.resize(globalCell.size());

            for (size_t i = 0; i < globalCell.size(); ++i) {
                // assume that values are given in mD = 9.86923e-16 m^2
                if (homogeneous_)
                    permZ_[i] = 9.86923e-16;
                else
                    permZ_[i] = 9.86923e-16*eclVector[globalCell[i]];
            }
        }
        else {
            std::cout << "PERMZ not found, set equal to PERMX..." << std::endl;
            permZ_ = permX_;
        }

        // parameters for the Van Genuchten law
        materialParams_.setVgAlpha(0.00045);
        materialParams_.setVgn(7.3);
    }

    PermeabilityType permeability(const Element& element,
                                      const SubControlVolume& scv,
                                      const ElementSolutionVector& elemSol) const
    {
        int eIdx = GridCreator::grid().leafGridView().indexSet().index(element);

        Tensor K(0);
        K[0][0] = K[1][1] = permX_[eIdx];
        K[2][2] = permZ_[eIdx];

        return K;
    }

    /*!
     * \brief Returns the porosity \f$[-]\f$
     *
     * \param element The finite element
     * \param fvGeometry The finite volume geometry of the element
     * \param scvIdx The local index of the sub-control volume
     */
    Scalar porosity(const Element& element,
                        const SubControlVolume& scv,
                        const ElementSolutionVector& elemSol) const
    {
        int eIdx = GridCreator::grid().leafGridView().indexSet().index(element);

        return porosity_[eIdx];
    }

    /*!
     * \brief Function for defining the parameters needed by constitutive relationships (kr-sw, pc-sw, etc.).
     *
     * \param element The current element
     * \param scv The sub-control volume inside the element.
     * \param elemSol The solution at the dofs connected to the element.
     * \return the material parameters object
     */
    const MaterialLawParams& materialLawParams(const Element& element,
                                                const SubControlVolume& scv,
                                                const ElementSolutionVector& elemSol) const
    {
        return materialParams_;
    }

private:
    MaterialLawParams materialParams_;
    std::vector<Scalar> porosity_;
    std::vector<Scalar> permX_;
    std::vector<Scalar> permZ_;
    bool homogeneous_;
};

} // end namespace
#endif

