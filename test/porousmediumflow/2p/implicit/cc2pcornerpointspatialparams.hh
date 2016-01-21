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
SET_TYPE_PROP(CC2PCornerPointSpatialParams, SpatialParams, Dumux::CC2PCornerPointSpatialParams<TypeTag>);

// Set the material Law
SET_PROP(CC2PCornerPointSpatialParams, MaterialLaw)
{
private:
    // define the material law which is parameterized by effective
    // saturations
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef RegularizedVanGenuchten<Scalar> EffectiveLaw;
public:
    // define the material law parameterized by absolute saturations
    typedef EffToAbsLaw<EffectiveLaw> type;
};
}

template<class TypeTag>
class CC2PCornerPointSpatialParams : public ImplicitSpatialParams<TypeTag>
{
    typedef ImplicitSpatialParams<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, GridCreator) GridCreator;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename Grid::ctype CoordScalar;

    enum {
        dim=GridView::dimension,
        dimWorld=GridView::dimensionworld
    };

    typedef Dune::FieldVector<CoordScalar,dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar,dimWorld,dimWorld> DimWorldMatrix;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;

public:
    //get the material law from the property system
    typedef typename GET_PROP_TYPE(TypeTag, MaterialLaw) MaterialLaw;
    typedef typename MaterialLaw::Params MaterialLawParams;

    /*!
     * \brief The constructor
     *
     * \param gridView The grid view
     */
    CC2PCornerPointSpatialParams(const GridView& gridView)
    : ParentType(gridView)
    {
        homogeneous_ = GET_RUNTIME_PARAM_FROM_GROUP(TypeTag, bool, Problem, Homogeneous);

        const std::vector<int>& globalCell = GridCreator::grid().globalCell();

        if (GridCreator::deck()->hasKeyword("PORO")) {
            std::cout << "Found PORO..." << std::endl;
            std::vector<double> eclVector = GridCreator::deck()->getKeyword("PORO")->getRawDoubleData();
            porosity_.resize(globalCell.size());

            for (size_t i = 0; i < globalCell.size(); ++i) {
                if (homogeneous_)
                    porosity_[i] = 0.2;
                else
                    porosity_[i] = eclVector[globalCell[i]];
            }
        }

        if (GridCreator::deck()->hasKeyword("PERMX")) {
            std::cout << "Found PERMX..." << std::endl;
            std::vector<double> eclVector = GridCreator::deck()->getKeyword("PERMX")->getRawDoubleData();
            permX_.resize(globalCell.size());

            for (size_t i = 0; i < globalCell.size(); ++i) {
                // assume that values are given in mD = 9.86923e-16 m^2
                if (homogeneous_)
                    permX_[i] = 9.86923e-16;
                else
                    permX_[i] = 9.86923e-16*eclVector[globalCell[i]];
            }
        }

        if (GridCreator::deck()->hasKeyword("PERMZ")) {
            std::cout << "Found PERMZ..." << std::endl;
            std::vector<double> eclVector = GridCreator::deck()->getKeyword("PERMZ")->getRawDoubleData();
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

    /*!
     * \brief Returns the scalar intrinsic permeability \f$[m^2]\f$
     *
     * \param element The finite element
     * \param fvGeometry The finite volume geometry of the element
     * \param scvIdx The local index of the sub-control volume
     */
    const DimWorldMatrix intrinsicPermeability(const Element &element,
                                       const FVElementGeometry &fvGeometry,
                                       int scvIdx) const
    {
        int eIdx = GridCreator::grid().leafGridView().indexSet().index(element);

        DimWorldMatrix K(0);
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
    Scalar porosity(const Element &element,
                    const FVElementGeometry &fvGeometry,
                    int scvIdx) const
    {
        int eIdx = GridCreator::grid().leafGridView().indexSet().index(element);

        return porosity_[eIdx];
    }

    /*!
     * \brief Returns the parameter object for the Brooks-Corey material law
     *
     * \param element The finite element
     * \param fvGeometry The finite volume geometry of the element
     * \param scvIdx The local index of the sub-control volume
     */
    const MaterialLawParams& materialLawParams(const Element &element,
                                                const FVElementGeometry &fvGeometry,
                                                int scvIdx) const
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

