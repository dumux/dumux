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
 * \brief The spatial parameters for the El1P2CProblem which uses the
 *        linear elastic one-phase two-component model
 */
#ifndef DUMUX_ELONEPTWOCSPARAMETERS_HH
#define DUMUX_ELONEPTWOCSPARAMETERS_HH

#include <dumux/material/spatialparams/implicitspatialparams1p.hh>
#include <dumux/material/fluidmatrixinteractions/2p/linearmaterial.hh>
#include <dumux/material/fluidmatrixinteractions/2p/regularizedbrookscorey.hh>
#include <dumux/material/fluidmatrixinteractions/2p/efftoabslaw.hh>

namespace Dumux
{

//forward declaration
template<class TypeTag>
class El1P2CSpatialParams;

namespace Properties
{
// The spatial parameters TypeTag
NEW_TYPE_TAG(El1P2CSpatialParams);

// Set the spatial parameters
SET_TYPE_PROP(El1P2CSpatialParams, SpatialParams, Dumux::El1P2CSpatialParams<TypeTag>);

}

/*!
 * \ingroup ElOnePTwoCBoxModel
 * \brief The spatial parameters for the El1P2CProblem which uses the
 *        linear elastic one-phase two-component model
 */
template<class TypeTag>
class El1P2CSpatialParams : public ImplicitSpatialParamsOneP<TypeTag>
{
    typedef ImplicitSpatialParamsOneP<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, Grid) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename Grid::ctype CoordScalar;
    enum {
        dim=GridView::dimension,
        dimWorld=GridView::dimensionworld,
    };

    typedef Dune::FieldVector<CoordScalar,dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<Scalar,dim,dim> DimMatrix;

    typedef typename GET_PROP_TYPE(TypeTag, FVElementGeometry) FVElementGeometry;
    typedef typename GridView::template Codim<0>::Entity Element;

public:

    El1P2CSpatialParams(const GridView &gridView)
    : ParentType(gridView)
    {
        // intrinsic permeabilities [m^2]
        K_ = Scalar(0.0);
        for (int i = 0; i < dim; i++){
            K_[i][i] = 1.E-12; //[m²]
        }

        // porosities [-]
        phi_ = 0.2;

        // rock density [kg/m^3]
        solidDensity_ = 2650.0;

        // Young's modulus [Pa]
        E_ = 6.e9;
        // Poisson's ratio [-]
        nu_ = 0.2;
        // Lame parameters [Pa]
        lambda_ = (E_ * nu_) / ((1 + nu_)*(1 - 2 * nu_));
        mu_ = E_ / (2 * (1 + nu_));
     }

    ~El1P2CSpatialParams()
    {}

    /*!
     * \brief Apply the intrinsic permeability tensor \f$[m^2]\f$ to a pressure
     *        potential gradient.
     *
     * \param element The current finite element
     * \param fvGeometry The current finite volume geometry of the element
     * \param scvIdx The local index of the sub-control volume where
     *                    the porosity needs to be defined
     */

    const DimMatrix intrinsicPermeability(const Element &element,
                                       const FVElementGeometry &fvGeometry,
                                       int scvIdx) const
    {
        return K_;
    }

    /*!
     * \brief Define the porosity \f$[-]\f$ of the soil
     *
     * \param element The finite element
     * \param fvGeometry The finite volume geometry
     * \param scvIdx The local index of the sub-control volume where
     *                    the porosity needs to be defined
     */
    double porosity(const Element &element,
                    const FVElementGeometry &fvGeometry,
                    int scvIdx) const
    {
        return phi_;
    }

    /*!
     * \brief Define the porosity \f$[-]\f$ of the soil
     *
     * \param globalPos The global position of the vertex
     */
    double porosity(const GlobalPosition& globalPos) const
    {
        return phi_;
    }

    /*!
     * \brief Define the density \f$[kg/m^3]\f$ of the rock
     *
     * \param element The finite element
     * \param scvIdx The local index of the sub-control volume where
     *                    the porosity needs to be defined
     */
    const Scalar rockDensity(const Element &element,
                                        int scvIdx) const
    {
        return solidDensity_;
    }

    /*!
     * \brief Define the density \f$[kg/m^3]\f$ of the rock
     *
     * \param globalPos The global position of the vertex
     */
    const Scalar rockDensity(const GlobalPosition &globalPos) const
    {
        return solidDensity_;
    }

    /*!
     * \brief Define the Lame parameters \f$[Pa]\f$ linear elastic rock
     *
     * \param element The finite element
     * \param fvGeometry The finite volume geometry
     * \param scvIdx The local index of the sub-control volume where
     *                    the porosity needs to be defined
     */
    const Dune::FieldVector<Scalar,2> lameParams(const Element &element,
                                           const FVElementGeometry &fvGeometry,
                                           int scvIdx) const
    {
        Dune::FieldVector<Scalar, 2> param;

        param[0] = lambda_;
        param[1] = mu_;

        return param;
    }

    /*!
     * \brief Define if the model should apply two-point approximation 
     * instead of a box approximation for the fluxes
     *
     * \param element The finite element
     * \param vertexI first point for the two-point flux approximation
     * \param vertexJ second point for the two-point flux approximation
     */
    bool useTwoPointGradient(const Element &element,
                             int vertexI,
                             int vertexJ) const
    {
        return false;
    }

    /*!
     * \brief Return dispersivity (not needed here
     *
     * \param element The finite element
     * \param fvGeometry The finite volume geometry
     * \param scvIdx The local index of the sub-control volume where
     */
    Dune::FieldVector<Scalar,dim> dispersivity(const Element &element,
                    const FVElementGeometry &fvGeometry,
                    int scvIdx) const
    {
        return Dune::FieldVector<Scalar,dim>(0);
    }

private:
    Dune::FieldMatrix<Scalar,dim,dim> K_;
    Scalar layerBottom_;
    Scalar solidDensity_;
    Scalar phi_;
    Scalar lambda_;
    Scalar mu_;
    Scalar E_;
    Scalar nu_;
    static constexpr Scalar eps_ = 3e-6;
    int episode_;

};
}
#endif
