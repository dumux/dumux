// $Id$
/*****************************************************************************
 *   Copyright (C) 2010 by Markus Wolff                                      *
 *   Copyright (C) 2010 by Klaus Mosthaf                                     *
 *   Copyright (C) 2007-2008 by Bernd Flemisch                               *
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
#ifndef DUMUX_LENSSPATIALPARAMETERS_1P2C_HH
#define DUMUX_LENSSPATIALPARAMETERS_1P2C_HH

#include <dumux/material/spatialparameters/boxspatialparameters.hh>
#include <dumux/boxmodels/1p2c/1p2cmodel.hh>

/**
 * @file
 * @brief Class for defining spatial parameters
 * @author Bernd Flemisch, Klaus Mosthaf, Markus Wolff
 */

namespace Dumux
{

/** \todo Please doc me! */

template<class TypeTag>
class LensSpatialParameters1p2c : public BoxSpatialParameters<TypeTag>
{
    typedef BoxSpatialParameters<TypeTag> ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Grid)) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename Grid::ctype CoordScalar;

//    typedef typename GET_PROP_TYPE(TypeTag, PTAG(OnePTwoCIndices)) Indices;

    enum {
        dim=GridView::dimension,
        dimWorld=GridView::dimensionworld,
//
//        // indices of the primary variables
//        contiEqIdx = Indices::contiEqIdx,
//        transEqIdx = Indices::transEqIdx
    };

    typedef Dune::FieldVector<CoordScalar,dimWorld> GlobalPosition;
    typedef Dune::FieldMatrix<CoordScalar,dimWorld,dimWorld> FieldMatrix;

    typedef typename GridView::template Codim<0>::Entity Element;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(FVElementGeometry)) FVElementGeometry;


public:

    LensSpatialParameters1p2c(const GridView& gridView)
        : ParentType(gridView),
          lensK_(0),
          outerK_(0)
    {
        Dumux::InterfaceSoilProperties interfaceSoilProps("interface1p2c.xml");

        lensPorosity_ = interfaceSoilProps.ISP_FinePorosity;
        outerPorosity_ = interfaceSoilProps.ISP_CoarsePorosity;

        longitudinalDispersivity_ = 1.0e-5;
        transverseDispersivity_ = 1.0e-6;
    //Remark: The example is very bad. The numerical diffusion is very high, so that the dispersion/diffusion coefficients nearly do not have any influence.


        for(int i = 0; i < dim; i++){
            lensK_[i][i] = interfaceSoilProps.ISP_FinePermeability;
            outerK_[i][i] = interfaceSoilProps.ISP_CoarsePermeability;
        }
    }

    /*!
     * \brief Apply the intrinsic permeability tensor to a pressure
     *        potential gradient.
     *
     * \param element The current finite element
     * \param fvElemGeom The current finite volume geometry of the element
     * \param scvIdx The index sub-control volume face where the
     *                      intrinsic velocity ought to be calculated.
     */
    const FieldMatrix& intrinsicPermeability(const Element &element,
                                 const FVElementGeometry &fvElemGeom,
                                 int scvIdx) const
    {
        const GlobalPosition &globalPos = fvElemGeom.subContVol[scvIdx].global;
        if (isInLens_(globalPos))
            return lensK_;
        return outerK_;
    }

    Scalar porosity(const Element &element,
                    const FVElementGeometry &fvElemGeom,
                    int scvIdx) const
    {
        const GlobalPosition &globalPos = fvElemGeom.subContVol[scvIdx].global;
        if (isInLens_(globalPos))
            return lensPorosity_;
        return outerPorosity_;
    }

    //! Set the bounding box of the fine-sand lens
    void setLensCoords(const GlobalPosition& lensLowerLeft,
                       const GlobalPosition& lensUpperRight)
    {
        lensLowerLeft_ = lensLowerLeft;
        lensUpperRight_ = lensUpperRight;
    }

    /*!
     * \brief Define the tortuosity \f$[?]\f$.
     *
     * \param element The finite element
     * \param fvElemGeom The finite volume geometry
     * \param scvIdx The local index of the sub-control volume where
     */
    Scalar tortuosity(const Element &element,
                    const FVElementGeometry &fvElemGeom,
                    int scvIdx) const
    {
        return 0;
    }

    /*!
     * \brief Define the dispersivity \f$[m]\f$.
     *
     * \param element The finite element
     * \param fvElemGeom The finite volume geometry
     * \param scvIdx The local index of the sub-control volume where
     */
    Dune::FieldVector<Scalar,2> dispersivity(const Element &element,
            const FVElementGeometry &fvElemGeom,
            int scvIdx) const
    {
        Dune::FieldVector<Scalar,2> values (0);
        values[0] = longitudinalDispersivity_; // alpha_l
        values[1] = transverseDispersivity_; // alpha_t

        return values;
    }

    bool useTwoPointGradient(const Element &elem,
                             int vertexI,
                             int vertexJ) const
    {
        return false;
    }


private:
    bool isInLens_(const GlobalPosition &pos) const
    {
        for (int i = 0; i < dim; ++i) {
            if (pos[i] < lensLowerLeft_[i] || pos[i] > lensUpperRight_[i])
                return false;
        }
        return true;
    }

    GlobalPosition lensLowerLeft_;
    GlobalPosition lensUpperRight_;

    FieldMatrix lensK_;
    FieldMatrix outerK_;
    Scalar lensPorosity_;
    Scalar outerPorosity_;
    Scalar longitudinalDispersivity_;
    Scalar transverseDispersivity_;
};

} // end namespace
#endif
