/*****************************************************************************
 *   Copyright (C) 2008-2009 by Markus Wolff                                 *
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
/*!
 * \file
 *
 * \brief spatial parameters for the GRUWA-equivalent exercise
 */
#ifndef GROUNDWATER_SPATIALPARAMETERS_HH
#define GROUNDWATER_SPATIALPARAMETERS_HH

#include <dumux/material/spatialparameters/fvspatialparameters.hh>

namespace Dumux
{

/*!
 * \ingroup IMPETtests
 * \brief spatial parameters for the test problem for diffusion models.
 */
template<class TypeTag>
class GroundwaterSpatialParams: public FVSpatialParameters<TypeTag>
{
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Grid)) Grid;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename Grid::ctype CoordScalar;

    enum
        {dim=Grid::dimension, dimWorld=Grid::dimensionworld, numEq=1};
    typedef typename Grid::Traits::template Codim<0>::Entity Element;

    typedef Dune::FieldVector<CoordScalar, dimWorld> GlobalPosition;
    typedef Dune::FieldVector<CoordScalar, dim> LocalPosition;
    typedef Dune::FieldMatrix<Scalar,dim,dim> FieldMatrix;


public:

    const FieldMatrix& intrinsicPermeabilityAtPos (const GlobalPosition& globalPos) const
    {
        if (lenses_.size())
        {
            for (int lensnumber = 0; lensnumber != lenses_.size(); lensnumber++)
            {
                if ((lenses_[lensnumber].lowerLeft[0] < globalPos[0])
                  && (lenses_[lensnumber].lowerLeft[1] < globalPos[1])
                  && (lenses_[lensnumber].upperRight[0] > globalPos[0])
                  && (lenses_[lensnumber].upperRight[1] > globalPos[1]))
                {
                    return lenses_[lensnumber].permeability;
                }


            }
        }
        // Schleife über alle Linsen
        return permeability_;
    }

    double porosity(const Element& element) const
    {
        return porosity_;
    }

    GroundwaterSpatialParams(const GridView& gridView)
    : FVSpatialParameters<TypeTag>(gridView), permeability_(0)
    {
        Dumux::InterfaceSoilProperties interfaceSoilProps("interface_groundwater.xml");
        porosity_ = interfaceSoilProps.porosity;
        permeability_[0][0] = interfaceSoilProps.permeability;
        permeability_[0][1] = 0;
        permeability_[1][0] = 0;
        permeability_[1][1] = interfaceSoilProps.permeability;

        lenses_ = interfaceSoilProps.lenses;
    }

private:
    mutable FieldMatrix permeability_;
    Scalar porosity_;
    std::vector <Lens> lenses_;
};

} // end namespace
#endif
