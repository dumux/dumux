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

#include <dumux/material/spatialparameters/fvspatialparameters1p.hh>

namespace Dumux
{

struct Lens
{
    Dune::FieldMatrix<double,2,2> permeability;
    //double permeability;
    Dune::FieldVector<double,2> lowerLeft;
    Dune::FieldVector<double,2> upperRight;
};


/*!
 * \ingroup IMPETtests
 * \brief spatial parameters for the test problem for diffusion models.
 */
template<class TypeTag>
class GroundwaterSpatialParams: public FVSpatialParametersOneP<TypeTag>
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
    typedef typename GET_PROP(TypeTag, PTAG(ParameterTree)) Params;

public:

    const FieldMatrix& intrinsicPermeabilityAtPos (const GlobalPosition& globalPos) const
    {
        if (lenses_.size())
        {
            for (int lensnumber = 0; lensnumber != lenses_.size(); lensnumber++)
            {
                if ((lenses_[lensnumber].lowerLeft[0] <= globalPos[0])
                  && (lenses_[lensnumber].lowerLeft[1] <= globalPos[1])
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
    : FVSpatialParametersOneP<TypeTag>(gridView), permeability_(0)
    {
        delta_=1e-3;
        porosity_ = 0.2;

        Scalar permFactor = 0.001/(1000*9.81);

        permeability_[0][0] = Params::tree().template get<double>("SpatialParameters.permeability")*permFactor;
        permeability_[1][1] = permeability_[0][0];
        permeability_[0][1] = 0;
        permeability_[1][0] = 0;

        //Lenses:
        std::vector<double> lenses = Params::tree().template get<std::vector<double>>("SpatialParameters.lenses");
        int NumberOfLenses = std::trunc(lenses.size()/5);

        for (int lensCount=0; lensCount<NumberOfLenses ; lensCount++)
        {
            Lens tempLens;
            tempLens.lowerLeft[0]=lenses[lensCount*5];
            tempLens.upperRight[0]=lenses[lensCount*5+1];
            tempLens.lowerLeft[1]=lenses[lensCount*5+2];
            tempLens.upperRight[1]=lenses[lensCount*5+3];
            tempLens.permeability[0][0]=lenses[lensCount*5+4]*permFactor;
            tempLens.permeability[1][1] = tempLens.permeability[0][0];
            tempLens.permeability[0][1] = 0;
            tempLens.permeability[1][0] = 0;
            lenses_.push_back(tempLens);
        }
    }

private:
    mutable FieldMatrix permeability_;
    Scalar porosity_;
    double delta_;
    std::vector <Lens> lenses_;
};

} // end namespace
#endif
