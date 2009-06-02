// $Id$
/*****************************************************************************
 *   Copyright (C) 2009 by Karin Erbertseder                                 *
 *   Copyright (C) 2009 by Andreas Lauser                                    *
 *   Copyright (C) 2008 by Bernd Flemisch                                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
#ifndef DUNE_TISSUE_SOILPROPERTIES_HH
#define DUNE_TISSUE_SOILPROPERTIES_HH

#include <dumux/material/property_baseclasses.hh>

namespace Dune
{

/** 
 * \brief Definition of the properties of the human tissue 
 *
 * \todo do not derive from Matrix2p, since this is the base class for
 *       two-phase models
 */
template<class Grid, class Scalar>
class TissueSoil: public Matrix2p<Grid, Scalar>
{
public:
    typedef    typename Grid::Traits::template Codim<0>::Entity Element;
    typedef typename Grid::ctype CoordScalar;
    enum
        {    dim=Grid::dimension, numEq=2};

    const FieldMatrix<CoordScalar,dim,dim> &K (const Dune::FieldVector<CoordScalar,dim>& globalPos, const Element& element, const Dune::FieldVector<CoordScalar,dim>& localPos) const
    {
        if (10<globalPos[0] && globalPos[0]<12 && 10<globalPos[1] && globalPos[1]<12)
            return permloc_;                        //tumor tissue
        else
            return permlocWell_;                    //healthy tissue
    }

    double porosity(const Dune::FieldVector<CoordScalar,dim> &globalPos, const Element &element, const Dune::FieldVector<CoordScalar,dim> &localPos) const
    {
        if (10<globalPos[0] && globalPos[0]<12 && 10<globalPos[1] && globalPos[1]<12)
            return 0.31;                        //tumor tissue
        else
            return 0.13;                    //healthy tissue
    }

    double tortuosity (const Dune::FieldVector<CoordScalar,dim>& globalPos, const Element& element, const Dune::FieldVector<CoordScalar,dim>& localPos) const
    {
        if (10<globalPos[0] && globalPos[0]<12 && 10<globalPos[1] && globalPos[1]<12)
            return 0.706;                        //tumor tissue
        else
            return 0.280;                    //healthy tissue
    }

    virtual double Sr_w(const Dune::FieldVector<CoordScalar,dim>& globalPos, const Element& element, const Dune::FieldVector<CoordScalar,dim>& localPos, const double T = 283.15) const
    {
        return 0;
    }

    virtual double Sr_n(const Dune::FieldVector<CoordScalar,dim>& globalPos, const Element& element, const Dune::FieldVector<CoordScalar,dim>& localPos, const double T = 283.15) const
    {
        return 0;
    }

    virtual std::vector<double> paramRelPerm(const Dune::FieldVector<CoordScalar,dim>& globalPos, const Element& element, const Dune::FieldVector<CoordScalar,dim>& localPos, const double T = 283.15) const
    {
        std::vector<double> param(0);
        return param;
    }

    TissueSoil()
        :Matrix2p<Grid,Scalar>()
    {
        permloc_ = 0;
        permlocWell_ = 0;
        for (int i = 0; i < dim; i++)
            permloc_[i][i] = 2.142e-11;             //tumor tissue
        for (int i = 0; i < dim; i++)
            permlocWell_[i][i] = 4.424e-12;        //healthy tissue
    }

private:
    Dune::FieldMatrix<CoordScalar,dim,dim> permloc_, permlocWell_;
};

} // end namespace

#endif
