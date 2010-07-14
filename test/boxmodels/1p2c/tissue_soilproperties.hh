// $Id: tissue_soilproperties.hh 3783 2010-06-24 11:33:53Z bernd $
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
#ifndef DUMUX_TISSUE_SOILPROPERTIES_HH
#define DUMUX_TISSUE_SOILPROPERTIES_HH

#include <dumux/old_material/property_baseclasses.hh>

namespace Dumux
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

    const Dune::FieldMatrix<CoordScalar,dim,dim> &K (const Dune::FieldVector<CoordScalar,dim>& globalPos, const Element& element, const Dune::FieldVector<CoordScalar,dim>& localPos) const
    {
        if (isTumor_(globalPos))
            return permTumor_;
        else
            return permTissue_;
    }

    double porosity(const Dune::FieldVector<CoordScalar,dim> &globalPos, const Element &element, const Dune::FieldVector<CoordScalar,dim> &localPos) const
    {
        if (isTumor_(globalPos))
            return 0.31;        // tumor tissue
        else
            return 0.13;        // healthy tissue
    }

    double tortuosity (const Dune::FieldVector<CoordScalar,dim>& globalPos, const Element& element, const Dune::FieldVector<CoordScalar,dim>& localPos) const
    {
        if (isTumor_(globalPos))
            return 0.706;       // tumor tissue
        else
            return 0.280;       // healthy tissue
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

    bool useTwoPointGradient(const Element &elem,
                             int vertexI,
                             int vertexJ) const
    {
        bool inTumor = isTumor_(elem.geometry().corner(0));
        int n = elem.template count<dim>();
        for (int i = 1; i < n; ++i) {
            if ((inTumor && !isTumor_(elem.geometry().corner(i))) ||
                (!inTumor && isTumor_(elem.geometry().corner(i))))
                return true;
        }
        return false;
    }

    TissueSoil()
        :Matrix2p<Grid,Scalar>()
    {
        permTumor_ = 0;
        permTissue_ = 0;
        for (int i = 0; i < dim; i++)
            permTumor_[i][i] = 2.142e-11;         //tumor tissue
        for (int i = 0; i < dim; i++)
            permTissue_[i][i] = 4.424e-12;        //healthy tissue
    }

private:
    bool isTumor_(const Dune::FieldVector<CoordScalar,dim>& globalPos) const
    {
        if(globalPos[0] > 0.99 && globalPos[0] < 2.01 &&
           globalPos[1] > 0.99 && globalPos[1] < 3.01)
            return true;
        return false;
    }

    Dune::FieldMatrix<CoordScalar,dim,dim> permTumor_;
    Dune::FieldMatrix<CoordScalar,dim,dim> permTissue_;
};

} // end namespace

#endif
