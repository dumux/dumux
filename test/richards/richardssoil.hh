/*****************************************************************************
 *   Copyright (C) 2009 by Onur Dogan                                        *
 *   Copyright (C) 2009 by Andreas Lauser                                    *
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
#ifndef DUNE_NEW_RICHARDS_SOIL_HH
#define DUNE_NEW_RICHARDS_SOIL_HH

#include <dumux/material/property_baseclasses.hh>
#include <dumux/material/relperm_pc_law.hh>
#include <dumux/material/matrixproperties.hh>
#include <dumux/material/twophaserelations.hh>

namespace Dune
{

/** \todo Please doc me! */
template<class Grid, class ScalarT>
class RichardsSoil : public Matrix2p<Grid, ScalarT>
{
public:
    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef ScalarT Scalar;
    typedef typename Grid::ctype CoordScalar;
    enum {dim=Grid::dimension, dimWorld=Grid::dimensionworld};

    typedef Dune::FieldVector<CoordScalar,dim>      LocalPosition;
    typedef Dune::FieldVector<CoordScalar,dimWorld> GlobalPosition;

    RichardsSoil():Matrix2p<Grid,Scalar>()
    {
        Kout_ = 0.;
        for(int i = 0; i < dim; i++)
            Kout_[i][i] = 5e-10;
    }

    ~RichardsSoil()
    {}

    const FieldMatrix<CoordScalar,dim,dim> &K (const GlobalPosition &x, const Element& e, const LocalPosition &xi) const
    {
        return Kout_;
    }

    double porosity(const GlobalPosition &x, const Element& e, const LocalPosition &xi) const
    {
        return 0.4;
    }

    double Sr_w(const GlobalPosition &x, const Element& e, const LocalPosition &xi, const double T) const
    {
        return 0.05;
    }

    double Sr_n(const GlobalPosition &x, const Element& e, const LocalPosition &xi, const double T) const
    {
        return 0.0;
    }

    /* ATTENTION: define heat capacity per cubic meter! Be sure, that it corresponds to porosity!
     * Best thing will be to define heatCap = (specific heatCapacity of material) * density * porosity*/
    double heatCap(const GlobalPosition &x, const Element& e, const LocalPosition &xi) const
    {
        return     790 /* spec. heat cap. of granite */
            * 2700 /* density of granite */
            * porosity(x, e, xi);
    }

    double heatCond(const GlobalPosition &x, const Element& e, const LocalPosition &xi, const double sat) const
    {
        static const double lWater = 0.6;
        static const double lGranite = 2.8;
        double poro = porosity(x, e, xi);
        double lsat = pow(lGranite, (1-poro)) * pow(lWater, poro);
        double ldry = pow(lGranite, (1-poro));
        return ldry + sqrt(sat) * (ldry - lsat);
    }

    std::vector<double> paramRelPerm(const GlobalPosition &x, const Element& e, const LocalPosition &xi, const double T) const
    {
        // example for Brooks-Corey parameters
        std::vector<double> param(2);
        param[0] = 0; // pCMin
        param[1] = 1e5; // pCMax

        return param;
    }

    typename Matrix2p<Grid,Scalar>::modelFlag relPermFlag(const GlobalPosition &x, const Element& e, const LocalPosition &xi) const
    {
        return Matrix2p<Grid,Scalar>::linear;
    }

private:
    FieldMatrix<Scalar,dim,dim> Kout_;
};

} // end namespace

#endif
