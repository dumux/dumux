// $Id$
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
#ifndef DUNE_ONEP_SOIL_HH
#define DUNE_ONEP_SOIL_HH

namespace Dune
{

/**
 * Defines soil properties for single phase flow problem
 *
 */


//////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////--SOIL--//////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////////////////////////////

/** \todo Please doc me! */

template<class Grid, class ScalarT>
class OnePSoil
{
public:
    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef ScalarT Scalar;
    typedef typename Grid::ctype CoordScalar;
    enum {dim=Grid::dimension, dimWorld=Grid::dimensionworld};

    typedef Dune::FieldVector<CoordScalar,dim>      LocalPosition;
    typedef Dune::FieldVector<CoordScalar,dimWorld> GlobalPosition;

    OnePSoil()
    {
        Kout_ = 0.;
        for(int i = 0; i < dim; i++)
            Kout_[i][i] = 1.0e-10;
    }

    ~OnePSoil()
    {}

    const FieldMatrix<CoordScalar,dim,dim> &K (const GlobalPosition &x, const Element& e, const LocalPosition &xi) const
    {
        return Kout_;
    }

    double porosity(const GlobalPosition &x, const Element& e, const LocalPosition &xi) const
    {
        return 0.4;
    }

private:
    FieldMatrix<Scalar,dim,dim> Kout_;
};

} // end namespace

#endif
