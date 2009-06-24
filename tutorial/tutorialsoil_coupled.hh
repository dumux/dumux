// $Id$
/*****************************************************************************
 *   Copyright (C) 2008-2009 by Melanie Darcis                               *
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
#ifndef TUTORIALSOIL_COUPLED_HH
#define TUTORIALSOIL_COUPLED_HH

#include <dumux/material/matrixproperties.hh>

namespace Dune
{

template<class Grid, class Scalar>
class TutorialSoil: public HomogeneousSoil<Grid, Scalar> /*@\label{tutorial-coupled:tutorialsoil}@*/
{
public:
    typedef    typename Grid::Traits::template Codim<0>::Entity Element;

    enum{dim=Grid::dimension};

    // method returning the intrinsic permeability tensor K depending
    // on the position within the domain
    const FieldMatrix<Scalar,dim,dim> &K (const FieldVector<Scalar,dim>& globalPos, const Element& element, /*@\label{tutorial-coupled:permeability}@*/
                                  const FieldVector<Scalar,dim>& localPos) const
    {
        return K_;
    }

    // method returning the porosity of the porous matrix depending on
    // the position within the domain
    double porosity(const FieldVector<Scalar,dim>& globalPos, const Element& element, /*@\label{tutorial-coupled:porosity}@*/
                    const FieldVector<Scalar,dim>& localPos) const
    {
        return 0.2;
    }

    // method returning the residual saturation of the wetting fluid
    // depending on the position within the domain and on the
    // temperature
    double Sr_w(const FieldVector<Scalar,dim>& globalPos, const Element& element, /*@\label{tutorial-coupled:srw}@*/
                const FieldVector<Scalar,dim>& localPos, const double T = 283.15) const
    {
        return 0;
    }

    // method returning the residual saturation of the non-wetting
    // fluid depending on the position within the domain and on the
    // temperature
    double Sr_n(const FieldVector<Scalar,dim>& globalPos, const Element& element, /*@\label{tutorial-coupled:srn}@*/
                const FieldVector<Scalar,dim>& localPos, const double T = 283.15) const
    {
        return 0;
    }

    // method returning the parameters of the capillary pressure and
    // the relative permeability functionms depending on the position
    // within the domain and on the temperature
    std::vector<double> paramRelPerm(const FieldVector<Scalar,dim>& globalPos, const Element& element, /*@\label{tutorial-coupled:parameters}@*/
                                     const FieldVector<Scalar,dim>& localPos, const double T = 283.15) const
    {
        std::vector<double> param(2);

        //linear law parameters
        param[0] = 0; // minimal capillary pressure
        param[1] = 0; // maximal capillary pressure

        //Brooks-Corey parameters
        //        param[0] = 2; // lambda
        //        param[1] = 0.; // entry-pressure

        return param;
    }

    // method returning the kind of relation used for the calculation
    // of the capillary pressure and the relative permeabilities
    // depending on the position within the domain
    typename Matrix2p<Grid,Scalar>::modelFlag relPermFlag(const FieldVector<Scalar,dim>& globalPos, const Element& element, /*@\label{tutorial-coupled:flags}@*/
                                                   const FieldVector<Scalar,dim>& localPos) const
    {
        return Matrix2p<Grid,Scalar>::linear; //flag types defined in
    }                                   //dumux/material/property_baseclasses.hh

    TutorialSoil()
        :HomogeneousSoil<Grid,Scalar>(), K_(0)
    {
        for(int i = 0; i < dim; i++)
            K_[i][i] = 1e-7;
    }

private:
    FieldMatrix<Scalar,dim,dim> K_;
};
} // end namespace
#endif
