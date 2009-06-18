// $Id$
/*****************************************************************************
 *   Copyright (C) <YEARS> by <ADD_AUTHOR_HERE>                              *
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
#ifndef TUTORIAL_SOILPROPERTIES
#define TUTORIAL_SOILPROPERTIES

#include <dumux/material/matrixproperties.hh>

namespace Dune
{

/** \todo Please doc me! */

template<class Grid, class Scalar>
class TutorialSoil: public HomogeneousSoil<Grid, Scalar> /*@\label{tutorial-decoupled:tutorialsoil}@*/
{
public:
    typedef    typename Grid::Traits::template Codim<0>::Entity Element;

    enum{dim=Grid::dimension, dimWorld=Grid::dimensionworld};
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldVector<Scalar, dim> LocalPosition;
    typedef Dune::FieldMatrix<Scalar,dim,dim> FieldMatrix;

    // function returning the intrinsic permeability tensor K
    // depending on the position within the domain
    const FieldMatrix &K(const GlobalPosition& globalPos, const Element& element, /*@\label{tutorial-decoupled:permeability}@*/
                         const LocalPosition& localPos)
    {
        return K_;
    }

    // function returning the porosity of the porous matrix
    // depending on the position within the domain
    double porosity(const GlobalPosition& globalPos, const Element& element, /*@\label{tutorial-decoupled:porosity}@*/
                    const LocalPosition& localPos) const
    {
        return 0.2;
    }

    // function returning the residual saturation of the wetting fluid
    // depending on the position within the domain and on the temperature
    double Sr_w(const GlobalPosition& globalPos, const Element& element, /*@\label{tutorial-decoupled:srw}@*/
                const LocalPosition& localPos, const double T = 283.15) const
    {
        return 0;
    }

    // function returning the residual saturation of the non-wetting fluid
    // depending on the position within the domain and on the temperature
    double Sr_n(const GlobalPosition& globalPos, const Element& element, /*@\label{tutorial-decoupled:srn}@*/
                const LocalPosition& localPos, const double T = 283.15) const
    {
        return 0;
    }

    // function returning the parameters of the capillary pressure
    // and the relative permeability functions
    // depending on the position within the domain and on the temperature
    std::vector<double> paramRelPerm(const GlobalPosition& globalPos, const Element& element, /*@\label{tutorial-decoupled:parameters}@*/
                                     const LocalPosition& localPos, const double T = 283.15) const
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

    // function returning the kind of relation used for the calculation of the capillary
    // pressure and the relative permeabilities depending on the position within the domain
    typename Matrix2p<Grid,Scalar>::modelFlag relPermFlag(const GlobalPosition& globalPos, const Element& e, /*@\label{tutorial-decoupled:flags}@*/
                                                          const LocalPosition& localPos) const
    {
        return Matrix2p<Grid,Scalar>::linear; //flag types defined in
    }                                   //dumux/material/property_baseclasses.hh

    TutorialSoil()
        :HomogeneousSoil<Grid,Scalar>(),K_(0)
    {
        for(int i = 0; i < dim; i++)
            K_[i][i] = 1e-7;
    }

private:
    FieldMatrix K_;

};
} // end namespace
#endif
