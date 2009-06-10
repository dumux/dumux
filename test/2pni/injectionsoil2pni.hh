//$Id$
/*****************************************************************************
 *   Copyright (C) 2007-2008 by Melanie Darcis                               *
 *   Copyright (C) 2007-2008 by Bernd Flemisch                               *
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
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
#ifndef DUNE_INJECTIONSOIL2PNI_HH
#define DUNE_INJECTIONSOIL2PNI_HH

/*!
 * \file
 *
 * \brief This file contains the soil parameters which are required for
 *        the non-isothermal two-phase injection problem
 */

namespace Dune {
/*!
 * \ingroup TwoPNIBoxProblems
 * \brief This file contains the soil parameters which are required for
 *        the non-isothermal two-phase injection problem
 *
 * This means pressure and concentration gradients, phase densities at
 * the intergration point, etc.
 */

template<class Grid, class Scalar>
class InjectionSoil: public Matrix2p<Grid, Scalar> {
public:
	enum {
		// Grid and world dimension
		dim = Grid::dimension,
		dimWorld = Grid::dimensionworld,
	};
typedef	typename Grid::Traits::template Codim<0>::Entity Element;
	typedef Dune::FieldVector<Scalar, dim> LocalPosition;
	typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

	/*!
	 * \brief Matrix permeability
	 */
	virtual const FieldMatrix<Scalar,dim,dim> &K (const GlobalPosition &globalPos, const Element& element, const LocalPosition &localPos) const
	{
		if ((globalPos[0]>= innerLowerLeft_[0] && globalPos[0] <= innerUpperRight_[0])
				&& (globalPos[1]>= innerLowerLeft_[1] && globalPos[1] <= innerUpperRight_[1]))
		return Kin_;
		else
		return Kout_;
	}

	/*!
	 * \brief Matrix porosity
	 */
	virtual double porosity(const GlobalPosition &globalPos, const Element& element, const LocalPosition &localPos) const
	{
		return 0.15;
	}

	/*!
	 * \brief Residual saturation of wetting phase
	 */
	virtual double Sr_w(const GlobalPosition &globalPos, const Element& element, const LocalPosition &localPos, const double T) const
	{
		return 0.2;
	}

	/*!
	 * \brief Residual saturation of non-wetting phase
	 */
	virtual double Sr_n(const GlobalPosition &globalPos, const Element& element, const LocalPosition &localPos, const double T) const
	{
		return 0.05;
	}

	/*!
	 * \brief Modified heat capacity of porous matrix
	 * 		  volume fraction of rock and rock density are already included
	 */
	virtual double heatCap(const GlobalPosition &x, const Element& e, const LocalPosition &xi) const
	{
		return (800 // spec. heat cap. of sediment
				* 2650 // density of sediment
				* (1-porosity(x, e, xi)));
	}

	/*!
	 * \brief Heat conductivity of porous medium
	 * 		  dependent of wetting fluid saturation
	 */
	virtual double heatCond(const GlobalPosition &x, const Element& e, const LocalPosition &xi, const double sat) const
	{
		static const double ldry = 0.32;
		static const double lsat = 2.7;
		double l_wurz, Sw;
		Sw = sat;
		if(Sw<0.0) Sw = 0.0; // regularisation
		if(Sw>1.) Sw = 1.; // regularisation

		l_wurz = ldry + sqrt(Sw)*(lsat-ldry);

		return(l_wurz);
	}

	/*!
	 * \brief Flag for capillary pressure and relative permeability
	 * 		  function (e.g. Van Genuchten, Brooks Corey ..)
	 */
	virtual typename Matrix2p<Grid,Scalar>::modelFlag relPermFlag(const GlobalPosition &globalPos, const Element& element, const LocalPosition &localPos) const
	{
		return Matrix2p<Grid,Scalar>::brooks_corey;
	}

	/*!
	 * \brief Parameters for capillary pressure and relative permeability
	 * 		  function
	 */
	virtual std::vector<double> paramRelPerm(const GlobalPosition &globalPos, const Element& element, const LocalPosition &localPos, const double T) const
	{
		// example for Brooks-Corey parameters
		std::vector<double> param(2);

		param[0] = 2.; // lambda
		param[1] = 10000.; // entry-pressure

		return param;
	}

	InjectionSoil()
        : Matrix2p<Grid,Scalar>()
	{
		Kin_ = Kout_ = 0;
		for(int i = 0; i < dim; i++)
		{
			Kin_[i][i] = 1e-18;
			Kout_[i][i] = 2e-14;
		}
	}

	~InjectionSoil()
	{}

    //! Set the bounding box of the fine-sand lens
    void setLensCoords(const FieldVector<Scalar,dim>& innerLowerLeft,
                       const FieldVector<Scalar,dim>& innerUpperRight)
    {
        innerLowerLeft_ = innerLowerLeft;
        innerUpperRight_ = innerUpperRight;
    }


private:
	FieldMatrix<Scalar,dim,dim> Kin_;
	FieldMatrix<Scalar,dim,dim> Kout_;
	FieldVector<Scalar,dim> innerLowerLeft_;
	FieldVector<Scalar,dim> innerUpperRight_;
};

} // end namespace
#endif

