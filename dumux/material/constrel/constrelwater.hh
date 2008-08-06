#ifndef DUNE_CONSTRELWATER_HH
#define DUNE_CONSTRELWATER_HH

#include "constrelco2.hh"

/************************************************************************/
/*									*/
/*	Calculation of the mass density of pure water                   */
/*									*/
/*	after IAPWS 1997 (only region 1)				*/
/*									*/
/*	validity: temperature: 273.15 K <= T <= 1073.15 K		*/
/*		  pressure : p <= 100 MPa (1000 bar, 1.0E8Pa)		*/
/*									*/
/*	(IAPWS: The international Association for the properties	*/
/*		of water and steam)					*/ 
/*									*/
/*	(http://www.iapws.org)						*/
/*									*/
/************************************************************************/

double I[35] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 1.0, 1.0, 1.0, 
		1.0, 1.0, 2.0, 2.0, 2.0, 2.0, 2.0, 3.0, 3.0, 3.0, 4.0, 4.0, 
		4.0, 5.0, 8.0, 8.0, 21.0, 23.0, 29.0, 30.0, 31.0, 32.0};
double J[35] = {0.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0, 4.0, 5.0, -9.0, -7.0, -1.0, 
		0.0, 1.0, 3.0, -3.0, 0.0, 1.0, 3.0, 17.0, -4.0, 0.0, 6.0, 
		-5.0, -2.0, 10.0, -8.0, -11.0, -6.0, -29.0, -31.0, -38.0, -39.0, -40.0, -41.0};
double n[35] = {0.0, 0.14632971213167, -0.84548187169114, -0.37563603672040E1, 0.33855169168385E1, 
		-0.95791963387872, 0.15772038513228, -0.16616417199507E-1, 0.81214629983568E-3, 
		0.28319080123804E-3, -0.60706301565874E-3, -0.18990068218419E-1, -0.32529748770505E-1, 
		-0.21841717175414E-1, -0.52838357969930E-4, -0.47184321073267E-3, -0.30001780793026E-3, 
		0.47661393906987E-4, -0.44141845330846E-5, -0.72694996295794E-15, -0.31679644845054E-4, 
		-0.28270797985312E-5, -0.85205128120103E-9, -0.22425281908000E-5, -0.65171222895601E-6, 
		-0.14341729937924E-12, -0.40516996860117E-6, -0.12734301741641E-8, -0.17424871230634E-9, 
		-0.68762131295531E-18, 0.14478307828521E-19, 0.26335781662795E-22, -0.11947622640071E-22, 
		0.18228094581404E-23, -0.93537087292458E-25};



class ConstrelWater
{
public: 
	
	double mass_density_water_IAPWS (double Temp, double pw) const
	{
		double gam_pi;
		double pi, p_star;
		double R;
		double v;	/* specific volume */
		double rho;	/* mass density */

	        /* Regularisierung */
	        if (Temp < 273.15) Temp = 273.15;
	        if (Temp > 1073.15) Temp = 1073.15;
	        if (pw < 1.0) pw = 1.0;
	        if (pw > 1.0E8) pw = 1.0E8;

		R = 461.526;	   /* [J/kg K] specific gas constant for ordinary water */
		p_star = 16.53E6;	/* [Pa] */ 

		pi = pw / p_star;

		gam_pi = gamma_pi(pi, Temp);

		v = pi * gam_pi * R * Temp / pw;	
		
		rho = 1 / v;	/* kg / m^3 */
		
		return(rho);  /* unit: [kg/m^3] */
	}
	double gamma_pi(double pi, double Temp) const
	{
	        int i;
	        double T_star;
	        double tau;
	        double gam_pi;

	        T_star = 1386.0;        /* [K] */

	        tau = T_star / Temp;    /* reduced temperature */

	        gam_pi = 0.0;

	        for (i = 1; i <= 34; i++)
	        {
	                gam_pi = gam_pi - n[i] * I[i] * pow((7.1 - pi), (I[i] - 1)) * pow((tau - 1.222), J[i]);
	        }

	        return(gam_pi);
	}


	/************************************************************************/
	/*									*/
	/*	Calculation of the enthalpy of pure water                       */
	/*									*/
	/*	after IAPWS 1997 (only region 1)				*/
	/*									*/
	/*	validity: temperature: 273.15 K <= T <= 1073.15 K		*/
	/*		  pressure : p <= 100 MPa (1000 bar, 1.0E8Pa)		*/
	/*									*/
	/*	(IAPWS: The international Association for the properties	*/
	/*		of water and steam)					*/ 
	/*									*/
	/*	(http://www.iapws.org)						*/
	/*									*/
	/************************************************************************/

	double enthalpy_water (double Temp, double pw) const
	{
		
		double gam_tau;
		double tau, T_star;
		double R;
		double h; /* specific enthalpy */
		
	        /* Regularisierung */
	        if (Temp < 273.15) Temp = 273.15;
	        if (Temp > 623.15) Temp = 623.15;
	        if (pw < 1.0) pw = 1.0;
	        if (pw > 1.0E8) pw = 1.0E8;

		R = 461.526;	   /* [J/kg K] specific gas constant for ordinary water */
		T_star = 1386;	/* [K] */ 

		tau = T_star / Temp;	/* reduced temperature */
		
		gam_tau = gamma_tau(tau, pw);

		h = tau * gam_tau * R * Temp;	/* J/kg */
			
		return(h);
	}


	double gamma_tau (double tau, double pw) const
	{
	        int i;
	        double p_star;
	        double pi;
	        double gam_tau;


		/* Regularisierung */
	        if (pw < 1.0) pw = 1.0;
	        if (pw > 1.0E8) pw = 1.0E8;


	        p_star = 1.653E7;        /* [Pa] */

	        pi = pw / p_star;    /* reduced pressure */

	        gam_tau = 0.0;

	        for (i = 1; i <= 34; i++)
	        {
	                gam_tau = gam_tau + n[i] * pow((7.1 - pi), I[i]) *J[i]*
			pow((tau - 1.222), (J[i]-1));
	        }

	        return(gam_tau);
	}
	
	double gamma_tau_1(double Temp, double pw)
	{
		static const double I[35]={0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 1., 1., 1., 1., 1., 2., 2.,
				 2., 2., 2., 3., 3., 3., 4., 4., 4., 5., 8., 8., 21., 23., 29., 30., 31., 32.};
		
	  static const double J[35]={0.,-2.,-1.,0.,1.,2.,3.,4.,5.,-9.,-7.,-1.,0.,1.,3.,-3.,0.,1.,3.,17.,
			-4.,0.,6.,-5.,-2.,10.,-8.,-11.,-6.,-29.,-31.,-38.,-39.,-40.,-41.};

	  static const double n[35]={0.,0.14632971213167, -0.84548187169114, -0.37563603672040E1, 0.33855169168385E1,
			-0.95791963387872, 0.15772038513228, -0.16616417199501E-1, 0.81214629983568E-3,
			0.28319080123804E-3, -0.60706301565874E-3, -0.18990068218419E-1, -0.32529748770505E-1,
			-0.21841717175412E-1, -0.52838357969930E-4, -0.47184321073267E-3, -0.30001780793026E-3,
			0.47661393906987E-4, -0.44141845330846E-5, -0.72694996297594E-15, -0.31679644845054E-4,
			-0.28270797985312E-5, -0.85205128120103E-9, -0.22425281908000E-5, -0.65171222895601E-6,
			-0.14341729937924E-12, -0.40516996860117E-6, -0.12734301741641E-8, -0.17424871230634E-9,
			-0.68762131295531E-18, 0.14478307828521E-19, 0.26335781662795E-22, -0.11947622640071E-22,
			0.18228094581404E-23, -0.93537087292458E-25};

	   double pr = 16.53E6;		// reduced pressure [Pa]
	   double Tr = 1386.; 		// reduced Temperature [K]
	   double tau = Tr/Temp;
	   double pi = pw/pr;

	   double gamma_tau;
	   int i;
		
	  /*Calculation*/
	   gamma_tau = 0.0;


		for (i=1; i<=34; i++){
			gamma_tau =  gamma_tau + n[i] * pow((7.1 - pi),I[i]) * J[i] * pow((tau - 1.222),(J[i] -1.));
		}

	  return (gamma_tau);
	}


	/**************************************************************************/
	/*                                                                        */
	/* Computation of the mass density of pure water after                    */
	/* Batzle & Wang (1992)                                                   */
	/* equation given by Adams & Bachu in Geofluids (2002) 2, 257-271         */
	/*                                                                        */
	/* comment: this function might not be quite as exact as the IAPWS        */
	/*          formulation, but it comes close and is a lot faster           */
	/*                                                                        */
	/**************************************************************************/

	double mass_density_water_Batzle (double Temp, double pw) const
	{
	    double TempC, pMPa, rhow;

	    TempC = Temp - 273.15;
	    pMPa  = pw/1.0E6;

	    rhow = (1 + 1.0E-6*(-80*TempC - 3.3*pow(TempC,2) + 0.00175*pow(TempC,3)
	             + 489*pMPa - 2*TempC*pMPa + 0.016*pow(TempC,2)*pMPa
	             - 1.3E-5*pow(TempC,3)*pMPa - 0.333*pow(pMPa,2)  
	             - 0.002*TempC*pow(pMPa,2))) * 1000.0 ;

	    return(rhow); /* unit: [kg/m^3] */
	}

	/***************************************************************************/
	/*                                                                         */ 
	/* Computation of the viscosity of pure water: IAPWS, 2003                 */ 
	/* (revised release on the IAPS formulation 1985 for the viscosity of      */ 
	/*  ordinary water substance, www.iapws.org)                               */ 
	/*                                                                         */ 
	/* validity: 273.15 K <= Temp <= 423.15 K                                  */ 
	/*           p <= 5000 bar                                                 */ 
	/*                                                                         */ 
	/***************************************************************************/

	double viscosity_water (double Temp, double pw) const
	{
	  double visco;
	  double mu, mu0, mu1, mu_star;
	  double T_bar, T_star;
	  double rho, rho_bar, rho_star;
	  double sum, H1, H2, H3;
	  double H[6][7];
	  int i,j;

	  rho_star = 317.763;   /* rho at critical point [kg/m^3] */
	  T_star = 647.226;     /* Temp at critical point [K] */
	  mu_star = 55.071E-6;   /* viscosity at critical point [Pa s] */

	  rho = mass_density_water_Batzle (Temp, pw);

	  rho_bar = rho/rho_star;
	  T_bar   = Temp/T_star;

	  /* computation of mu0 */

	  H1 = 0.978197;
	  H2 = 0.579829;
	  H3 = -0.202354;

	  sum = 1.0 + H1/T_bar + H2/(T_bar*T_bar) + H3/(T_bar*T_bar*T_bar);

	  mu0 = sqrt(T_bar)/sum;
	 
	  /* computation of mu1 */  

	  for(i=0; i<=5; i++)
	  {
	    for(j=0; j<=6; j++)
	    {
	      H[i][j]=0.0;
	    }
	  }
	  H[0][0] =  0.5132047;
	  H[1][0] =  0.3205656;
	  H[4][0] = -0.7782567;
	  H[5][0] =  0.1885447;
	  H[0][1] =  0.2151778;
	  H[1][1] =  0.7317883;
	  H[2][1] =  1.2410440;
	  H[3][1] =  1.4767830;
	  H[0][2] = -0.2818107;
	  H[1][2] = -1.0707860;
	  H[2][2] = -1.2631840;
	  H[0][3] =  0.1778064;
	  H[1][3] =  0.4605040;
	  H[2][3] =  0.2340379;
	  H[3][3] = -0.4924179;
	  H[0][4] = -0.0417661;
	  H[3][4] =  0.1600435;
	  H[1][5] = -0.01578386;
	  H[3][6] = -0.003629481;

	  sum = 0.0;

	  for (i=0; i<=5; i++)
	  {
	    for (j=0; j<=6; j++)
	    {
	      sum += H[i][j] * pow((1/T_bar-1), i) * pow((rho_bar-1), j);
	    }
	  }

	  mu1 = exp(rho_bar*sum);

	  mu = mu0*mu1;

	  visco = mu*mu_star;

	  return visco;
	}

	/********************************************************************/
	/*                                                                  */
	/* Density of water depending on the amount of dissolved CO2        */
	/* Garcia (2001)                                                    */
	/*                                                                  */
	/********************************************************************/

	double mass_density_waterCO2 (double Temp, double pw, double x_CO2_w) const
	{
		/* x_CO2_w : mole fraction of CO2 in water phase */

		double rho_pure, rho_inv, rho;
		double xww;
		double T, M_CO2, M_H2O, M_T, V_phi;

		M_CO2 = 0.044;
		M_H2O = 0.018;

		T = Temp - 273.15;		/* T : temperature in °C */

		rho_pure = mass_density_water_Batzle (Temp, pw);

		xww = 1.0 - x_CO2_w;

		M_T = M_H2O * xww + M_CO2 * x_CO2_w;

		V_phi = (37.51 - 9.585E-2*T + 8.74E-4*pow(T,2) - 5.044E-7*pow(T,3)) / 1.0E6;

		rho_inv = x_CO2_w * V_phi/M_T + M_H2O * xww / (rho_pure * M_T);

		rho = 1/rho_inv;

		return(rho);
	}

	/*********************************************************************/
	/* CITATION??????????????????????????????????*************************/
	double lambda_w (double T) const
	{
	double l[11], t[11];
	double downT, downl, upT, upl, lambda;
	int i;

	l[1] = 0.564; 
	l[2] = 0.574;
	l[3] = 0.584;
	l[4] = 0.597;
	l[5] = 0.618;	/*W/mK*/
	l[6] = 0.627;
	l[7] = 0.645;
	l[8] = 0.651;
	l[9] = 0.670;
	l[10] = 0.682;

	t[1] = 273.15;
	t[2] = 278.15;
	t[3] = 283.15;
	t[4] = 293.15;
	t[5] = 303.15;
	t[6] = 313.15;
	t[7] = 323.15;
	t[8] = 333.15;
	t[9] = 353.15;
	t[10] = 373.15;


	for (i=1; i<=10; i++)
	{
		if (t[i] <= T && T <= t[i+1])
		{
		upT = t[i+1];
		downT=t[i];
		upl = l[i+1];
		downl=l[i];
		}
	}

	lambda = downl + ((T-downT)/(upT-downT))*(upl-downl);
	return(lambda);
	}
	/*********************************************************************/


	double lambda_pm (double T, double p, double Sw, double phi) const
	{
		double Sg, l_s, l_w, l_g, l_pm, l_dry, l_sat, llow, lhigh; 
		ConstrelCO2 co2;

		Sg = 1 - Sw;
		l_s = 4.5; /*W/mK*/
		l_w = lambda_w(T);
		l_g = co2.lambda_CO2 (T, p);
		l_pm = l_s*(1- phi) + l_w*Sw*phi + l_g* Sg*phi;
		
	/*wet*/
	/*Hashin-Shtrikman bounds for water and solid*/ /*Hartmann et al 2004*/
	/*HS-*/	llow = l_w + ( (1-phi) / (  ( 1/(l_s-l_w) )+(phi/(3*l_w))  ) );
	/*HS+*/	lhigh = l_s + ( phi / (  ( 1/(l_w-l_s) )+(phi/(3*l_s))  ) );
		
		l_sat = 2 / ( 1/llow + 1/lhigh); /*harmonic averaging of the two
						bounds*/
		
	/*dry*/
	/*Hashin-Shtrikman bounds for gas and solid*/
	/*HS-*/	llow = l_g + ( (1-phi) / (  ( 1/(l_s-l_g) )+(phi/(3*l_g))  ) );
	/*HS+*/	lhigh = l_s + ( phi / (  ( 1/(l_g-l_s) )+(phi/(3*l_s))  ) );
		
		l_dry = 2 / ( 1/llow + 1/lhigh); /*harmonic averaging of the two
						bounds*/
						
		l_pm = l_dry + sqrt(Sw)*(l_sat - l_dry); /*Somerton*/

	return(l_pm);
	}




};
#endif
