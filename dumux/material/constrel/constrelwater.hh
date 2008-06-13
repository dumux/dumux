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
	        double I[35], J[35];
	        double n[35];
	        double T_star;
	        double tau;
	        double gam_pi;

	        T_star = 1386.0;        /* [K] */

	        tau = T_star / Temp;    /* reduced temperature */

	/* list of the numerical values of the coefficients and exponents */

	I[1] = 0.0;     J[1] = -2.0;      n[1] =  0.14632971213167;
	I[2] = 0.0;     J[2] = -1.0;      n[2] = -0.84548187169114;
	I[3] = 0.0;     J[3] =  0.0;      n[3] = -0.37563603672040E1;
	I[4] = 0.0;     J[4] =  1.0;      n[4] =  0.33855169168385E1;
	I[5] = 0.0;     J[5] =  2.0;      n[5] = -0.95791963387872;
	I[6] = 0.0;     J[6] =  3.0;      n[6] =  0.15772038513228;
	I[7] = 0.0;     J[7] =  4.0;      n[7] = -0.16616417199507E-1;
	I[8] = 0.0;     J[8] =  5.0;      n[8] =  0.81214629983568E-3;
	I[9] = 1.0;     J[9] = -9.0;      n[9] =  0.28319080123804E-3;
	I[10] = 1.0;    J[10] = -7.0;     n[10] = -0.60706301565874E-3;
	I[11] = 1.0;    J[11] = -1.0;     n[11] = -0.18990068218419E-1;
	I[12] = 1.0;    J[12] =  0.0;     n[12] = -0.32529748770505E-1;
	I[13] = 1.0;    J[13] =  1.0;     n[13] = -0.21841717175414E-1;
	I[14] = 1.0;    J[14] =  3.0;     n[14] = -0.52838357969930E-4;
	I[15] = 2.0;    J[15] = -3.0;     n[15] = -0.47184321073267E-3;
	I[16] = 2.0;    J[16] =  0.0;     n[16] = -0.30001780793026E-3;
	I[17] = 2.0;    J[17] =  1.0;     n[17] =  0.47661393906987E-4;
	I[18] = 2.0;    J[18] =  3.0;     n[18] = -0.44141845330846E-5;
	I[19] = 2.0;    J[19] = 17.0;     n[19] = -0.72694996295794E-15;
	I[20] = 3.0;    J[20] = -4.0;     n[20] = -0.31679644845054E-4;
	I[21] = 3.0;    J[21] =  0.0;     n[21] = -0.28270797985312E-5;
	I[22] = 3.0;    J[22] =  6.0;     n[22] = -0.85205128120103E-9;
	I[23] = 4.0;    J[23] = -5.0;     n[23] = -0.22425281908000E-5;
	I[24] = 4.0;    J[24] = -2.0;     n[24] = -0.65171222895601E-6;
	I[25] = 4.0;    J[25] = 10.0;     n[25] = -0.14341729937924E-12;
	I[26] = 5.0;    J[26] = -8.0;     n[26] = -0.40516996860117E-6;
	I[27] = 8.0;    J[27] = -11.0;    n[27] = -0.12734301741641E-8;
	I[28] = 8.0;    J[28] = -6.0;     n[28] = -0.17424871230634E-9;
	I[29] = 21.0;   J[29] = -29.0;    n[29] = -0.68762131295531E-18;
	I[30] = 23.0;   J[30] = -31.0;    n[30] =  0.14478307828521E-19;
	I[31] = 29.0;   J[31] = -38.0;    n[31] =  0.26335781662795E-22;
	I[32] = 30.0;   J[32] = -39.0;    n[32] = -0.11947622640071E-22;
	I[33] = 31.0;   J[33] = -40.0;    n[33] =  0.18228094581404E-23;
	I[34] = 32.0;   J[34] = -41.0;    n[34] = -0.93537087292458E-25;

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
		
		gam_tau = gamma_tau (tau, pw);

		h = tau * gam_tau * R * Temp;	/* J/kg */
			
		return(h);
	}


	double gamma_tau (double tau, double pw) const
	{
	        int i;
	        double I[35], J[35];
	        double n[35];
	        double p_star;
	        double pi;
	        double gam_tau;


		/* Regularisierung */
	        if (pw < 1.0) pw = 1.0;
	        if (pw > 1.0E8) pw = 1.0E8;


	        p_star = 1.653E7;        /* [Pa] */

	        pi = pw / p_star;    /* reduced pressure */

	/* list of the numerical values of the coefficients and exponents */

	I[1] = 0.0;     J[1] = -2.0;      n[1] =  0.14632971213167;
	I[2] = 0.0;     J[2] = -1.0;      n[2] = -0.84548187169114;
	I[3] = 0.0;     J[3] =  0.0;      n[3] = -0.37563603672040E1;
	I[4] = 0.0;     J[4] =  1.0;      n[4] =  0.33855169168385E1;
	I[5] = 0.0;     J[5] =  2.0;      n[5] = -0.95791963387872;
	I[6] = 0.0;     J[6] =  3.0;      n[6] =  0.15772038513228;
	I[7] = 0.0;     J[7] =  4.0;      n[7] = -0.16616417199507E-1;
	I[8] = 0.0;     J[8] =  5.0;      n[8] =  0.81214629983568E-3;
	I[9] = 1.0;     J[9] = -9.0;      n[9] =  0.28319080123804E-3;
	I[10] = 1.0;    J[10] = -7.0;     n[10] = -0.60706301565874E-3;
	I[11] = 1.0;    J[11] = -1.0;     n[11] = -0.18990068218419E-1;
	I[12] = 1.0;    J[12] =  0.0;     n[12] = -0.32529748770505E-1;
	I[13] = 1.0;    J[13] =  1.0;     n[13] = -0.21841717175414E-1;
	I[14] = 1.0;    J[14] =  3.0;     n[14] = -0.52838357969930E-4;
	I[15] = 2.0;    J[15] = -3.0;     n[15] = -0.47184321073267E-3;
	I[16] = 2.0;    J[16] =  0.0;     n[16] = -0.30001780793026E-3;
	I[17] = 2.0;    J[17] =  1.0;     n[17] =  0.47661393906987E-4;
	I[18] = 2.0;    J[18] =  3.0;     n[18] = -0.44141845330846E-5;
	I[19] = 2.0;    J[19] = 17.0;     n[19] = -0.72694996295794E-15;
	I[20] = 3.0;    J[20] = -4.0;     n[20] = -0.31679644845054E-4;
	I[21] = 3.0;    J[21] =  0.0;     n[21] = -0.28270797985312E-5;
	I[22] = 3.0;    J[22] =  6.0;     n[22] = -0.85205128120103E-9;
	I[23] = 4.0;    J[23] = -5.0;     n[23] = -0.22425281908000E-5;
	I[24] = 4.0;    J[24] = -2.0;     n[24] = -0.65171222895601E-6;
	I[25] = 4.0;    J[25] = 10.0;     n[25] = -0.14341729937924E-12;
	I[26] = 5.0;    J[26] = -8.0;     n[26] = -0.40516996860117E-6;
	I[27] = 8.0;    J[27] = -11.0;    n[27] = -0.12734301741641E-8;
	I[28] = 8.0;    J[28] = -6.0;     n[28] = -0.17424871230634E-9;
	I[29] = 21.0;   J[29] = -29.0;    n[29] = -0.68762131295531E-18;
	I[30] = 23.0;   J[30] = -31.0;    n[30] =  0.14478307828521E-19;
	I[31] = 29.0;   J[31] = -38.0;    n[31] =  0.26335781662795E-22;
	I[32] = 30.0;   J[32] = -39.0;    n[32] = -0.11947622640071E-22;
	I[33] = 31.0;   J[33] = -40.0;    n[33] =  0.18228094581404E-23;
	I[34] = 32.0;   J[34] = -41.0;    n[34] = -0.93537087292458E-25;

	        gam_tau = 0.0;

	        for (i = 1; i <= 34; i++)
	        {
	                gam_tau = gam_tau + n[i] * pow((7.1 - pi), I[i]) *J[i]*
			pow((tau - 1.222), (J[i]-1));
	        }

	        return(gam_tau);
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
