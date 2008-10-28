// $Id$

#ifndef DUNE_CONSTRELBRINE_HH
#define DUNE_CONSTRELBRINE_HH

class ConstrelBrine {
public:
	/**************************************************************************/
	/*                                                                        */
	/* Computation of the mass density of brine Batzle & Wang (1992)          */
	/* equation given by Adams & Bachu in Geofluids (2002) 2, 257-271         */
	/*                                                                        */
	/**************************************************************************/

	double mass_density_brine(double Temp, double pw, double Xsw) const {
		/* Xsw : mass fraction of salt in water [-] */

		double rhoBrine, rhow;
		double TempC, pMPa;
		ConstrelWater water;

		TempC = Temp - 273.15;
		pMPa = pw/1.0E6;

		rhow = water.mass_density_water_Batzle (Temp, pw);

		rhoBrine = rhow + Xsw*(0.668 + 0.44*Xsw + 1.0E-6
				*(300*pMPa - 2400*pMPa*Xsw+ TempC*(80.0 - 3*TempC - 3300*Xsw
						- 13*pMPa + 47*pMPa*Xsw)))*1000.;

		return (rhoBrine); /* unit: [kg/m^3] */
	}
	/***********************************************************************/
	/*                                                                     */
	/* Total brine density with dissolved CO2                              */
	/* rho_{b,CO2} = rho_w + contribution(salt) + contribution(CO2)        */
	/*                                                                     */
	/***********************************************************************/

	double mass_density_brine_CO2(double Temp, double pw, double S,
			double x_CO2_w) const {
		double rho_brine, rho_pure, rho_wCO2, contribCO2, rho_brineCO2;
		ConstrelWater water;
		/* S : salinity as a mass fraction [-] */
		/* x_CO2_w : mole fraction of CO2 in water phase */

		rho_brine = this->mass_density_brine (Temp, pw, S);
		rho_pure = water.mass_density_water_Batzle (Temp, pw);
		rho_wCO2 = water.mass_density_waterCO2 (Temp, pw, x_CO2_w);
		contribCO2 = rho_wCO2 - rho_pure;

		rho_brineCO2 = rho_brine + contribCO2;

		return (rho_brineCO2);
	}

	/******************************************************************/
	/*                                                                */
	/* Viscosity of brine after Batzle & Wang (1992)                  */
	/* Cited by Bachu and Adams: "Equations of State for basin        */
	/* geofluids" (2002)                                              */
	/*                                                                */
	/******************************************************************/

	double viscosity_brine(double Temp, double S) const {
		double T_C;
		double A;
		double mu_brine;

		if(Temp <= 275.) // regularisation
		{
			Temp = 275;
		}

		T_C = Temp - 273.15;

		A = (0.42*pow((pow(S, 0.8)-0.17), 2) + 0.045)*pow(T_C, 0.8);

		mu_brine = 0.1 + 0.333*S + (1.65+91.9*S*S*S)*exp(-A);

		mu_brine = mu_brine/1000.0; /* unit: Pa s */

		return mu_brine;
	}

	double enthalpy_brine(double T, double p, double X) const{

		double theta, h_NaCl;
		double f[4], a[4][3], m, h_ls, h_ls1, d_h;
		double X_lSAT, delta_h;
		int i, j;
		double R, hw;
		ConstrelWater water;

		R = 8.314E-3; /*kJ/(mol*K)*/
		theta = T - 273.15;

		/*Numerische Koeffizienten nach PALLISER*/
		f[0] = 2.63500E-1;
		f[1] = 7.48368E-6;
		f[2] = 1.44611E-6;
		f[3] = -3.80860E-10;

		/*Numerische Koeffizienten nach MICHAELIDES f"ur Enthalpie Brine*/
		a[0][0] = -9633.6;
		a[2][0] = -0.90963;
		a[0][1] = -4080.0;
		a[2][1] = -0.36524;
		a[0][2] = +286.49;
		a[2][2] = +0.249667E-1;
		a[1][0] = +166.58;
		a[3][0] = +0.17965E-2;
		a[1][1] = +68.577;
		a[3][1] = +0.71924E-3;
		a[1][2] = -4.6856;
		a[3][2] = -0.4900E-4;

		X_lSAT = f[0] + f[1]*theta + f[2]*pow(theta, 2) + f[3]*pow(theta, 3);
		/*Regularisierung*/
		if (X>X_lSAT) {
			X = X_lSAT;
		}

		hw = water.enthalpy_water (T, p) /1E3; /* kJ/kg */

		/*DAUBERT UND DANNER*/
		/*U=*/h_NaCl = (3.6710E4*T + 0.5*(6.2770E1)*T*T - ((6.6670E-2)/3)*T*T*T
				+((2.8000E-5)/4)*pow(T, 4))/(58.44E3)- 2.045698e+02; /* kJ/kg */

		m = (1E3/58.44)*(X/(1-X));
		i = 0;
		j = 0;
		d_h = 0;

		for (i = 0; i<=3; i++) {
			for (j=0; j<=2; j++) {
				d_h = d_h + a[i][j] * pow(theta, i) * pow(m, j);
			}
		}

		delta_h = (4.184/(1E3 + (58.44 * m)))*d_h;

		/* Enthalpy of brine */

		h_ls1 =(1-X)*hw + X*h_NaCl + X*delta_h; /* kJ/kg */

		h_ls = h_ls1*1E3; /*J/kg*/

		return (h_ls);
	}

	double enthalpy_brine_CO2(double T, double p, double X, double X_CO2_w) const{
		/* X_CO2_w : mass fraction of CO2 in brine */

		/* same function as enthalpy_brine, only extended by CO2 content */

		double theta, h_NaCl;
		double f[4], a[4][3], m, h_ls, h_ls1, d_h;
		double X_lSAT, delta_h;
		int i, j;
		double delta_hCO2, hg, R, hw;
		ConstrelWater water;
		ConstrelCO2 co2;

		R = 8.314E-3; /*kJ/(mol*K)*/
		theta = T - 273.15;

		/*Numerische Koeffizienten nach PALLISER*/
		f[0] = 2.63500E-1;
		f[1] = 7.48368E-6;
		f[2] = 1.44611E-6;
		f[3] = -3.80860E-10;

		/*Numerische Koeffizienten nach MICHAELIDES f"ur Enthalpie Brine*/
		a[0][0] = -9633.6;
		a[2][0] = -0.90963;
		a[0][1] = -4080.0;
		a[2][1] = -0.36524;
		a[0][2] = +286.49;
		a[2][2] = +0.249667E-1;
		a[1][0] = +166.58;
		a[3][0] = +0.17965E-2;
		a[1][1] = +68.577;
		a[3][1] = +0.71924E-3;
		a[1][2] = -4.6856;
		a[3][2] = -0.4900E-4;

		X_lSAT = f[0] + f[1]*theta + f[2]*pow(theta, 2) + f[3]*pow(theta, 3);
		/*Regularisierung*/
		if (X>X_lSAT) {
			X = X_lSAT;
		}

		hw = water.enthalpy_water (T, p) /1E3; /* kJ/kg */

		/*DAUBERT UND DANNER*/
		/*U=*/h_NaCl = (3.6710E4*T + 0.5*(6.2770E1)*T*T - ((6.6670E-2)/3)*T*T*T
				+((2.8000E-5)/4)*pow(T, 4))/(58.44E3)- 2.045698e+02; /* kJ/kg */

		m = (1E3/58.44)*(X/(1-X));
		i = 0;
		j = 0;
		d_h = 0;

		for (i = 0; i<=3; i++) {
			for (j=0; j<=2; j++) {
				d_h = d_h + a[i][j] * pow(theta, i) * pow(m, j);
			}
		}

		delta_h = (4.184/(1E3 + (58.44 * m)))*d_h;

		/* Enthalpy of brine without CO2 */

		h_ls1 =(1-X)*hw + X*h_NaCl + X*delta_h; /* kJ/kg */

		/* Enthalpy of CO2 */
		hg = co2.enthalpy(T, p)/1E3- delta_hCO2;

		/* Enthalpy of brine with dissolved CO2 */
		h_ls = (h_ls1*(1-X_CO2_w)+ hg*X_CO2_w)*1E3; /*J/kg*/

		return (h_ls);
	}
};

#endif
