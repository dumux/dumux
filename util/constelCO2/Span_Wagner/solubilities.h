/************************************************************************/
/*									*/
/*	In this file the functions for the CO2-solubility in brine	*/
/*	and water solubility in the CO2-phase are implemented after	*/
/*	Duan & Sun, Chemical Geology 193 (2003) 257-271			*/
/*									*/
/************************************************************************/

#ifndef DOUBLE
#define DOUBLE double
#define INT int
#endif

DOUBLE XCO2inGasphase(DOUBLE Temp, DOUBLE pg);
DOUBLE pure_water_pressure (DOUBLE Temp);

DOUBLE ComputeA(DOUBLE Temp, DOUBLE pg);
DOUBLE ComputeB(DOUBLE Temp, DOUBLE pg);
DOUBLE ComputeC(DOUBLE Temp, DOUBLE pg);
DOUBLE ComputeZ(DOUBLE Temp, DOUBLE pg);
DOUBLE FugacityCO2(DOUBLE Temp, DOUBLE pg);

DOUBLE SolCO2inWater(DOUBLE Temp, DOUBLE pg, DOUBLE X_NaCl);

double iterateVr(DOUBLE T, DOUBLE p);

DOUBLE HeatOfSolutionCO2inWater (DOUBLE Temp);

DOUBLE Henry_coeff (DOUBLE Temp, DOUBLE molality);
DOUBLE Henry_Bando(DOUBLE Temp, DOUBLE p, DOUBLE Xsw);



