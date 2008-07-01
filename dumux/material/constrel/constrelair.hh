//from MUFTE
class ConstrelAir
{
	
public:
		
	  /** @brief mass density of gas phase, calculated with ideal gas law
	   * @param T Temperature \f$ \left[ K \right] \f$
	   * @param p Pressure \f$ \left[ Pa \right] \f$
	   * @param  x \f$ \left[ - \right] \f$
	   * @return mass density \f$ \left[ \frac{kg}{m^3} \right] \f$
	   */
	double rho_idGG_molar (double Temp, double pg) const
	{
	    double value;
	    const double RU = 8.31451; 			// univ. gas constant [J/(mol K)] 
	    const double EPS = 1e-18;
	    
	    if(Temp<250.) Temp=250.;    // ACHTUNG Regularisierung
	    if(Temp>500.) Temp=500.;    // ACHTUNG Regularisierung
	    if(pg<EPS) pg=EPS;          // ACHTUNG Regularisierung
	    if(pg>1.E8) pg=1.E8;        // ACHTUNG Regularisierung
	
	    value=pg/(RU*Temp);
	
	    return(value);
	} 

	double rho_idGG_mass (double Temp, double pg, const double molarMass) const
	{
	    double value;
	    const double RU = 8.31451; 			// univ. gas constant [J/(mol K)] 
//	    const double EPS = 1e-18;
	    
//	    if(Temp<250.) Temp=250.;    // ACHTUNG Regularisierung
//	    if(Temp>500.) Temp=500.;    // ACHTUNG Regularisierung
//	    if(pg<EPS) pg=EPS;          // ACHTUNG Regularisierung
//	    if(pg>1.E8) pg=1.E8;        // ACHTUNG Regularisierung
	
	    value=pg/(RU*Temp);
	    value*=molarMass;
	
	    return(value);
	} 

	
	/*** viscosity of air ***/
	double viscosity_air (double Temp) const
	{
        double r;

        if(Temp<250.) Temp = 250.;      /* ACHTUNG Regularisierung */
        if(Temp>500.) Temp = 500.;      /* ACHTUNG Regularisierung */
  
        r = 1.496*1.E-6*pow(Temp,1.5)/(Temp+120.);

        return (r);
	}
	
	/*** saturation pressure curve ***/
	double pwsat(double Temp) const
	{
	    double K_0,K_1,K_2,K_3;
	    const double k0=-4.05968210;
	    const double k1=5.13225550;
	    const double k2=-1.18424070;
	    const double k3=0.117795920;
	    const double k4=-0.00515764200;
	    const double k5=-0.00146895370;
	    const double k6=0.000536228180;
	    const double k7=0.000124553990;
	    const double k8=-0.0000491542880;
	    const double k9=0.0000463025650;
	    const double k10=0.0000153013340;
	    const double k11=-0.0000209545300;
	    double u,T0,T1,T2,T3,T4,T5,T6,T7,T8,T9,T10,T11,theta;
	    double p_s,Temp_c1,p_c1,hoch;
	    const double EPS = 1e-18;
	
	    K_0=2.0;
	    K_1=0.95;
	    K_2=1.45220712;
	    K_3=-0.848789530;
	
	    Temp_c1=647.3;           /* kritische Temperatur [Kelvin] */
	    p_c1=22120000;           /* kritischer Druck [Pa] */
	
	
	    theta=Temp/Temp_c1;     
	    if(theta<0.4)    theta=0.4;     /* ACHTUNG Regularisierung */
	    if(theta>1.05)    theta=1.05;    /* ACHTUNG Regularisierung */
	
	    u=(K_0*pow((1/(theta+EPS) - K_1),0.4) - K_2)/K_3;
	    if(u<-10.0) u=-10.0;     if(u>10.0) u=10.0; /* ACHTUNG Regularisierung */
	
	    T0=1;
	    T1=u;
	    T2=2*u*T1-T0;
	    T3=2*u*T2-T1;
	    T4=2*u*T3-T2;
	    T5=2*u*T4-T3;
	    T6=2*u*T5-T4;
	    T7=2*u*T6-T5;
	    T8=2*u*T7-T6;
	    T9=2*u*T8-T7;
	    T10=2*u*T9-T8;
	    T11=2*u*T10-T9;
	
	    hoch=k0*T0+k1*T1+k2*T2+k3*T3+k4*T4+k5*T5+k6*T6+k7*T7+k8*T8+k9*T9+k10*T10+k11*T11;
	    if(hoch<-40.0) hoch=-40.0;      /* ACHTUNG Regularisierung */
	    if(hoch>40.0) hoch=40.0;
	
	    p_s=exp(hoch)*p_c1;
	
	
//	    if(isnan(p_s))
//	    {
//	        sprintf(buf,"isnan pwsat \n");
//	        p_s = 0.0;
//	    }
	
	    return (p_s);
	
	}
	
	/*** saturation pressure curve ***/
	double pwsat_antoine (double Temp) const
	{
		/*** SOURCE: Luedecke, Luedecke (2000): Thermodynamik, ***/
		/*** Physikalisch-chemische Grundlagen der thermischen ***/
		/*** Verfahrenstechnik; Berlin, Springer               ***/
	
		/* Antoine Equation: log10 (p) = A - B / (C + T)         */
		/*                   units: p in mbar, t in Grad Celsius */
		/*                   validity: 1 <= T <= 100             */
		/* Antoine constants A, B, C after Gmehling et al. 1977  */

		/* Antoine constants */
		const double A = 8.19621;
		const double B = 1730.63;
		const double C = 233.436;		
		double Celsius;		/* temperature in Grad Celsius */
		double exponent;	/* log10 (p) = exponent = A-B/(C+T) */
		double pwsat;
	
		if (Temp < 274.15) Temp = 274.15; /* ACHTUNG Regularisierung */
		if (Temp > 373.15) Temp = 373.15; /* ACHTUNG Regularisierung */
			
		Celsius = Temp - 273.15;
		
		exponent = A - (B / (Celsius + C));
		
		pwsat = pow (10.0, exponent) *100; /*1mbar = 100Pa*/
		
		return(pwsat);
	}
	
	/*** Partialdruck von Wasserdampf in der Gasphase ***/
	double partial_pressure_gas_w (double pg, double Temp, double Sw, double Xwg) const
	{
	    double pwg;
	    const double EPS = 1e-18;
	
	    if(Temp<250.) Temp = 250.;    /* ACHTUNG Regularisierung */
	    if(Temp>500.) Temp = 500.;    /* ACHTUNG Regularisierung */
	    if(pg>1.E8)   pg=1.E8;        /* ACHTUNG Regularisierung */
	    if(pg<1.E-8)  pg=1.E-8;       /* ACHTUNG Regularisierung */
	
	    if (Sw<EPS)
	    {
	        pwg=Xwg*pg;   /* pwgS0 */
	    }
	    else pwg=pwsat(Temp);   /* pwgS */
	
	
//	    if(isnan(pwg)){
//	        sprintf(buf,"isnan partial_pressure_gas_w \n");
//	        pwg = 0.0;
//	    }
	
	    return(pwg);
	}
	
	/*** viscosity of water vapour ***/
	double visco_w_vap (double Temp) const
	{
	    double Tc,Pc,Zc,M,Tr,mu_r;
	    double Fp0,xi,eta_xi,r;
	    const double EPS=1e-18;
	
	    if(Temp<250.) Temp = 250.;    /* ACHTUNG Regularisierung */
	    if(Temp>500.) Temp = 500.;    /* ACHTUNG Regularisierung */
	
	    Tc = 647.3;            /* [K] */
	    Pc = 221.2;            /* [bar] */
	    Zc = 0.231;
	    M = 18.015;            /* [g/mol] */
	    Tr = Temp/Tc;
	    if(Tr<EPS) Tr=EPS;     /* Regularisierung */
	    mu_r = 0.0897;         /* polarity factor */
	    Fp0 = 1. + 0.221*(0.96+0.1*(Tr-0.7));
	    xi = 3.334E-3;
	    eta_xi = (0.807*pow(Tr,0.618) - 0.357*exp((-0.449)*Tr)
	            + 0.34*exp((-4.058)*Tr) + 0.018) * Fp0;
	    r = eta_xi/xi;         /* [1.E-6 P] */
	    r = r/1.E7;            /* [Pa s] */
	
//	    if(isnan(r)){
//	        sprintf(buf,"isnan visco_w_vap \n");
//	        r = 0.0;
//	    }
	
	    return (r);
	}
	
	
	/*** Enthalpie von stark ueberhitztem Dampf ***/
	double hsteam(double p_gw, double Temp) const
	{
	    double A[6];        /* Data A */
	    double C[7];        /* Data C */
	    double B[3][8];     /* Data B */
	    double SB[5];       /* Data SB */
	    double Z[3][8];     /* Data Z */
	    double TK,PK,EL1,X,BEL,BELP,SUM,SC,d__1,S2,v,HK,D,H;
	    double S[3],R[3];
	    double b,z;
	    const double EPS = 1e-18;
	
	    if(Temp<270.)  Temp=270.;        /* ACHTUNG Regularisierung */
	    if(Temp>500.)  Temp=500.;        /* ACHTUNG Regularisierung */
	    if(p_gw<1.E-8) p_gw=1.E-8;       /* ACHTUNG Regularisierung */
	    if(p_gw>1.E8)  p_gw=1.E8;        /* ACHTUNG Regularisierung */
	
	    /* Initialisierung der Data-Arrays */
	    A[0]=16.83599274;
	    A[1]=28.56067796;
	    A[2]=-54.38923329;
	    A[3]=0.4330662834;
	    A[4]=-0.6547711697;
	    A[5]=8.565182058E-2;
	    C[0]=1.936587558E2;
	    C[1]=-1.388522425E3;
	    C[2]=4.126607219E3;
	    C[3]=-6.508211677E3;
	    C[4]=5.745984054E3;
	    C[5]=-2.693088365E3;
	    C[6]=5.235718623E2;
	    B[0][0]=6.670375918E-2;
	    B[0][1]=1.388983801;
	    B[0][2]=0.;
	    B[0][3]=8.390104328E-2;
	    B[0][4]=2.614670893E-2;
	    B[0][5]=-3.373439453E-2;
	    B[0][6]=4.520918904E-1;
	    B[0][7]=1.069036614E-1;
	    B[1][0]=0.;
	    B[1][1]=-5.975336707E-1;
	    B[1][2]=-8.847535804E-2;
	    B[1][3]=0.;
	    B[1][4]=5.958051609E-1;
	    B[1][5]=-5.159303373E-1;
	    B[1][6]=2.075021122E-1;
	    B[1][7]=1.190610271E-1;
	    B[2][0]=-9.867174132E-2;
	    B[2][1]=0.;
	    B[2][2]=1.683998803E-1;
	    B[2][3]=-5.809438001E-2;
	    B[2][4]=0.;
	    B[2][5]=6.552390126E-3;
	    B[2][6]=5.710218649E-4;
	    B[2][7]=0.;
	    SB[0]=7.633333333E-1;
	    SB[1]=4.006073948E-1;
	    SB[2]=8.636081627E-2;
	    SB[3]=-8.532322921E-1;
	    SB[4]=3.460208861E-1;
	    Z[0][0]=13.;
	    Z[0][1]=3.;
	    Z[0][2]=0.;
	    Z[0][3]=18.;
	    Z[0][4]=2.;
	    Z[0][5]=1.;
	    Z[0][6]=18.;
	    Z[0][7]=10.;
	    Z[1][0]=0.;
	    Z[1][1]=25.;
	    Z[1][2]=14.;
	    Z[1][3]=0.;
	    Z[1][4]=32.;
	    Z[1][5]=28.;
	    Z[1][6]=24.;
	    Z[1][7]=12.;
	    Z[2][0]=11.;
	    Z[2][1]=0.;
	    Z[2][2]=24.;
	    Z[2][3]=18.;
	    Z[2][4]=0.;
	    Z[2][5]=24.;
	    Z[2][6]=14.;
	    Z[2][7]=0.;
	
	    /* Berechnung */
	    TK=Temp/647.3;
	    PK=p_gw/2.212E7;
	    EL1=4.260321148;
	    X=exp(SB[0]*(1-TK));
	    BEL=15.74373327-34.17061978*TK+19.31380707*TK*TK;
	    BELP=-34.17061978+2*19.31380707*TK;
	    SUM=EL1*(TK/PK);
	
	    for (int j=1; j<=5; ++j)
	    {
	        SC=0.;
	        for (int k=1; k<=3; ++k)
	        {
	            b=B[k-1][j-1]; 
	            z=Z[k-1][j-1]; 
	            SC=b*pow(X,z)+SC;
	        }
	        d__1=(double) (j-1);
	        SUM=SUM-j*pow(PK,d__1)*SC;
	    }
	
	    S[0]=1/(pow(PK,4)+EPS) + SB[1]*pow(X,14);
	    S[1]=1/(pow(PK,5)+EPS) + SB[2]*pow(X,19);
	    S[2]=1/(pow(PK,6)+EPS) + SB[3]*pow(X,54) + SB[4]*pow(X,27);
	
	    for (int j=6; j<=8; ++j)
	    {
	        S2=0.;
	        for (int k=1; k<=3; ++k)
	        {
	            S2=S2+B[k-1][j-1]*pow(X,(Z[k-1][j-1]));
	        }
	        d__1= (double) (1-j);
	        SUM=SUM-(S2/pow(S[j-6],2))*(j-2)*pow(PK,d__1);
	    }
	    S2=0.; 
	
	    for (int k=1; k<=7; ++k)
	    {
	        d__1= (double) (k-1);
	        S2=C[k-1]*pow(X,d__1)+S2;
	    }
	
	    SUM=SUM+11*pow((PK/BEL),10)*S2;
	    v=SUM*0.00317;            /* specific volume */
	    D=1/v;                /* density !! */
	
	    HK=0.;
	    HK=A[0]*TK;
	
	    for (int k=1; k<=5; ++k)
	    {
	        d__1= (double) (k-1);
	        HK=HK-A[k]*(k-2)*pow(TK,d__1);
	    }
	
	    for (int j=1; j<=5; ++j)
	    {
	        S2=0.;
	        for (int k=1; k<=3; ++k)
	        {
	            S2=S2+B[k-1][j-1]*(1+Z[k-1][j-1]*SB[0]*TK)*pow(X,(Z[k-1][j-1]));
	        }
	        d__1= (double) j;
	        HK=HK-pow(PK,d__1)*S2;
	    }
	
	    R[0]=14*SB[1]*pow(X,14);
	    R[1]=19*SB[2]*pow(X,19);
	    R[2]=54*SB[3]*pow(X,54)+27*SB[4]*pow(X,27);
	
	    for (int j=6; j<=8; ++j)
	    {
	        S2=0.;
	        for (int k=1; k<=3; ++k)
	        {
	            S2=S2+B[k-1][j-1]*pow(X,(Z[k-1][j-1]))*(1/S[j-6])*
	                   (1+(Z[k-1][j-1]*SB[0]*TK-SB[0]*TK*R[j-6])/S[j-6]);
	        }
	        HK=HK-S2;
	    }
	    S2=0.;
	
	    for (int k=1; k<=7; ++k)
	    {
	        d__1= (double) (k-1);
	        S2=S2+(1+TK*(10*(BELP/BEL)+(k-1)*SB[0]))*C[k-1]*pow(X,d__1);
	    }
	
	    HK=HK+PK*pow((PK/BEL),10)*S2;
	
	    H=HK*70120.4;          /* specific enthalpy */
	
	    /* U=H-p_gw*v; */            /* specific Internal energy */
		
//	    if(isnan(H)){
//	        sprintf(buf,"isnan hsteam \n");
//	        H = 0.0;
//	    }
	
	    return (H);
	}
	
	
	/* Interpolationspolynome fuer leicht ueberhitzten */
	/* Dampf fuer 100/50/25/5/1 bar                    */ 
	double h_supst100(double T) const
	{
	    double h,x;
	    x=T/100;
	
	    h=(((-115262.537121*x+1865074.74575)*x-11287756.0179)*x+30630020.8681)*x-28677461.4108;
	    /* h=-28677461.4108+30630020.8681*x-11287756.0179*pow(x,2)
	         +1865074.74575*pow(x,3)-115262.537121*pow(x,4); */
		
	    return(h);
	}
	
	//	/* Enthalpie von Wasserdampf in der Gasphase */
	//	double enth_gw(double Sw, double p_gw, double Temp)
	//	{
	//	    double T_sat=0.,epsilon,tsa,tsa_neu,g,g_strich,Celsius;
	//	    double pwsat_tsa;
	//	    int i;
	//	    int crit_I,crit_II=0.;
	//	    double a,h;
	//	    double h0,DT_sat,Dh1,Dh2,Dh;
	//	    double wicht;
	//	    const double EPS = 1e-18;
	//	
	//	    if(Temp<250.) Temp=250.;    /* ACHTUNG Regularisierung */
	//	    if(Temp>500.) Temp=500.;    /* ACHTUNG Regularisierung */
	//	    if(p_gw<EPS) p_gw=EPS;    /* ACHTUNG Regularisierung */
	//	    if(p_gw>1.E8) p_gw=1.E8;     /* ACHTUNG Regularisierung */
	//	    if(Sw<0.0) Sw=0.0;        /* ACHTUNG Regularisierung */
	//	    if(Sw>1.0) Sw=1.0;        /* ACHTUNG Regularisierung */
	//	
	//	    Celsius=Temp-273.15;
	//	
	//	    // Nassdampf, leicht oder stark ueberhitzter Dampf ??
	//	
	//	    if (Sw>EPS) crit_I=1;   /* Nassdampf ! */
	//	    else
	//	    {
	//	        /*** Ermittlung von T_sat(p_gw) ***/
	//			/* Antoine-Gleichung ist hinreichend genau, */
	//			/* auch wenn pwsat als Funktion gewaehlt wurde */
	//
	////			T_sat = inverse_pwsat_antoine(p_gw);
	//	
	//			/*** Ermittlung von T_sat(p_gw) ***/
	//			/* Suche T_sat iterativ: 					*/
	//			/* T[i+1] = T[i]+ pwsat(T[i])/pwsat'(T[i])	*/
	//			/* Bestimme pwsat'(T) numerisch				*/
	//			/*********************************************
	//	        tsa=Temp;
	//	        epsilon=1.E-2;
	//	        for(i=1;i<100;++i)
	//	        {
	//	            pwsat_tsa=pwsat(tsa);
	//	            g=pwsat_tsa-p_gw;
	//	            g_strich=(pwsat(tsa+epsilon)-pwsat_tsa)/epsilon;
	//	            if (g_strich <= epsilon) g_strich = epsilon;
	//	            tsa_neu=tsa-g/g_strich;
	//	            a=tsa-tsa_neu;
	//	            if (a<0) a=-a;
	//	            if (a<epsilon)
	//	            {
	//	              T_sat=tsa;
	//	              break;
	//	            }
	//	            tsa=tsa_neu;
	//	        }
	//			*********************************************/
	//	
	//	        /*if (Celsius/(T_sat-273.15)<1.8) crit_I=2;*/   /* leicht ueberhitzt */
	//	        if (Celsius < 1.8*(T_sat-273.15)) crit_I=2;   /* leicht ueberhitzt */
	//	        else crit_I=3;   /* stark ueberhitzt */
	//	    }
	//	
	//	    switch(crit_I)
	//	    {
	//	        case 1:
	//	        {
	//	            h=hsat(Temp);
	//	
	////	            if(isnan(h)){
	////	            sprintf(buf,"isnan enth_gw \n");
	////	            h = 0.0;
	////	            }
	//	            
	//	            return (h);
	//	        }
	//	        case 2:
	//	        {
	//	            /* Abfrage des Interpolationsbereichs */
	//	            if (p_gw<=1.E5) crit_II=1;
	//	            if ((p_gw>1.E5)&&(p_gw<=5.E5)) crit_II=2;
	//	            if ((p_gw>5.E5)&&(p_gw<=25.E5)) crit_II=3;
	//	            if ((p_gw>25.E5)&&(p_gw<=50.E5)) crit_II=4;
	//	            if ((p_gw>50.E5)&&(p_gw<=100.E5)) crit_II=5;
	//	
	//	            switch(crit_II)
	//	            {
	//	                case 1:
	//	                {
	//	                    h0=hsat(T_sat);
	//	                    h=h0-1.759*(T_sat-Temp);
	//	                    wicht=(1.E5-p_gw)/1.E5;
	//	                    DT_sat=Temp-T_sat;
	//	                    Dh1=h_supst1(DT_sat+99.632)-h_supst1(99.632);
	//	                    h=h*wicht+(1-wicht)*(h0+Dh1);
	//	
	//	                    return (h);
	//	                  }
	//	                  case 2:
	//	                  {
	//	                    h0=hsat(T_sat);
	//	                    DT_sat=Temp-T_sat;
	//	                    Dh1=h_supst1(DT_sat+99.632)-h_supst1(99.632);
	//	                    Dh2=h_supst5(DT_sat+151.866)-h_supst5(151.866);
	//	                    Dh=(T_sat-273.15-99.632)*Dh2/52.231  
	//	                        + (151.866-(T_sat-273.15))*Dh1/52.231;
	//	                    h=h0+Dh;
	//	
	//	                    return (h);
	//	                  }
	//	                  case 3:
	//	                  {
	//	                    h0=hsat(T_sat);
	//	                    DT_sat=Temp-T_sat;
	//	                    Dh1=h_supst5(DT_sat+151.866)-h_supst5(151.866);
	//	                    Dh2=h_supst25(DT_sat+223.989)-h_supst25(223.989);
	//	                    Dh=(T_sat-273.15-151.866)*Dh2/72.123 
	//	                        + (223.989-(T_sat-273.15))*Dh1/72.123;
	//	                    h=h0+Dh;
	//	
	//	                    return (h);
	//	                  }    
	//	                  case 4:
	//	                  {
	//	                    h0=hsat(T_sat);
	//	                    DT_sat=Temp-T_sat;
	//	                    Dh1=h_supst25(DT_sat+223.989)-h_supst25(223.989);
	//	                    Dh2=h_supst50(DT_sat+263.977)-h_supst50(263.977);
	//	                    Dh=(T_sat-273.15-223.989)*Dh2/39.988 
	//	                        + (263.977-(T_sat-273.15))*Dh1/39.988;
	//	                    h=h0+Dh;
	//	
	//	                    return (h);
	//	                  }
	//	                  case 5:
	//	                  {
	//	                    h0=hsat(T_sat);
	//	                    DT_sat=Temp-T_sat;
	//	                    Dh1=h_supst50(DT_sat+263.977)-h_supst50(263.977);
	//	                    Dh2=h_supst100(DT_sat+311.031)-h_supst100(311.031);
	//	                    Dh=(T_sat-273.15-263.977)*Dh2/47.054 
	//	                        + (311.031-(T_sat-273.15))*Dh1/47.054;
	//	                    h=h0+Dh;
	//	
	//	                    return (h);
	//	                  }
	//	            }
	//	        }
	//	        case 3:
	//	        {
	//	            h=hsteam(p_gw,Temp);
	//	
	//	
	////	            if(isnan(h)){
	////	                sprintf(buf,"isnan enth_gw \n");
	////	                h = 0.0;
	//	            }
	//	
	//	              return (h);
	//	        }
	//	    }
	//	    return(h); /* just for the compiler */
	//	}
		
	//	/*** spezifische Enthalpie der Gasphase EINHEITEN ! ***/
	//	double sp_enth2p2cni_g (double Temp, double Xag, double Xwg, double
	//	pg, double Sw)
	//	{
	//	   double omega_gw, omega_ag; /* mass fractions (notwendig?) */
	//	   double h_g, h_ga, h_gw, C_va, p_gw, m_w;
	//	   const double EPS)=1e-18;
	//	
	//	   if(Temp<EPS)  Temp=EPS;      /* ACHTUNG Regularisierung */
	//	   if(Temp>500.) Temp=500.;     /* ACHTUNG Regularisierung */
	//	   if(pg<EPS)    pg=EPS;        /* ACHTUNG Regularisierung */
	//	   if(pg>1.E8)   pg=1.E8;       /* ACHTUNG Regularisierung */
	//	   if(Sw<0.0)    Sw = 0.0;      /* ACHTUNG Regularisierung */
	//	   if(Sw>1.)     Sw = 1.;       /* ACHTUNG Regularisierung */
	//	   if(Xwg<0.0)   Xwg = 0.0;     /* ACHTUNG Regularisierung */
	//	   if(Xwg>1.)    Xwg = 1.;      /* ACHTUNG Regularisierung */
	//	   if(Xag<0.0)   Xag = 0.0;     /* ACHTUNG Regularisierung */
	//	   if(Xag>1.)    Xag = 1.;      /* ACHTUNG Regularisierung */
	//	
	//	   omega_ag = (Xag*molecular_weight_air())/(Xag*molecular_weight_air()
	//	               + Xwg*molecular_weight_water());
	//	   omega_gw = (Xwg*molecular_weight_water())/(Xag*molecular_weight_air()
	//	               + Xwg*molecular_weight_water());
	//	
	//	   p_gw = pg*Xwg;
	//	   h_gw = enth_gw(Sw,p_gw,Temp);
	//	   C_va = 733.;   /* W"armekapazit"at von Luft [J/(kg * K)] */
	//	   h_ga = C_va*(Temp-273.15) + Rluft*Temp; 
	//	   h_g = omega_ag*h_ga+omega_gw*h_gw;
	//	
	//	
	//	    if(isnan(h_g)){
	//	        sprintf(buf,"isnan sp_enthalpy_g \n");
	//	        h_g = 0.0;
	//	    }
	//	   return(h_g);
	//	}


};


