// $Id$

#ifndef DUNE_CONSTRELAIR_HH
#define DUNE_CONSTRELAIR_HH

#include <dumux/exceptions.hh>

namespace Dune {

/** \todo Please doc me! */

class ConstrelAir
{
static const double eps_ = 1e-18;

public:

      /** @brief molar density of gas phase, calculated with ideal gas law
       * @param temperature Temperature \f$ \left[ K \right] \f$
       * @param p Pressure \f$ \left[ Pa \right] \f$
       * @return molar density \f$ \left[ \frac{mol}{m^3} \right] \f$
       */
    double rho_idGG_molar (double temperature, double pg) const
    {
        double rho;
        const double RU = 8.31451;          // univ. gas constant [J/(mol K)]

        if(temperature < 273.15 || temperature > 500.){
            DUNE_THROW(Dune::NumericalProblem,
                "ConstrelAir: Temperature " << temperature << " out of range at " << __FILE__ << ":" << __LINE__);
        }
        else if (pg < eps_ || pg > 1.0E8) {
            DUNE_THROW(Dune::NumericalProblem,
                "ConstrelAir: Pressure " << pg << " out of range at " << __FILE__ << ":" << __LINE__);
        }

        rho=pg/(RU*temperature);

        return(rho);
    }

      /** @brief mass density of gas phase, calculated with ideal gas law
       * @param temperature Temperature \f$ \left[ K \right] \f$
       * @param p Pressure \f$ \left[ Pa \right] \f$
       * @return mass density \f$ \left[ \frac{kg}{m^3} \right] \f$
       */
    double rho_idGG_mass (double temperature, double pg, const double molarMass) const
    {
        double rho;
        const double RU = 8.31451;          // univ. gas constant [J/(mol K)]

        if(temperature < 273.15 || temperature > 500.)
        {
            DUNE_THROW(Dune::NumericalProblem,
                "ConstrelAir: Temperature " << temperature << " out of range at " << __FILE__ << ":" << __LINE__);
        }
        else if (pg < eps_ || pg > 1.0E8) {
            DUNE_THROW(Dune::NumericalProblem,
                "ConstrelAir: Pressure " << pg << " out of range at " << __FILE__ << ":" << __LINE__);
        }

        rho=molarMass*pg/(RU*temperature);

        return(rho);
    }


    /*** viscosity of air ***/
    double viscosity_air (double temperature) const
    {
        double r;

        if(temperature < 273.15 || temperature > 500.)
        {
            DUNE_THROW(Dune::NumericalProblem,
                "ConstrelAir: Temperature " << temperature << " out of range at " << __FILE__ << ":" << __LINE__);
        }

        r = 1.496*1.E-6*pow(temperature,1.5)/(temperature+120.);

        return (r);
    }

    /*** saturation pressure curve ***/
    double pwsat(double temperature) const
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

        K_0=2.0;
        K_1=0.95;
        K_2=1.45220712;
        K_3=-0.848789530;

        Temp_c1=647.3;           /* kritische Temperatur [Kelvin] */
        p_c1=22120000;           /* kritischer Druck [Pa] */


        theta=temperature/Temp_c1;
        if(theta < 0.4 || theta > 1.05)
        {
            DUNE_THROW(Dune::NumericalProblem,
                "ConstrelAir: Theta " << theta << " out of range at " << __FILE__ << ":" << __LINE__);
        }
//        if(theta<0.4)    theta=0.4;     /* ACHTUNG Regularisierung */
//        if(theta>1.05)    theta=1.05;    /* ACHTUNG Regularisierung */

        u=(K_0*pow((1/(theta+eps_) - K_1),0.4) - K_2)/K_3;
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

        return (p_s);

    }

    /*** saturation pressure curve ***/
    double pwsat_antoine (double temperature) const
    {
        /*** SOURCE: Luedecke, Luedecke (2000): Thermodynamik, ***/
        /*** Physikalisch-chemische Grundlagen der thermischen ***/
        /*** Verfahrenstechnik; Berlin, Springer               ***/

        /* Antoine Equation: log10 (p) = A - B / (C + T)         */
        /*                   units: p in mbar, t in Grad Celsius */
        /*                   validity: 1 <= temperature <= 100   */
        /* Antoine constants A, B, C after Gmehling et al. 1977  */

        /* Antoine constants */
        const double A = 8.19621;
        const double B = 1730.63;
        const double C = 233.436;
        double Celsius;     /* temperature in Grad Celsius */
        double exponent;    /* log10 (p) = exponent = A-B/(C+T) */
        double pwsat;

        if(temperature < 274.15 || temperature > 373.15)
        {
            DUNE_THROW(Dune::NumericalProblem,
                "ConstrelAir: Temperature " << temperature << " out of range at " << __FILE__ << ":" << __LINE__);
        }
//        if (temperature < 274.15) temperature = 274.15; /* ACHTUNG Regularisierung */
//        if (temperature > 373.15) temperature = 373.15; /* ACHTUNG Regularisierung */

        Celsius = temperature - 273.15;

        exponent = A - (B / (Celsius + C));

        pwsat = pow (10.0, exponent) *100; /*1mbar = 100Pa*/

        return(pwsat); // [Pa]
    }

    // Partial pressure of water vapor in the gas phase
    double partial_pressure_gas_w (double pg, double temperature, double Sw, double Xwg) const
    {
        double pwg;

        if(temperature < 250. || temperature > 500.)
        {
            DUNE_THROW(Dune::NumericalProblem,
                "ConstrelAir: Temperature " << temperature << " out of range at " << __FILE__ << ":" << __LINE__);
        }
        else if (pg < eps_ || pg > 1.0E8) {
            DUNE_THROW(Dune::NumericalProblem,
                "ConstrelAir: Pressure " << pg << " out of range at " << __FILE__ << ":" << __LINE__);
        }
//        if(temperature<250.) temperature = 250.;    /* ACHTUNG Regularisierung */
//        if(temperature>500.) temperature = 500.;    /* ACHTUNG Regularisierung */
//        if(pg>1.E8)   pg=1.E8;        /* ACHTUNG Regularisierung */
//        if(pg<1.E-8)  pg=1.E-8;       /* ACHTUNG Regularisierung */

        if (Sw<eps_)
        {
            pwg=Xwg*pg;   /* pwgS0 */
        }
        else pwg=pwsat(temperature);   /* pwgS */

        return(pwg);
    }

    // viscosity of water vapour
    double visco_w_vap (double temperature) const
    {
        double Tc,Pc,Zc,M,Tr,mu_r;
        double Fp0,xi,eta_xi,r;

        if(temperature < 250. || temperature > 500.)
        {
            DUNE_THROW(Dune::NumericalProblem,
                "ConstrelAir: Temperature " << temperature << " out of range at " << __FILE__ << ":" << __LINE__);
        }

        Tc = 647.3;            /* [K] */
        Pc = 221.2;            /* [bar] */
        Zc = 0.231;
        M = 18.015;            /* [g/mol] */
        Tr = temperature/Tc;
        if(Tr<eps_) Tr=eps_;     /* Regularisierung */
        mu_r = 0.0897;         /* polarity factor */
        Fp0 = 1. + 0.221*(0.96+0.1*(Tr-0.7));
        xi = 3.334E-3;
        eta_xi = (0.807*pow(Tr,0.618) - 0.357*exp((-0.449)*Tr)
                + 0.34*exp((-4.058)*Tr) + 0.018) * Fp0;
        r = eta_xi/xi;         /* [1.E-6 P] */
        r = r/1.E7;            /* [Pa s] */

//      if(isnan(r)){
//          sprintf(buf,"isnan visco_w_vap \n");
//          r = 0.0;
//      }

        return (r);
    }

    /****************************************************************************/
    /****************************************************************************/
    /*  specific enthalpies and heat capacities (thermo properties)             */
    /****************************************************************************/
    /****************************************************************************/

    /* sources: Falta 1990, IFC 1967   */


    /* specific enthalpy of the water phase [J/kg]  (  -> 1 J = 1 (kg*m2)/s2  ) */
    double sp_enthalpy_w (double temperature, double pw)
    {
    double A[23],SA[12];
    double TKR,PNMR,Y,ZP,Z,PAR1,PAR2,PAR3,PAR4,PAR5;
    double VMKR,v,D,YD,SNUM,PRT1,PRT2,PRT3,PRT4,PRT5,ENTR,H,U;
    int i;

    if(temperature < 250. || temperature > 500.)
    {
        DUNE_THROW(Dune::NumericalProblem,
            "ConstrelAir: Temperature " << temperature << " out of range at " << __FILE__ << ":" << __LINE__);
    }
    else if (pw < eps_ || pw > 1.0E8) {
        DUNE_THROW(Dune::NumericalProblem,
            "ConstrelAir: Pressure " << pw << " out of range at " << __FILE__ << ":" << __LINE__);
    }

    A[0]=6.824687741E3;
    A[1]=-5.422063673E2;
    A[2]=-2.096666205E4;
    A[3]=3.941286787E4;
    A[4]=-6.733277739E4;
    A[5]=9.902381028E4;
    A[6]=-1.093911774E5;
    A[7]=8.590841667E4;
    A[8]=-4.511168742E4;
    A[9]=1.418138926E4;
    A[10]=-2.017271113E3;
    A[11]=7.982692717;
    A[12]=-2.616571843E-2;
    A[13]=1.522411790E-3;
    A[14]=2.284279054E-2;
    A[15]=2.421647003E2;
    A[16]=1.269716088E-10;
    A[17]=2.074838328E-7;
    A[18]=2.17402035E-8;
    A[19]=1.105710498E-9;
    A[20]=1.293441934E1;
    A[21]=1.308119072E-5;
    A[22]=6.047626338E-14;
    SA[0]=8.438375405E-1;
    SA[1]=5.362162162E-4;
    SA[2]=1.720;
    SA[3]=7.342278489E-2;
    SA[4]=4.975858870E-2;
    SA[5]=6.53715430E-1;
    SA[6]=1.15E-6;
    SA[7]=1.5108E-5;
    SA[8]=1.4188E-1;
    SA[9]=7.002753165;
    SA[10]=2.995284926E-4;
    SA[11]=2.04E-1;

            /* Berechnung */
    TKR=temperature/647.3;
    PNMR=pw/2.212E7;
    Y=1-SA[0]*TKR*TKR-SA[1]/pow(TKR,6);
    ZP=(SA[2]*Y*Y-2*SA[3]*TKR+2*SA[4]*PNMR);
    if (ZP<0.) return (-1); /* Temperature out of range */

    Z=Y+pow(ZP,0.5);
    PAR1=A[11]*SA[4]/pow(Z,(5/17));
    PAR2=A[12]+A[13]*TKR+A[14]*TKR*TKR+A[15]*pow((SA[5]-TKR),10)+A[16]/(SA[6]
                                                                 +pow(TKR,19));
    PAR3=(A[17]+2*A[18]*PNMR+3*A[19]*PNMR*PNMR)/(SA[7]+pow(TKR,11));
    PAR4=A[20]*pow(TKR,18)*(SA[8]+TKR*TKR)*(-3/pow((SA[9]+PNMR),4) + SA[10]);
    PAR5=3*A[21]*(SA[11]-TKR)*PNMR*PNMR+4*A[22]*pow(PNMR,3)/pow(TKR,20);
    VMKR=PAR1+PAR2-PAR3-PAR4+PAR5;

    v=VMKR*3.17E-3;         /* specific volume */
    D=1/v;                  /* density */

    YD=-2*SA[0]*TKR+6*SA[1]/pow(TKR,7);
    SNUM=0.;

    for (i=1;i<=10;++i)
    {
      SNUM=SNUM+(i-2)*A[i]*pow(TKR,(i-1));
    }

    PRT1=A[11]*(Z*(17*(Z/29-Y/12)+5*TKR*YD/12)+SA[3]*TKR
                          -(SA[2]-1)*TKR*Y*YD)/pow(Z,(5/17));
    PRT2=PNMR*(A[12]-A[14])*TKR*TKR+A[15]*(9*TKR+SA[5])*pow((SA[5]-TKR),9)
                      +(A[16]*(20*pow(TKR,19)+SA[6]))/pow((SA[6]+pow(TKR,19)),2);
    PRT3=(12*pow(TKR,11)+SA[7])/pow((SA[7]+pow(TKR,11)),2) * (A[17]*PNMR
                                          +A[18]*PNMR*PNMR+A[19]*pow(PNMR,3));
    PRT4=A[20]*pow(TKR,18)*(17*SA[9]+19*TKR*TKR)*(1/pow((SA[9]+PNMR),3)
                                           +SA[10]*PNMR);
    PRT5=A[21]*SA[11]*pow(PNMR,3)+21*A[22]/pow(TKR,22) * pow(PNMR,4);
    ENTR=A[0]*TKR-SNUM+PRT1+PRT2-PRT3+PRT4+PRT5;
    H=ENTR*70120.4;         /* specific enthalpy */
    U=H-pw*v;               /* specific internal energy */

    return (H);
    }


    // Enthalpy of saturated vapor
    double hsat(double temperature) const
    {
    double h,x;
    x=(temperature-273.15)/100;          /* in Celsius/100 !! */

    h=2500514.22+188815.35*x-24027.4*pow(x,2)+29367.67*pow(x,3)-25873.63*pow(x,4)
                             +7886.38*pow(x,5)-982.6*pow(x,6);

//      if(isnan(h)){
//          sprintf(buf,"isnan hsat \n");
//      h = 0.0;
//      }

    return(h);
    }

    // Enthalpy of superheated steam
    double hsteam(double p_gw, double temperature) const
    {
    static double A[6];         /* Data A */
    static double C[7];     /* Data C */
    static double B[3][8];      /* Data B */
    static double SB[5];        /* Data SB */
    static double Z[3][8];      /* Data Z */
    double TK,PK,EL1,X,BEL,BELP,SUM,SC,d__1,S2,v,HK,D,H,U;
    double S[3],R[3];
    int j,k;
    double b,z;

    if(temperature < 270. || temperature > 500.)
    {
        DUNE_THROW(Dune::NumericalProblem,
            "ConstrelAir: Temperature " << temperature << " out of range at " << __FILE__ << ":" << __LINE__);
    }
    else if (p_gw < eps_ || p_gw > 1.0E8) {
        DUNE_THROW(Dune::NumericalProblem,
            "ConstrelAir: Pressure " << p_gw << " out of range at " << __FILE__ << ":" << __LINE__);
    }

    /* Initialization of data-arrays */
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

    /* Calculation */
    TK=temperature/647.3;
    PK=p_gw/2.212E7;
    EL1=4.260321148;
    X=exp(SB[0]*(1-TK));
    BEL=15.74373327-34.17061978*TK+19.31380707*pow(TK,2);
    BELP=-34.17061978+2*19.31380707*TK;
    SUM=EL1*(TK/PK);

    for (j=1;j<=5;++j)
    {
      SC=0.;
      for (k=1;k<=3;++k)
      {
        b=B[k-1][j-1];
        z=Z[k-1][j-1];
        SC=b*pow(X,z)+SC;
      }
      d__1=(double) (j-1);
      SUM=SUM-j*pow(PK,d__1)*SC;
    }

    S[0]=1/(pow(PK,4)+eps_) + SB[1]*pow(X,14);
    S[1]=1/(pow(PK,5)+eps_)+ SB[2]*pow(X,19);
    S[2]=1/(pow(PK,6)+eps_) + SB[3]*pow(X,54) + SB[4]*pow(X,27);

    for (j=6;j<=8;++j)
    {
      S2=0.;
      for (k=1;k<=3;++k)
      {
        S2=S2+B[k-1][j-1]*pow(X,(Z[k-1][j-1]));
      }
      d__1= (double) (1-j);
      SUM=SUM-(S2/pow(S[j-6],2))*(j-2)*pow(PK,d__1);
    }
    S2=0.;

    for (k=1;k<=7;++k)
    {
      d__1= (double) (k-1);
      S2=C[k-1]*pow(X,d__1)+S2;
    }

    SUM=SUM+11*pow((PK/BEL),10)*S2;
    v=SUM*0.00317;          /* specific volume */
    D=1/v;              /* density !! */

    HK=0.;
    HK=A[0]*TK;

    for (k=1;k<=5;++k)
    {
      d__1= (double) (k-1);
      HK=HK-A[k]*(k-2)*pow(TK,d__1);
    }

    for (j=1;j<=5;++j)
    {
      S2=0.;
      for (k=1;k<=3;++k)
      {
        S2=S2+B[k-1][j-1]*(1+Z[k-1][j-1]*SB[0]*TK)*pow(X,(Z[k-1][j-1]));
      }
      d__1= (double) j;
      HK=HK-pow(PK,d__1)*S2;
    }

    R[0]=14*SB[1]*pow(X,14);
    R[1]=19*SB[2]*pow(X,19);
    R[2]=54*SB[3]*pow(X,54)+27*SB[4]*pow(X,27);

    for (j=6;j<=8;++j)
    {
      S2=0.;
      for (k=1;k<=3;++k)
      {
        S2=S2+B[k-1][j-1]*pow(X,(Z[k-1][j-1]))*(1/S[j-6])*
                             (1+(Z[k-1][j-1]*SB[0]*TK-SB[0]*TK*R[j-6])/S[j-6]);
      }
      HK=HK-S2;
    }

    S2=0.;

    for (k=1;k<=7;++k)
    {
      d__1= (double) (k-1);
      S2=S2+(1+TK*(10*(BELP/BEL)+(k-1)*SB[0]))*C[k-1]*pow(X,d__1);
    }

    HK=HK+PK*pow((PK/BEL),10)*S2;

    H=HK*70120.4;           /* specific enthalpy */

    U=H-p_gw*v;         /* specific internal energy */


//      if(isnan(H)){
//          sprintf(buf,"isnan hsteam \n");
//      H = 0.0;
//      }

    return (H);
    }


    /* Interpolation polynoms for slightly overheated */
    /* steam, for 100/50/25/5/1 bar                    */
    double h_supst100(double temperature) const
    {
    double h,x;
    x=temperature/100;

    h=-28677461.4108+30630020.8681*x-11287756.0179*pow(x,2)
             +1865074.74575*pow(x,3)-115262.537121*pow(x,4);

    return(h);
    }


    double h_supst50(double temperature) const
    {
    double h,x;
    x=temperature/100;

    h=1188555.125+806071*x-37758.307*x*x-22191.355*x*x*x+3148.348*pow(x,4);


    return(h);
    }


    double h_supst25(double temperature) const
    {
    double h,x;
    x=temperature/100;

    h=239668.627+2761061.512*x-1170664.429*x*x+242749.898*pow(x,3)
              -18933.3687*pow(x,4);


    return(h);
    }

    double h_supst5(double temperature) const
    {
    double h,x;
    x=temperature/100;

    h=2287056.901+417531.975*x-106064.44*x*x+23491.494*pow(x,3)-1921.521*pow(x,4);

    return(h);
    }


    double h_supst1(double temperature) const
    {
    double h,x;
    x=temperature/100;

    h=2456564.043+241284.295*x-29300.9264*x*x+8110.1514*pow(x,3)-757.5625*pow(x,4);

    return(h);
    }



    /* Enthalpy of water vapor in the gasphase */
    double enthalpy_vapor(double Sw, double p_gw, double temperature) const
    {
    double T_sat=0.,epsilon,tsa,tsa_neu,g,g_strich,Celsius;
    int i;
    int crit_I,crit_II=0.;
    double a,h;
    double h0,DT_sat,Dh1,Dh2,Dh;
    double wicht;

    if(temperature < 250. || temperature > 500.)
    {
        DUNE_THROW(Dune::NumericalProblem,
            "ConstrelAir: Temperature " << temperature << " out of range at " << __FILE__ << ":" << __LINE__);
    }
    else if (p_gw < eps_ || p_gw > 1.0E8) {
        DUNE_THROW(Dune::NumericalProblem,
            "ConstrelAir: Pressure " << p_gw << " out of range at " << __FILE__ << ":" << __LINE__);
    }
    else if (Sw < 0.0 || Sw > 1.0) {
        DUNE_THROW(Dune::NumericalProblem,
            "ConstrelAir: Water saturation " << Sw << " out of range at " << __FILE__ << ":" << __LINE__);
    }

    Celsius=temperature-273.15;

        /* Nassdampf, leicht oder stark ueberhitzter Dampf ?? */

    if (Sw > eps_) crit_I=1;   /* Nassdampf ! */
    else
    {
        /* Ermittlung von T_sat(p_gw) */
      tsa=temperature;
      epsilon=1.E-2;
      for(i=1;i<100;++i)
      {
        g=pwsat(tsa)-p_gw;
        g_strich=(pwsat(tsa+epsilon)-pwsat(tsa))/epsilon;
        if (g_strich <= epsilon) g_strich = epsilon;
        tsa_neu=tsa-g/g_strich;
        a=tsa-tsa_neu;
        if (a<0) a=-a;
        if (a<epsilon)
        {
          T_sat=tsa;
          break;
        }
        tsa=tsa_neu;
      }


      if (Celsius/(T_sat-273.15)<1.8) crit_I=2;   /* leicht ueberhitzt */
      else crit_I=3;   /* stark ueberhitzt */
    }

    switch(crit_I)
    {
      case 1:
      {
        h=hsat(temperature);

//      if(isnan(h)){
//          sprintf(buf,"isnan enth_gw \n");
//      h = 0.0;
//      }

        return (h);
      }
      case 2:
      {
        /* Abfrage des Interpolationsbereichs */
        if (p_gw<=1.E5) crit_II=1;
        if ((p_gw>1.E5)&&(p_gw<=5.E5)) crit_II=2;
        if ((p_gw>5.E5)&&(p_gw<=25.E5)) crit_II=3;
        if ((p_gw>25.E5)&&(p_gw<=50.E5)) crit_II=4;
        if ((p_gw>50.E5)&&(p_gw<=100.E5)) crit_II=5;

        switch(crit_II)
        {
          case 1:
          {
        h0=hsat(T_sat);
        h=h0-1.759*(T_sat-temperature);
        wicht=(1.E5-p_gw)/1.E5;
        DT_sat=temperature-T_sat;
        Dh1=h_supst1(DT_sat+99.632)-h_supst1(99.632);
        h=h*wicht+(1-wicht)*(h0+Dh1);

        return (h);
          }
          case 2:
          {
        h0=hsat(T_sat);
        DT_sat=temperature-T_sat;
        Dh1=h_supst1(DT_sat+99.632)-h_supst1(99.632);
        Dh2=h_supst5(DT_sat+151.866)-h_supst5(151.866);
        Dh=(T_sat-273.15-99.632)*Dh2/52.231
            + (151.866-(T_sat-273.15))*Dh1/52.231;
        h=h0+Dh;

        return (h);
          }
          case 3:
          {
            h0=hsat(T_sat);
            DT_sat=temperature-T_sat;
            Dh1=h_supst5(DT_sat+151.866)-h_supst5(151.866);
            Dh2=h_supst25(DT_sat+223.989)-h_supst25(223.989);
            Dh=(T_sat-273.15-151.866)*Dh2/72.123
            + (223.989-(T_sat-273.15))*Dh1/72.123;
            h=h0+Dh;

        return (h);
          }
          case 4:
          {
            h0=hsat(T_sat);
            DT_sat=temperature-T_sat;
            Dh1=h_supst25(DT_sat+223.989)-h_supst25(223.989);
            Dh2=h_supst50(DT_sat+263.977)-h_supst50(263.977);
            Dh=(T_sat-273.15-223.989)*Dh2/39.988
            + (263.977-(T_sat-273.15))*Dh1/39.988;
            h=h0+Dh;

            return (h);
          }
          case 5:
          {
            h0=hsat(T_sat);
            DT_sat=temperature-T_sat;
            Dh1=h_supst50(DT_sat+263.977)-h_supst50(263.977);
            Dh2=h_supst100(DT_sat+311.031)-h_supst100(311.031);
            Dh=(T_sat-273.15-263.977)*Dh2/47.054
            + (311.031-(T_sat-273.15))*Dh1/47.054;
            h=h0+Dh;

            return (h);
          }
        }
      }
      case 3:
      {
        h=hsteam(p_gw,temperature);

//      if(isnan(h)){
//          sprintf(buf,"isnan enth_gw \n");
//      h = 0.0;
//      }

      return (h);
      }
    }
      return(1);
    }

    /* spezifische Enthalpie der Gasphase EINHEITEN ! */
//  double sp_enth2p2cni_g (double temperature, double Xag, double Xwg, double pg, double Sw)
    double sp_enth2p2cni_g (double temperature, double pg, double Xwg) const
    {
       double omega_gw, omega_ag; // mass fractions (notwendig?)
       double h_g, h_ga, h_gw, C_va, p_gw; //, m_w;
       const double Rluft = 287.2;   // gas const. for air [J/(kg K)]
       const double molecular_weight_water = 0.018016;
       const double molecular_weight_air = 0.02896;
       double Xag = 1.0 - Xwg;

       if(temperature < 250. || temperature > 500.)
       {
           DUNE_THROW(Dune::NumericalProblem,
               "ConstrelAir: Temperature " << temperature << " out of range at " << __FILE__ << ":" << __LINE__);
       }
       else if (pg < eps_ || pg > 1.0E8) {
           DUNE_THROW(Dune::NumericalProblem,
               "ConstrelAir: Pressure " << pg << " out of range at " << __FILE__ << ":" << __LINE__);
       }
//       else if (Sw < 0.0 || Sw > 1.0) {
//           DUNE_THROW(Dune::NumericalProblem,
//               "ConstrelAir: Water saturation " << Sw << " out of range at " << __FILE__ << ":" << __LINE__);
//       }
       else if (Xwg < 0.0 || Xwg > 1.0) {
           DUNE_THROW(Dune::NumericalProblem,
               "ConstrelAir: Mass fraction Xwg " << Xwg << " out of range at " << __FILE__ << ":" << __LINE__);
       }
       else if (Xag < 0.0 || Xag > 1.0) {
           DUNE_THROW(Dune::NumericalProblem,
               "ConstrelAir: Massfraction Xag " << Xag << " out of range at " << __FILE__ << ":" << __LINE__);
       }

//     /*Regularisierung*/
//     if(temperature<eps_)
//      {
//       UserWriteF("sp_enth2p2cni_g: \t Regularisation temperature \t %18.8g -> %g K\n",temperature,eps_);
//      temperature=eps_;
//      }
//     if(temperature>500.)
//      {
//       UserWriteF("sp_enth2p2cni_g: \t Regularisation temperature \t %18.8g -> 500 K\n",temperature);
//      temperature=500.;
//      }
//     if(pg<eps_)
//      {
//      UserWriteF("sp_enth2p2cni_g: \t Regularisation pg \t %18.8g -> %g Pa\n",pg,eps_);
//      pg=eps_;
//      }
//     if(pg>1.E8)
//      {
//       UserWriteF("sp_enth2p2cni_g: \t Regularisation pg \t %18.8g -> 1.E8 Pa\n",pg);
//      pg=1.E8;
//      }
//     if(Sw<0.0)
//      {
//       UserWriteF("sp_enth2p2cni_g: \t Regularisation Sw \t %18.8g -> 0.0 \n",Sw);
//      Sw = 0.0;
//      }
//     if(Sw>1.)
//      {
//       UserWriteF("sp_enth2p2cni_g: \t Regularisation Sw \t %18.8g -> 1.0 \n",Sw);
//      Sw = 1.;
//      }
//     if(Xwg<0.0)
//      {
//       UserWriteF("sp_enth2p2cni_g: \t Regularisation Xwg \t %18.8g -> 0.0 \n",Xwg);
//      Xwg = 0.0;
//      }
//     if(Xwg>1.)
//      {
//       UserWriteF("sp_enth2p2cni_g: \t Regularisation Xwg \t %18.8g -> 1.0 \n",Xwg);
//      Xwg = 1.;
//      }
//     if(Xag<0.0)
//      {
//       UserWriteF("sp_enth2p2cni_g: \t Regularisation Xag \t %18.8g -> 0.0 \n",Xag);
//      Xag = 0.0;
//      }
//     if(Xag>1.)
//      {
//       UserWriteF("sp_enth2p2cni_g: \t Regularisation Xag \t %18.8g -> 1.0 \n",Xag);
//      Xag = 1.;
//      }
       omega_ag = (Xag*molecular_weight_air)/(Xag*molecular_weight_air
                   + Xwg*molecular_weight_water);
       omega_gw = (Xwg*molecular_weight_air)/(Xag*molecular_weight_air
                   + Xwg*molecular_weight_water);

       //HACK:
       double Sw = 1.0;

       p_gw = partial_pressure_gas_w(pg,temperature,Sw,Xwg);
       h_gw = enthalpy_vapor(Sw,p_gw,temperature); // enth. of steam
       C_va = 733.;   /* Waermekapazitaet von Luft [J/(kg * K)] */
       h_ga = C_va*(temperature-273.15) + Rluft*temperature; // enth. of gas
       h_g = omega_ag*h_ga+omega_gw*h_gw; //

//      /**   **/
//      if(isnan(h_g)){
//          sprintf(buf,"isnan sp_enthalpy_g \n");
//      h_g = 0.0;
//      }
       return(h_g);
    }
};

} // end namespace Dune

#endif
