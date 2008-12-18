
class Water
{
public:
  double viscosity (double T, double p)
  {
    return 1.; //[kg/(ms)]
  }
  double density (double T, double p)
  {
    return 1000.; // [kg/m^3]
  }
  double Sr()
  {
    return 0.0;
  }
  double henry (double T)
  {
    return (0.8942 + 1.47 * exp(-0.04394*T) )*1e-10; // [1/Pa]
  }
  double vaporPressure (double T)
  {
    return 1228.; // [Pa]
  }
  static const double molarMass = 0.018016; // [kg/mole]
};

class Air
{
public:

  double gasconstant ()
  {
	  return 286.991; /* individual gas const. for air [J/(kg K)] */
  }
  double viscosity ( double T, double p)
  {
   // return 17.75e-6;//[kg/(ms)] at 15°C
	  return 1.e-6;
  }
  double density ( double T, double p)
  {
    return 1.23; // [kg/m^3] at 15°C
  }
  double Sr()
  {
    return 0.1;
  }
  static const double molarMass = 0.02896; // [kg/mole]
};

class Uniform
{
public:
  double viscosity ( double T, double p)
  {
    return 1.0;//[kg/(ms)]
  }
  double density ( double T, double p)
  {
    return 1.0; // [kg/m^3]
  }
  double Sr()
  {
    return 0.0;
  }
  static const double molarMass = 1.0; // [kg/mole]
};


