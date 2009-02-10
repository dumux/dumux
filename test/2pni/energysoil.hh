//$Id: energysoil.hh 882 2008-12-04 09:05:55Z melanie $
#ifndef DUNE_TWOPHEATSOIL_HH
#define DUNE_TWOPHEATSOIL_HH

/**
 * @file
 * @brief  Class for defining an instance of a Matrix2p soil
 * @author Bernd Flemisch, Klaus Mosthaf
 */

namespace Dune
{

  template<class Grid, class Scalar>
  class TwoPHeatSoil: public Matrix2p<Grid,Scalar>
  {
  public:
      typedef typename Grid::Traits::template Codim<0>::Entity Element;
      typedef typename Grid::ctype CoordScalar;

      typedef BasicDomain<Grid, Scalar> ParentType;
      // the domain traits of the domain
      typedef typename ParentType::DomainTraits DomainTraits;
      typedef typename DomainTraits::LocalPosition LocalPosition;
      typedef typename DomainTraits::GlobalPosition GlobalPosition;

      enum {dim=Grid::dimension, numEq=1};

      // define PERMEABILITY tensor
      virtual const FieldMatrix<CoordScalar,dim,dim>& K (const GlobalPosition& globalPos, const Element& element, const LocalPosition& localPos)
      {
            return K_;
      }
      virtual double porosity(const GlobalPosition& globalPos, const Element& element, const LocalPosition& localPos) const
      {
          return 0.15;
      }

      virtual double Sr_w(const GlobalPosition& globalPos, const Element& element, const LocalPosition& localPos, const double T) const
      {
              return 0.2;
      }

      virtual double Sr_n(const GlobalPosition& globalPos, const Element& element, const LocalPosition& localPos, const double T) const
      {
          return 0.05;
      }

      virtual typename Matrix2p<Grid,Scalar>::modelFlag relPermFlag(const GlobalPosition& globalPos, const Element& element, const LocalPosition& localPos) const
      {
          return Matrix2p<Grid,Scalar>::brooks_corey;
      }

    virtual double heatCap(const GlobalPosition& globalPos, const Element& element, const LocalPosition& localPos) const
    {
        return     (800 /* spec. heat cap. of sediment */
                        * 2650 /* density of sediment */
                        * (1-porosity(globalPos, element, localPos)));
    }

    virtual double heatCond(const GlobalPosition& globalPos, const Element& element, const LocalPosition& localPos, const double sat) const
    {
        static const double ldry = 0.32;
        static const double lsat = 2.7;
        double l_wurz, Sw;
        Sw = sat;
        if(Sw<0.0) Sw = 0.0; /* ACHTUNG Regularisierung */
        if(Sw>1.) Sw = 1.; /* ACHTUNG Regularisierung */

        l_wurz = ldry + sqrt(Sw)*(lsat-ldry);

        if(isnan(l_wurz)) {
            std::cout <<"isnan heatcondwurzel \n"<<std::endl;
            l_wurz = 0.0;
        }
        return(l_wurz);
    }

      virtual std::vector<double> paramRelPerm(const GlobalPosition& globalPos, const Element& element, const LocalPosition& localPos, const double T) const
      {
          // example for Brooks-Corey parameters
        std::vector<double> param(5);

            param[0] = 2.;
            param[1] = 10000;

          return param;
      }

      TwoPHeatSoil()
      : Matrix2p<Grid,Scalar>()
      {
          K_ = 0;
          for(int i = 0; i < dim; i++)
          {
              K_[i][i] = 2e-14;
          }
      }

      ~TwoPHeatSoil()
      {}

  private:
      FieldMatrix<CoordScalar,dim,dim> K_;
   };

} // end namespace
#endif

