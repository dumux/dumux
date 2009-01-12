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

  template<class G, class RT>
  class TwoPHeatSoil: public Matrix2p<G,RT>
  {
  public:
      typedef typename G::Traits::template Codim<0>::Entity Entity;
      typedef typename G::ctype DT;
      enum {dim=G::dimension, m=1};

      // define PERMEABILITY tensor
      virtual FieldMatrix<DT,dim,dim> K (const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi)
      {
            return K_;
      }
      virtual double porosity(const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi) const
      {
          return 0.15;
      }

      virtual double Sr_w(const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi, const double T) const
      {
              return 0.2;
      }

      virtual double Sr_n(const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi, const double T) const
      {
          return 0.05;
      }

      virtual typename Matrix2p<G,RT>::modelFlag relPermFlag(const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi) const
      {
          return Matrix2p<G,RT>::brooks_corey;
      }

    virtual double heatCap(const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi) const
    {
        return     (800 /* spec. heat cap. of sediment */
                        * 2650 /* density of sediment */
                        * (1-porosity(x, e, xi)));
    }

    virtual double heatCond(const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi, const double sat) const
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

      virtual std::vector<double> paramRelPerm(const FieldVector<DT,dim>& x, const Entity& e, const FieldVector<DT,dim>& xi, const double T) const
      {
          // example for Brooks-Corey parameters
        std::vector<double> param(5);

            param[0] = 2.;
            param[1] = 10000;

          return param;
      }

      TwoPHeatSoil()
      : Matrix2p<G,RT>()
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
      FieldMatrix<DT,dim,dim> K_;
   };

} // end namespace
#endif

