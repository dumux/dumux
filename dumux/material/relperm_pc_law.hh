// $Id:$

/*****************************************************************************
* Copyright (C) 2008 by Jochen Fritz                                         *
* Institute of Hydraulic Engineering                                         *
* University of Stuttgart, Germany                                           *
* email: <givenname>.<name>@iws.uni-stuttgart.de                             *
*                                                                            *
* This program is free software; you can redistribute it and/or modify       *
* it under the terms of the GNU General Public License as published by       *
* the Free Software Foundation; either version 2 of the License, or          *
* (at your option) any later version, as long as this copyright notice       *
* is included in its original form.                                          *
*                                                                            *
* This program is distributed WITHOUT ANY WARRANTY.                          *
*****************************************************************************/

#ifndef DUNE_RELPERM_PC_LAW_HH
#define DUNE_RELPERM_PC_LAW_HH

#include <dumux/material/property_baseclasses.hh>

/**
 * \author Jochen Fritz
 */


namespace Dune
{
/*!\ingroup 2pRel
 * \brief Base class for relative permeability and capillary pressure - saturation relationships.
 * Derived classes of this base class are used as input for the TwoPhaseRelations class.
 * The functions get significant input from the member soil which is an object of Matrix2p.
 * Specifications of model parameters for the relative permeability and capillary pressure
 * functions have to be made in right order. Read further details in the derived classes!
 */
template<class Grid>
class RelPerm_pc {
public:
    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef typename Grid::ctype Scalar;
    enum {dim=Grid::dimension};
    /*! \brief the capillary pressure - saturation relation
     *
     *  \param saturationW the saturation of the wetting phase
     *  \param globalPos position in global coordinates
     *  \param element codim 0 entity for which the value is sought
     *  \param localPos position in local coordinates in element
     *  \param param standard vector containing the parameters for the material law
     *  \return the capillary pressure \f$ p_\text{c} (S_\text{w})\f$.
     */
    virtual double pC (double saturationW, const FieldVector<Scalar,dim>& globalPos, const Element& element, const FieldVector<Scalar,dim>& localPos,
                       const std::vector<double>& param, const double temperature=283.15) const = 0;

    /*! \brief the derivative of capillary pressure w.r.t. the saturation
     *
     *  \param saturationW the saturation of the wetting phase
     *  \param param standard vector containing the parameters for the material law
     *  \param temperature temperature
     *  \param globalPos position in global coordinates
     *  \param element codim 0 entity for which the value is sought
     *  \param localPos position in local coordinates in element
     *  \return the derivative \f$\text{d}p_\text{c}/\text{d}S_\text{element}\f$
     */
    virtual double dPdS (double saturationW, const FieldVector<Scalar,dim>& globalPos, const Element& element, const FieldVector<Scalar,dim>& localPos,
                         const std::vector<double>& param, const double temperature=283.15) const = 0;

    /*! \brief the wetting phase saturation w.r.t. the capillary pressure
     *
     *  \param pC the capillary pressure
     *  \param temperature temperature
     *  \param globalPos position in global coordinates
     *  \param element codim 0 entity for which the value is sought
     *  \param localPos position in local coordinates in element
     *  \return the wetting phase saturation
     */
    virtual double saturationW (double pC, const FieldVector<Scalar,dim>& globalPos, const Element& element, const FieldVector<Scalar,dim>& localPos, const double temperature=283.15) const = 0;

    /*! \brief the derivative of the saturation w.r.t. the capillary pressure
     *
     *  \param pC the capillary pressure
     *  \param temperature temperature
     *  \param globalPos position in global coordinates
     *  \param element codim 0 entity for which the value is sought
     *  \param localPos position in local coordinates in element
     *  \return the derivative \f$\text{d}S_w/\text{d}p_\text{c}\f$
     */
    virtual double dSdP (double pC, const FieldVector<Scalar,dim>& globalPos, const Element& element, const FieldVector<Scalar,dim>& localPos, const double temperature=283.15) const = 0;

    const bool isLinear() const
    {
        return linear_;
    }

    /*! \brief wetting phase relative permeability saturation relationship
     *
     *  \param saturationW the saturation of the wetting phase
     *  \param globalPos position in global coordinates
     *  \param element codim 0 entity for which the value is sought
     *  \param localPos position in local coordinates in element
     *  \return the wetting phase relative permeability
     */
    virtual double krw (const double saturationW, const FieldVector<Scalar,dim>& globalPos, const Element& element, const FieldVector<Scalar,dim>& localPos, double temperature=283.15) const = 0;

    /*! \brief nonwetting phase relative permeability saturation relationship
     *
     *  \param saturationN the saturation of the nonwetting phase
     *  \param globalPos position in global coordinates
     *  \param element codim 0 entity for which the value is sought
     *  \param localPos position in local coordinates in element
     *  \return the nonwetting phase relative permeability
     */
    virtual double krn (const double saturationN, const FieldVector<Scalar,dim>& globalPos, const Element& element, const FieldVector<Scalar,dim>& localPos, double temperature=283.15) const = 0;


    /** \brief relative permeability saturation relationship for both phases
     *  In many cases the relative permeabilities of both phases are needed at
     *  the same time. This function reduces unnecessary computational costs.
     *  \param saturationN the saturation of the nonwetting phase
     *  \param globalPos position in global coordinates
     *  \param element codim 0 entity for which the value is sought
     *  \param localPos position in local coordinates in element
     *  \return relative permeability vector: first entry wetting, seconde entry nonwetting phase
     */
    virtual std::vector<double> kr (const double saturationW, const FieldVector<Scalar,dim>& globalPos, const Element& element, const FieldVector<Scalar,dim>& localPos, double temperature=283.15) const = 0;

    /** \brief constructor
     *  \param soil a matrix property object.
     *  \param wP phase property object for the wetting Phase.
     *  \param nP phase property object for the nonwetting Phase.
     *  \param lin true specifies a linear model. Usually false. Only set true if you know what you are doing!
     */
    RelPerm_pc(const Matrix2p<Grid,double>& soil, const bool lin = false)
        : soil_(soil), linear_(lin)
    {
    }

    virtual ~RelPerm_pc()
    {
    }

protected:
    const Matrix2p<Grid,double>& soil_;
    const bool linear_;
};


/** \ingroup 2pRel
 * @brief Represents the linear relative permability - saturation relation \f$ k_{r\alpha}(S) = \frac{S-S_{r\alpha}}{1-S_{r\alpha}} \f$
 * Vector entries in Matrix2p::paramRelPerm must be in the order
 *         - minimum capillary pressure
 *         - maximum capillary pressure
 */
template<class Grid>
class LinearLaw : public RelPerm_pc<Grid>
{
public:
    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef typename Grid::ctype Scalar;
    enum {dim=Grid::dimension};

    double krw (double saturationW, const FieldVector<Scalar,dim>& globalPos, const Element& element, const FieldVector<Scalar,dim>& localPos, double temperature=283.15) const
    {
        double Sr_w = soil_.Sr_w(globalPos, element, localPos, temperature);
        double Sr_n = soil_.Sr_n(globalPos, element, localPos, temperature);
        return std::max(std::min((saturationW - Sr_w)/(1- Sr_w - Sr_n), 1.), 0.);
    }

    double krn (double saturationN, const FieldVector<Scalar,dim>& globalPos, const Element& element, const FieldVector<Scalar,dim>& localPos, double temperature=283.15) const
    {
        double Sr_w = soil_.Sr_w(globalPos, element, localPos, temperature);
        double Sr_n = soil_.Sr_n(globalPos, element, localPos, temperature);
        return std::max(std::min((saturationN - Sr_n)/(1- Sr_w - Sr_n), 1.), 0.);
    }

    std::vector<double> kr (const double saturationW, const FieldVector<Scalar,dim>& globalPos, const Element& element, const FieldVector<Scalar,dim>& localPos, double temperature=283.15) const
    {
        std::vector<double> kr(2);
        double Sr_w = soil_.Sr_w(globalPos, element, localPos, temperature);
        double Sr_n = soil_.Sr_n(globalPos, element, localPos, temperature);
        kr[0] = std::max(std::min((saturationW - Sr_w)/(1- Sr_w - Sr_n), 1.), 0.);
        kr[1] = std::max(std::min((1 - saturationW - Sr_n)/(1- Sr_w - Sr_n), 1.), 0.);
        return kr;
    }

    double pC (double saturationW, const FieldVector<Scalar,dim>& globalPos, const Element& element, const FieldVector<Scalar,dim>& localPos,
               const std::vector<double>& param, double temperature=283.15) const
    {
        double Swr = soil_.Sr_w(globalPos, element, localPos, temperature);
        double Snr = soil_.Sr_n(globalPos, element, localPos, temperature);

        if (saturationW > (1-Snr)) return param[0]; // min pc
        if (saturationW < Swr) return param[1]; // max pc

        return  param[0] + (param[1] - param[0]) * (1 - saturationW - Swr) / (1-Swr-Snr);
    }

    double dPdS (double saturationW, const FieldVector<Scalar,dim>& globalPos, const Element& element, const FieldVector<Scalar,dim>& localPos,
                 const std::vector<double>& param, double temperature=283.15) const
    {
        double Swr = soil_.Sr_w(globalPos, element, localPos, temperature);
        double Snr = soil_.Sr_n(globalPos, element, localPos, temperature);

        return (param[1] - param[0]) * (-1)/(1-Swr-Snr);
    }

    double saturationW (double pC, const FieldVector<Scalar,dim>& globalPos, const Element& element, const FieldVector<Scalar,dim>& localPos, double temperature=283.15) const
    {
        double Swr = soil_.Sr_w(globalPos, element, localPos, temperature);
        double Snr = soil_.Sr_n(globalPos, element, localPos, temperature);

        std::vector<double> param = soil_.paramRelPerm(globalPos, element, localPos, temperature);

        double Sw = 1-Snr - (pC-param[0])/(param[1]-param[0])*(1-Swr-Snr);
        if (Sw > (1-Snr)) return (1-Snr);
        if (Sw < Swr) return (Swr);

        return ( Sw);
    }

    double dSdP (double pC, const FieldVector<Scalar,dim>& globalPos, const Element& element, const FieldVector<Scalar,dim>& localPos, double temperature=283.15) const
    {
        double Swr = soil_.Sr_w(globalPos, element, localPos, temperature);
        double Snr = soil_.Sr_n(globalPos, element, localPos, temperature);

        std::vector<double> param = soil_.paramRelPerm(globalPos, element, localPos, temperature);

        return ( (1-Swr-Snr) / (param[0] - param[1])  ) ;
    }

    LinearLaw(const Matrix2p<Grid,double>& soil, bool lin = false)
        : RelPerm_pc<Grid>(soil, false), soil_(soil)
    {     }

private:
    double maxpc;
    double minpc;
    const Matrix2p<Grid,double>& soil_;
};


/*!\ingroup 2pRel
 * \brief van Genuchten mobility/saturation relation.
 *
 *  Employs the van Genuchten non-linear relative permeability/saturation relation.
 *  Vector entries in Matrix2p::paramRelPerm must be in the order
 *         - m
 *         - n
 *         - \f$ \epsilon \f$
 *         - \f$ \gamma \f$
 *         - \f$ \alpha \f$
 *
 */
template<class Grid>
class VanGenuchtenLaw : public RelPerm_pc<Grid>
{
public:

    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef double Scalar;
    enum {dim=Grid::dimension};

    double krw (double saturationW, const FieldVector<Scalar,dim>& globalPos, const Element& element, const FieldVector<Scalar,dim>& localPos, double temperature=283.15) const
    {
        double Swr = soil_.Sr_w(globalPos, element, localPos, temperature);
        double Snr = soil_.Sr_n(globalPos, element, localPos, temperature);

        double Se,krw,r;

        /* effective values */
        Se = (saturationW - Swr) / (1 - Snr - Swr);

        /* regularization */
        if(Se > 1.) return 1.;
        if(Se < machineEps_) Se = machineEps_;

        std::vector<double> param = soil_.paramRelPerm(globalPos, element, localPos, temperature);
        double m = param[0];
        double eps = param[2];

        /* compute value */
        r   = 1 - pow(Se, 1/m);
        krw = pow(Se, eps) * pow(1 - pow(r,m), 2);
        return(krw);
    }

    double krn (double saturationN, const FieldVector<Scalar,dim>& globalPos, const Element& element, const FieldVector<Scalar,dim>& localPos, double temperature=283.15) const
    {
        double Swr = soil_.Sr_w(globalPos, element, localPos, temperature);
        double Snr = soil_.Sr_n(globalPos, element, localPos, temperature);

        double Se, r;

        /* effective values */
        Se = (1 - saturationN - Swr) / (1 - Snr - Swr);

        /* effective Saturation Se has to be between 0 and 1! */
        if(Se > 1.) Se = 1.;
        if(Se < machineEps_) Se = machineEps_;

        std::vector<double> param = soil_.paramRelPerm(globalPos, element, localPos, temperature);
        double m = param[0];
        double gamma = param[3];

        /* compute value */
        r   = 1 - pow(Se, 1/m);
        return pow(1-Se, gamma) * pow(r, 2*m);
    }

    std::vector<double> kr (const double saturationW, const FieldVector<Scalar,dim>& globalPos, const Element& element, const FieldVector<Scalar,dim>& localPos, double temperature=283.15) const
    {
        std::vector<double> kr(2);
        // residual saturations
        double Srw = soil_.Sr_w(globalPos, element, localPos, temperature);
        double Srn = soil_.Sr_n(globalPos, element, localPos, temperature);
        // effective saturation
        double Se = (saturationW - Srw) / (1 - Srw - Srn);

        // regularization
        if(Se > 1.)
        {
            kr[0] = 1;
            kr[1] = 0;
        }
        if(Se < machineEps_) Se = machineEps_;

        // get Van Genuchten parameters
        std::vector<double> param = soil_.paramRelPerm(globalPos, element, localPos, temperature);
        double m = param[0];
        double eps = param[2];
        double gamma = param[3];

        // compute values
        double r   = 1 - pow(Se, 1/m);
        kr[0] = pow(Se, eps) * pow(1 - pow(r,m), 2);
        kr[1] = pow(1-Se, gamma) * pow(r, 2*m);
        return kr;
    }

    double pC (const double saturationW, const FieldVector<Scalar,dim>& globalPos, const Element& element, const FieldVector<Scalar,dim>& localPos,
               const std::vector<double>& param, double temperature=283.15) const
    {
        double r, x_, vgM;
        double pc, pc_prime, Se_regu;
        int asymptotic;

        double Swr = soil_.Sr_w(globalPos, element, localPos, temperature);
        double Snr = soil_.Sr_n(globalPos, element, localPos, temperature);

        double m = param[0];
        double n = param[1];
        if (m - (1-1/n) > 1e-5 || m - (1-1/n) < -1e-5)
        	std::cerr << "m and n do not fit!! \n";

        double alpha = param[4];

        double Se = (saturationW - Swr) / (1. - Swr - Snr);

        if(Se < 0.0) Se = 0.0;
        if(Se > 1.0) Se = 1.0;

        /* check size of S_e^(-1/m) */
        if (Se < epsPC_)
            asymptotic = 1;
        else if (Se > 1 - epsPC_)
            asymptotic = 0;
        else if (1.0 / pow(Se, 1.0/m) > 1000.0)
            asymptotic = 1;
        else asymptotic = 0;

        if (asymptotic) /* use pc_VanG = 1/alpha/pow(Se,1/(n-1)) */
        {
            if (Se > epsPC_)
                return(1.0 / (alpha * pow(Se, 1.0 / (n-1))));
            else /* regularize with tangent */
            {
//                pc  = 1.0 / (alpha * pow(epsPC_, 1.0 / (n-1)));
                pc  = pow(pow(epsPC_, -1.0/m) - 1.0, 1.0/n) / alpha;
                pc_prime = 1.0 / (pow(epsPC_, n/(n-1)) * alpha * (1-n) * (1-Swr-Snr));
                return((Se - epsPC_) * pc_prime + pc);
            }
        }
        else
        {    /* use correct van Genuchten curve */
            if (Se > epsPC_ && Se < 1 - epsPC_)
            {
                r = pow(Se, -1/m);
                x_ = r - 1;
                vgM = 1 - m;
                x_ = pow(x_, vgM);
                r = x_ / alpha;
                return(r);
            }
            else
            {
                /* value and derivative at regularization point */
                if (Se <= epsPC_) Se_regu = epsPC_;
                else Se_regu = 1 - epsPC_;
                pc       = pow(pow(Se_regu, -1/m) - 1, 1/n) / alpha;
                pc_prime = pow(pow(Se_regu, -1/m) - 1, 1/n-1) * pow(Se_regu, -1/m-1) * (-1/m) / alpha / (1-Snr-Swr);

                /* evaluate tangential */
                r        = pc;//(Se - Se_regu) * pc_prime + pc;
                if (r<0) r=0; // if Sw=1, pc with correct van Genuchten curve formulation is negative
                return(r);
            }
        }
    }

    double dPdS (double saturationW, const FieldVector<Scalar,dim>& globalPos, const Element& element, const FieldVector<Scalar,dim>& localPos,
                 const std::vector<double>& param, double temperature=283.15) const
    {
        double r, x_;

        double Swr = soil_.Sr_w(globalPos, element, localPos, temperature);
        double Snr = soil_.Sr_n(globalPos, element, localPos, temperature);

        double m = param[0];
        //            double n = param[1];
        double alpha = param[4];

        double Se = (saturationW - Swr) / (1. -Swr - Snr);

        /* regularization */
        if (Se <= 0.0) Se = 1.E-3;
        if (Se >= 1.0) Se = 1.0 - 1.E-5;

        /* compute value */
        r = pow(Se, -1/m);
        x_ = pow((r-1), (1-m));
        r = -(1-0.0) / alpha * x_ * (1-m) / m * r / (r-1) / Se / (1-Snr-Swr);
        return(r);
    }

    double saturationW (double pC, const FieldVector<Scalar,dim>& globalPos, const Element& element, const FieldVector<Scalar,dim>& localPos, double temperature=283.15) const
    {

        double QT1 = 0.001;
        double QT2 = 0.979;
        double Swr = soil_.Sr_w(globalPos, element, localPos, temperature);
        double Snr = soil_.Sr_n(globalPos, element, localPos, temperature);

        std::vector<double> param = soil_.paramRelPerm(globalPos, element, localPos, temperature);

        double m = param[0];
        double n = param[1];
        double alpha = param[4];

        /* check left tangent */
        double r = pow(QT1,-1/m);
        double x_ = pow(r-1,1/n);
        double pc = (x_/alpha);
        if (pC>=pc)
        {
            double pc_prime = (1-m)/alpha/pow(pow(QT1,-1/m)-1,m)/pow(QT1,1+1/m)*(-1/m);
            double Se = (pC-pc)/pc_prime+QT1;
            return (Se*(1.0-Swr-Snr)+Swr);
        }

        /* check right part */

        r = pow(QT2,-1/m);
        x_ = pow(r-1,1-m);
        pc = (x_/alpha);
        if (pC<=pc)
        {
            double Se = (1.0-pC/pc)*(1-QT2)+QT2;
            return (Se*(1.0-Swr-Snr)+Swr);
        }

        /* do inversion ... */

        double Se = 1.0/pow(pow(alpha*pC,n)+1.0,m);
        return (Se*(1.0-Swr-Snr)+Swr);

    }

    double dSdP (double pC, const FieldVector<Scalar,dim>& globalPos, const Element& element, const FieldVector<Scalar,dim>& localPos, double temperature=283.15) const
    {
        double QT2 = 0.979;
        double Swr = soil_.Sr_w(globalPos, element, localPos, temperature);
        double Snr = soil_.Sr_n(globalPos, element, localPos, temperature);

        std::vector<double> param = soil_.paramRelPerm(globalPos, element, localPos, temperature);

        double n = param[1];
        double alpha = param[4];

        double dsdSe = 1-Swr-Snr;

        double pcReg = this->pC (QT2, globalPos, element, localPos, param, temperature);
        double dswdpc_reg =-(n-1)*pow(alpha*pcReg,n)*pow(1+pow(alpha*pcReg,n),-2+1/n)/pcReg*dsdSe;
        double dswdpc = -(n-1)*pow(alpha*pC,n)*pow(1+pow(alpha*pC,n),-2+1/n)/pC*dsdSe;

        if (pC<pcReg)
        {
            dswdpc = (pC-pcReg)/pcReg*dswdpc_reg + dswdpc_reg;
        }

        return dswdpc;
    }

    VanGenuchtenLaw(const Matrix2p<Grid,double>& soil, bool lin = false)
        : RelPerm_pc<Grid>(soil, false), soil_(soil)
    {
    }

protected:
    static const double epsPC_ = 5e-6; // 1e-10; Alex:  //!< threshold for linearization of capillary pressure
    static const double machineEps_ = 1e-15;
    const Matrix2p<Grid,double>& soil_;
};


/*!\ingroup 2pRel
 * \brief Brooks-Corey mobility/saturation relation.
 *
 *  Employs the Brooks-Corey non-linear relative permeability/saturation relation, namely,
 *  \f{align*}
 *  S_\text{element} = \frac{S - S_\text{r}}{1 - S_\text{r}}, \quad
 *  \lambda_\text{w} = \mu_\text{w}^{-1} S_\text{element}^{\frac{2 + 3\lambda}{\lambda}},  \quad
 *  \lambda_\text{n} = \mu_\text{n}^{-1} \left( 1 - S_\text{element}^{\frac{2 + 2\lambda}{\lambda}}\right)
 *  \left( 1 - S_\text{element} \right)^2, \quad
 *  \lambda = \lambda_\text{w} + \lambda_\text{n}.
 *  \f}
 *
 *  Vector entries in Matrix2p::paramRelPerm must be in the order
 *         - \f$ lambda \f$
 *         - entry pressure \f$ p_c \f$
 */
template<class Grid>
class BrooksCoreyLaw : public RelPerm_pc<Grid>
{
public:
    typedef typename Grid::Traits::template Codim<0>::Entity Element;
    typedef typename Grid::ctype Scalar;
    enum {dim=Grid::dimension};

    double pC (double saturationW, const FieldVector<Scalar,dim>& globalPos, const Element& element, const FieldVector<Scalar,dim>& localPos,
               const std::vector<double>& param, const double temperature) const
    {
        //effective Saturation
        double Se = (saturationW - soil_.Sr_w(globalPos, element, localPos, temperature))
            /(1. - soil_.Sr_w(globalPos, element, localPos, temperature) - soil_.Sr_n(globalPos, element, localPos, temperature));

        double lambda = param[0];
        double p0 = param[1];

        if (Se > SweMin_)
            return brooksCorey_(Se, p0, lambda);
        else
        {
            // todo (0): better: use a splines between [SweMin/2, SweMin]
            double pC_0 = brooksCorey_(SweMin_/2, p0, lambda);
            double pC_1 = brooksCorey_(SweMin_, p0, lambda);

            double m = (pC_1 - pC_0)/(SweMin_/2);
            double b = pC_1 - m*(SweMin_);
            return m*Se + b;
        }
    }


    double dPdS (double saturationW, const FieldVector<Scalar,dim>& globalPos, const Element& element, const FieldVector<Scalar,dim>& localPos,
                 const std::vector<double>& param, const double temperature=283.15) const
    {
        double Swr = soil_.Sr_w(globalPos, element, localPos, temperature);
        double Snr = soil_.Sr_n(globalPos, element, localPos, temperature);

        double lambda = param[0];
        double p0 = param[1];

        double Se = (saturationW - Swr)
            /(1. - Swr - Snr);

        if(Se<0.0) Se = 0.0;
        if(Se>1.0) Se = 1.0;

        if (Se > SweMin_)
            return (-p0/lambda*pow(Se, -1.0/lambda-1)/(1-Snr-Swr));
        else
        {
            double dSedSwsquare=1/(1-Snr-Swr)/(1-Snr-Swr);
            return (-p0*dSedSwsquare/lambda/pow(SweMin_, (1+1/lambda)));
        }
    }

    double saturationW (double pC, const FieldVector<Scalar,dim>& globalPos, const Element& element, const FieldVector<Scalar,dim>& localPos, const double temperature=283.15) const
    {
        double lambda = soil_.paramRelPerm(globalPos, element, localPos, temperature)[0];
        double p0 = soil_.paramRelPerm(globalPos, element, localPos, temperature)[1];
        return ((1 - soil_.Sr_w(globalPos, element, localPos, temperature)) * pow(p0/pC, lambda) + soil_.Sr_n(globalPos, element, localPos, temperature));
    }

    double dSdP (double pC, const FieldVector<Scalar,dim>& globalPos, const Element& element, const FieldVector<Scalar,dim>& localPos, const double temperature=283.15) const
    {
        double lambda = soil_.paramRelPerm(globalPos, element, localPos, temperature)[0];
        double p0 = soil_.paramRelPerm(globalPos, element, localPos, temperature)[1];
        return (-(1 - soil_.Sr_w(globalPos, element, localPos, temperature))*pow(p0/pC, lambda-1)*lambda*pow(1.0/pC, 2));
    }


    double krw (double saturationW, const FieldVector<Scalar,dim>& globalPos, const Element& element, const FieldVector<Scalar,dim>& localPos, const double temperature=283.15) const
    {
        double Se = (saturationW - soil_.Sr_w(globalPos, element, localPos, temperature))
            /(1. - soil_.Sr_w(globalPos, element, localPos, temperature) - soil_.Sr_n(globalPos, element, localPos, temperature));

        //regularisation
        if (Se > 1) return 1.0;
        if (Se < epsKr_) return 0.0;

        double lambda = soil_.paramRelPerm(globalPos, element, localPos, temperature)[0];

        double exponent = (2. + 3*lambda) / lambda;
        return pow(Se, exponent);
    }


    double krn (double saturationN, const FieldVector<Scalar,dim>& globalPos, const Element& element, const FieldVector<Scalar,dim>& localPos, const double temperature=283.15) const
    {
        double Se = ((1.-saturationN) - soil_.Sr_w(globalPos, element, localPos, temperature))
            /(1. - soil_.Sr_w(globalPos, element, localPos, temperature) - soil_.Sr_n(globalPos, element, localPos, temperature));

        //regularisation
        if (Se > 1) return 0.0;
        if (Se < epsKr_) return 1.0;

        double lambda = soil_.paramRelPerm(globalPos, element, localPos, temperature)[0];
        double exponent = (2. + lambda) / lambda;
        return pow(1.-Se, 2) * ( 1. - pow(Se, exponent) );
    }

    std::vector<double> kr (const double saturationW, const FieldVector<Scalar,dim>& globalPos, const Element& element, const FieldVector<Scalar,dim>& localPos, const double temperature=283.15) const
    {
        std::vector<double> kr(2);
        double Srw = soil_.Sr_w(globalPos, element, localPos, temperature);
        double Srn = soil_.Sr_n(globalPos, element, localPos, temperature);
        double Se = (saturationW - Srw) / (1 - Srw - Srn);
        // regularization
        if (Se > 1)
        {
            kr[0] = 1;
            kr[1] = 0;
            return kr;
        }
        if (Se < epsKr_)
        {
            kr[0] = 0;
            kr[1] = 1;
            return kr;
        }

        double lambda = soil_.paramRelPerm(globalPos, element, localPos, temperature)[0];
        double exponent = (2. + 3*lambda) / lambda;
        kr[0] = pow(Se, exponent);
        exponent = (2. + lambda) / lambda;
        kr[1] = pow(1.-Se, 2) * ( 1. - pow(Se, exponent) );
        return kr;
    }

    BrooksCoreyLaw(const Matrix2p<Grid,double>& soil, bool lin = false)
        : RelPerm_pc<Grid>(soil, false), soil_(soil)
    {
    }

private:
    double brooksCorey_(double Swe, double p0, double lambda) const
    {
        return p0*pow(Swe, -1.0/lambda);
    }

    static const double SweMin_ = 0.10;
    static const double epsKr_ = 1e-15;
    const Matrix2p<Grid,double>& soil_;
};
}

#endif
