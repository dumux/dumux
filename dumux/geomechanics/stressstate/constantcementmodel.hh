// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   See the file COPYING for full copying permissions.                      *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 3 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the            *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 * \ingroup Components
 * \brief A much simpler (and thus potentially less buggy) version of
 *        pure water.
 */
#ifndef DUMUX_GEOMECHANICS_CONSTANT_CEMENT_MODEL_HH
#define DUMUX_GEOMECHANICS_CONSTANT_CEMENT_MODEL_HH


#include <cmath>
#include <dune/common/exceptions.hh>
#include <dumux/common/parameters.hh>
#include <dumux/geomechanics/lameparams.hh>

namespace Dumux {

/*!
 * \ingroup
 * \brief class for constant-cement model
 *        which calculates the effective mechanical moduli
 *        with two components properties and porosity.
 *        For details, check paper
 *        Per Avseth 2010, http://dx.doi.org/10.1190/1.3483770
 *
 * \tparam Scalar The type used for scalar values
 */
template <class Scalar>
class ConstantCementModel
{
    using LameParams = Dumux::LameParams<Scalar>;

public:
    /*!
     * \brief Construct a new Constant Cement Model object
     *
     */
    ConstantCementModel(const Scalar& phiCrit,
                        const Scalar& phib,
                        const Scalar& Ks,
                        const Scalar& Gs,
                        const Scalar& Kc,
                        const Scalar& Gc)
    : phiCrit_{phiCrit}
    , phib_{phib}
    , Ks_{Ks}
    , Gs_{Gs}
    , Kc_{Kc}
    , Gc_{Gc}
    {
        // calculate some fixed parameters
        nus_ = 0.5 * (Ks_/Gs_ - 2.0/3) / (Ks_/Gs_ + 1.0/3); //Eq. 29
        nuc_ = 0.5 * (Kc_/Gc_ - 2.0/3) / (Kc_/Gc_ + 1.0/3); //Eq. 28
        Mc_ = Kc_ + 4.0/3*Gc_; //Compressional modulus

        using std::pow;
        lambdan_ = 2 * Gc_ * (1 - nus_) * (1 - nuc_)
                 / M_PI / Gs_ / (1 - 2*nuc_); //Eq.25
        an_ = -0.024153 * pow(lambdan_, -1.3646); //Eq.18
        bn_ = 0.20405 * pow(lambdan_, -0.89008);  //Eq.19
        cn_ = 0.00024649 * pow(lambdan_, -1.9864);//Eq.20

        lambdaTau_ = Gc_ / M_PI / Gs_; // Eq.26
        aTau_ = - 1e-2 * (2.26*nus_* nus_ + 2.07*nus_ + 2.3)
              * pow(lambdaTau_, (0.079*nus_* nus_ + 0.1754*nus_ - 1.342)); // Eq.22
        bTau_ = (0.0573*nus_* nus_ + 0.0937*nus_ + 0.202)
              * pow(lambdaTau_, (0.0274*nus_* nus_ + 0.0529*nus_ - 0.8765)); //Eq.23
        cTau_ = 1e-4 * (9.654*nus_* nus_ + 4.945*nus_ + 3.1)
              * pow(lambdaTau_, (0.01867*  nus_* nus_ + 0.4011*nus_ - 1.8186)); //Eq.24

        // calculate moduli at critical porosity
        Scalar n = connectivity(phib_);
        Scalar vAlpha = alpha(phib_);
        Scalar vsn = sn(phib_, vAlpha);
        Scalar vsTau = sTau(phib_, vAlpha);
        Kb_ = contactCementBulkModulus(phib_, n, vsn);
        Gb_ = contactCementShearModulus(phib_, Kb_, n, vsTau);
    }

    /*!
     * \brief calculate the coordination number
     *
     * \param phi porosity
     * \return const Scalar
     */
    Scalar connectivity(const Scalar& phi) const
    {
        //return 9;
        return 20 - 34 * phi  + 14 * phi * phi; //Eq.10
    }

    /*!
     * \brief calculate parameter alpha
     *
     * \param phi
     * \return const Scalar
     */
    Scalar alpha(const Scalar& phi) const
    {
        using std::pow;
        return pow((2.0/3 * (phiCrit_-phi)) / (1-phiCrit_),0.5); //Eq. 27
    }

    /*!
     * \brief calculate parameter proportional to normal stiffness
     *
     * \param phi porosity
     * \param alpha
     * \return const Scalar
     */
    Scalar sn(const Scalar& phi, const Scalar& alpha) const
    {
        return an_ * alpha * alpha + bn_ * alpha + cn_; //Eq. 17
    }

    /*!
     * \brief calculate parameter proportional to tangential stiffness
     *
     * \param phi
     * \param alpha
     * \return const Scalar
     */
    Scalar sTau(const Scalar& phi, const Scalar& alpha) const
    {
        return aTau_ * alpha * alpha + bTau_ * alpha + cTau_; //Eq.21
    }

    /*!
     * \brief calculate bulk modulus after contact cement model
     * for porosity > critical porosity
     *
     * \param phi porosity
     * \param sn
     * \param n connectivity
     * \return const Scalar
     */
    Scalar contactCementBulkModulus(const Scalar& phi,
                                    const Scalar& n,
                                    const Scalar& sn) const
    {
        return n * (1-phiCrit_) * Mc_ * sn / 6; //Eq.15
    }

    /*!
     * \brief calculate shear modulus after contact cement model
     * for porosity > critical porosity
     *
     * \param phi porosity
     * \param Kb bulk modulus after contact cement model
     * \param n connectivity
     * \param sTau
     * \return const Scalar
     */
    Scalar contactCementShearModulus(const Scalar& phi,
                                    const Scalar& Kb,
                                    const Scalar& n,
                                    const Scalar& sTau) const
    {
        return 3 * Kb / 5 + 3 * n * (1-phiCrit_) * Gc_ * sTau / 20; //Eq. 16
    }

    LameParams effectiveLameModuli(const Scalar& phi) const
    {
        LameParams lameParams;
        Scalar K, G;

        if (phi>phiCrit_)
        {
            DUNE_THROW(Dune::RangeError,"porosity is greater than the critical porosity and the curve is undefined in this range.");
        }

        // return Kb Gb directly
        if(phi == phiCrit_ )
        {
            K = 0;
            G = 0;
        }
        else if (phi == phib_ )
        {
            K = Kb_;
            G = Gb_;
        }

        // using contact cement model
        else if (phi > phib_)
        {
            Scalar n = connectivity(phi);
            Scalar vAlpha = alpha(phi);
            Scalar vsn = sn(phi, vAlpha);
            Scalar vsTau = sTau(phi, vAlpha);
            K = contactCementBulkModulus(phi, n, vsn);
            G = contactCementShearModulus(phi, K, n, vsTau);
        }

        // using friable-sand model
        else
        {
            K = (phi/phib_)/(Kb_ + 4.0/3*Gb_)
              + (1 - phi/phib_)/(Ks_ + 4.0/3*Gb_);
            K = 1/K -4.0/3*Gb_;

            Scalar z = Gb_/6 * (9*Kb_ + 8*Gb_) / (Kb_ + 2*Gb_);

            G = (phi/phib_)/(Gb_ + z)
              + (1 - phi/phib_)/(Gs_+z);
            G = 1/G - z;
        }

        // std::cout<<"K = "<< K << "\n"
        //          <<"G = "<< G << std::endl;
        Scalar lambda = K - 2.0/3*G;
        lameParams.setLambda(lambda);
        lameParams.setMu(G);
        return lameParams;
    }

private:
    Scalar phiCrit_; //!< critical porosity for rock
    Scalar phib_; //!< well-sorted end-member porosity
    Scalar Kb_; //!< bulk modulus at critical porosity [Pa]
    Scalar Gb_; //!< shear modulus at critical porosity [Pa]

    Scalar Ks_; //!< bulk modulus of the grain material [Pa]
    Scalar Gs_; //!< shear modulus of the grain material [Pa]
    Scalar nus_; //!< Poisson's ratio of the grain material

    Scalar Kc_; //!< bulk modulus of the cement material [Pa]
    Scalar Gc_; //!< shear modulus of the cement material [Pa]
    Scalar nuc_; //!< Poisson's ratio of the cement material
    Scalar Mc_; //!< compressional modulus of the cement material [Pa]

    // parameters for calculation
    Scalar lambdan_;
    Scalar an_;
    Scalar bn_;
    Scalar cn_;
    Scalar lambdaTau_;
    Scalar aTau_;
    Scalar bTau_;
    Scalar cTau_;

};
} // end namespace Dumux

#endif
