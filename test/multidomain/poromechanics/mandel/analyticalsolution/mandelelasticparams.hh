// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
//
// SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
// SPDX-License-Identifier: GPL-3.0-or-later
//
//
// Created by yue on 7/25/23.
//

#ifndef MANDEL_MATERIAL_PORO_ELASTIC_PARAMETERS_HH
#define MANDEL_MATERIAL_PORO_ELASTIC_PARAMETERS_HH

#include <dumux/common/parameters.hh>

namespace Dumux{
/*!
 * \ingroup
 * \brief Auxilary structure for poro elastic parameters within the framework of
 * linear poro-elasticity.
 *
 * \tparam Scalar
 */
template<class Scalar>
struct MandelPoroElasticParameters{
    MandelPoroElasticParameters(){
        E_ = getParam<Scalar>("MaterialParameters.EModulus");
        nu_ = getParam<Scalar>("MaterialParameters.Nu");
        lambda_ = (E_ * nu_) / ((1+nu_) * (1-2*nu_));
        G_ = E_ / 2 /(1+nu_);
        Kb_ = E_ / 3 / (1 - 2*nu_);

        alpha_ = getParam<Scalar>("MaterialParameters.Alpha");
        M_ = getParam<Scalar>("MaterialParameters.MModulus");

        Kbu_ = Kb_ + alpha_*alpha_*M_;
        B_ = alpha_ * M_ / Kbu_;
        nuu_ = (3*Kbu_ - 2*G_) /2 /(3*Kbu_ + G_);
    };

    /*!
     * \brief return first lame parameter
     *
     * \return Scalar
     */
    const Scalar lambda() const
    {return lambda_; }

    /*!
     * \brief return the second lame parameter
     *
     * \return Scalar
     */
    const Scalar Nu() const
    {return nu_;}

    /*!
     * \brief return E Modulus
     *
     * \return Scalar
     */
    const Scalar E() const
    {return E_;}

    /*!
     * \brief return shear modulus
     *
     * \return Scalar
     */
    const Scalar G() const
    {return G_;}

    /*!
     * \brief return bulk modulus
     *
     * \return Scalar
     */
    const Scalar Kb() const
    {return Kb_;}

    /*!
     * \brief return Biot's coefficient
     *
     * \return Scalar
     */
    const Scalar& alpha() const
    {return alpha_;}

    /*!
     * \brief return Biot's Modulus
     *
     * \return Scalar
     */
    const Scalar& M() const
    {return M_;}

    /*!
     * \brief return Skempton pore pressure coeffcient
     *
     * \return Scalar
     */
    Scalar B() const
    {return B_;}

    /*!
     * \brief return undrained Possion's ratio
     *
     * \return Scalar
     */
    Scalar Nuu() const
    {return nuu_;}

    /*!
     * \brief return undrained bulk modulus
     *
     * \return Scalar
     */
    Scalar Kbu() const
    {return Kbu_;}

private:
    Scalar E_; // E modulus
    Scalar nu_; // Poisson's ratio
    Scalar G_; // shear modulus
    Scalar Kb_; // Bulk modulus
    Scalar lambda_; //first lame parameter

    Scalar alpha_; // Biot coefficient
    Scalar M_; // Biot Modulus

    Scalar Kbu_; // Undrained bulk modulus
    Scalar B_; // Skempton's coefficient
    Scalar nuu_; // Undrained Poisson's ratio
};
}//end namespace dumux
#endif
