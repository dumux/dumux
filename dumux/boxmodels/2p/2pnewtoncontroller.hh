// $Id: 2pnewtoncontroller.hh 3738 2010-06-15 14:01:09Z lauser $
/*****************************************************************************
 *   Copyright (C) 2008-2009 by Andreas Lauser                               *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software; you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation; either version 2 of the License, or       *
 *   (at your option) any later version, as long as this copyright notice    *
 *   is included in its original form.                                       *
 *                                                                           *
 *   This program is distributed WITHOUT ANY WARRANTY.                       *
 *****************************************************************************/
/*!
 * \file
 * \brief A newton controller for two-phase problems.
 *
 * This controller 'knows' what a 'physically meaningful' solution is
 * which allows the newton method to abort quicker if the solution is
 * way out of bounds.
 */
#ifndef DUMUX_2P_NEWTON_CONTROLLER_HH
#define DUMUX_2P_NEWTON_CONTROLLER_HH

#include <dumux/nonlinear/newtoncontroller.hh>

namespace Dumux {
/*!
 * \ingroup TwoPBoxModel
 * \brief A newton controller for two-phase problems.
 *
 * This controller 'knows' what a 'physically meaningful' solution is
 * which allows the newton method to abort quicker if the solution is
 * way out of bounds.
 */
template <class TypeTag>
class TwoPNewtonController : public NewtonController<TypeTag>
{
    typedef TwoPNewtonController<TypeTag>  ThisType;
    typedef NewtonController<TypeTag>      ParentType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(NewtonController)) Implementation;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Model)) Model;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(NewtonMethod)) NewtonMethod;

    typedef typename GET_PROP(TypeTag, PTAG(SolutionTypes)) SolutionTypes;
    typedef typename SolutionTypes::SolutionVector SolutionVector;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TwoPIndices)) Indices;
    enum {
        numEq = GET_PROP_VALUE(TypeTag, PTAG(NumEq)),
        pressureIdx = Indices::pressureIdx
    };

public:
    TwoPNewtonController()
    {};

#if 0
    /*!
     * \brief Update the error of the solution compared to the
     *        previous iteration.
     */
    void newtonUpdateRelError(const SolutionVector &uOld,
                              const SolutionVector &deltaU)
    {
        typedef typename SolutionVector::block_type FV;

        // get the maximum value of each primary in either the old or
        // the new solution
        Dune::FieldVector<Scalar, numEq> weight(1.0);

        // a change in pressure is considered to be less severe by a
        // factor of 10000 than a change saturation.
        weight[pressureIdx] = 1;// 1e-4;

        // calculate the relative error as the maximum relative
        // deflection in any degree of freedom.
        this->error_ = 0;
        int offenderVertIdx = -1;
        int offenderPVIdx = -1;
        Scalar offenderDelta = 0;
        for (int i = 0; i < int(uOld.size()); ++i) {
            for (int j = 0; j < FV::size; ++j) {
                // calculate the relative error at the current vertex
                // i and the current primary variable j
//              Scalar curErr = std::abs(deltaU[i][j] * weight[j]);
              Scalar curErr = std::abs(deltaU[i][j]/(1.0 + std::abs(uOld[i][j])));
#if 0
                // make sure that the specified tolerance is not below
                // machine precision!
                typedef std::numeric_limits<Scalar> limits;
                const Scalar machinePrec =
                    limits::epsilon()*100
                    * std::max<Scalar>(Scalar(1e10)/limits::max(),
                                       std::abs(uOld[i][j]));
                if (this->tolerance_ < machinePrec*weight[j])
                {
                    std::cerr << "Allowed tolerance ("
                              << this->tolerance_
                              << ") is below machine precision ("
                              << machinePrec*weight[j]
                              << ") at vertex " << i
                              << ", primary var " << j << "!\n";
                    if (curErr < machinePrec*weight[j])
                        // if the machine precision is reached, we
                        // accept the solution even if we're above the
                        // tolerance!
                        curErr = this->tolerance_/100;
                };
#endif

                if (this->error_ < curErr) {
                    offenderVertIdx = i;
                    offenderPVIdx = j;
                    offenderDelta = deltaU[i][j];
                    this->error_ = curErr;
                };
            }
        };

        this->endIterMsg() << ", worst offender: vertex " << offenderVertIdx
                           << ", primary var " << offenderPVIdx
                           << ", delta: " << offenderDelta;

        this->model().gridView().comm().max(this->error_);
    }
#endif
};
}

#endif
