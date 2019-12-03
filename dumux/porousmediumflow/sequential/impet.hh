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
#ifndef DUMUX_IMPET_HH
#define DUMUX_IMPET_HH

/**
 * @file
 * @brief  IMPET scheme
 */

#include "impetproperties.hh"

namespace Dumux
{
/**
 * \ingroup IMPET
 * \brief IMplicit Pressure Explicit Transport (IMPET) scheme for the solution of weakly coupled diffusion-transport formulations.
 *
 * The model implements the sequential equations of two-phase flow.
 * These equations can be derived from the two-phase flow equations shown for the two-phase box model (TwoPBoxModel).
 * The first equation to solve is a pressure equation of elliptic character.
 * The second one is a transport equation (e.g. for saturation, concentration,...), which can be hyperbolic or parabolic.
 *
 * As the equations are only weakly coupled they do not have to be solved simultaneously
 * but can be solved sequentially. First the pressure equation is solved implicitly,
 * second the transport equation can be solved explicitly. This solution procedure is called IMPES algorithm
 * (IMplicit Pressure Explicit Saturation) for immiscible flow or IMPEC algorithm
 * (IMplicit Pressure Explicit Concentration) for miscible flow.
 *
 * \tparam TypeTag The Type Tag
 */
template<class TypeTag> class IMPET
{
    using Problem = GetPropType<TypeTag, Properties::Problem>;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;

    using SolutionTypes = GetProp<TypeTag, Properties::SolutionTypes>;
    using ElementMapper = typename SolutionTypes::ElementMapper;
    using TransportSolutionType = GetPropType<TypeTag, Properties::TransportSolutionType>;

    enum IterationType
        {
            noIter,
            iterToNumIter,
            iterToConverged,
        };

public:
    using SolutionType = typename SolutionTypes::ScalarSolution;

    //! Set initial solution and initialize parameters
    void initialize()
    {
        //initial values of transported quantity
        problem_.transportModel().initialize();
        //call function with true to get a first initialisation of the pressure field
        problem_.pressureModel().initialize();

        //update constitutive functions
        problem_.pressureModel().updateMaterialLaws();

        Dune::dinfo << "IMPET: initialization done." << std::endl;

        return;
    }

    //! Calculate the update.
    /*!
     *  \param  t         time
     *  \param dt         time step size
     *  \param updateVec  vector for the update values
     *
     *  Calculates the new pressure and velocity and determines the time step size
     *  and the update of the transported quantity for the explicit time step.
     */
    void update(const Scalar t, Scalar& dt, TransportSolutionType& updateVec)
    {
        if (iterFlag_ == noIter)
        {
            //update pressure
            problem_.pressureModel().update();

            //calculate defect of transported quantity
            problem_.transportModel().update(t, dt, updateVec, true);

            dt *= cFLFactor_;
        }
        else if (iterFlag_ == iterToNumIter || iterFlag_ == iterToConverged)
        {
            bool converg = false;
            int iter = 0;
            int iterTot = 0;

            // the method is valid for any transported quantity.
            TransportSolutionType transValueOldIter;
            problem_.transportModel().getTransportedQuantity(transValueOldIter);
            TransportSolutionType updateOldIter(transValueOldIter);
            updateOldIter = 0;
            TransportSolutionType transportedQuantity(transValueOldIter);
            TransportSolutionType updateHelp(transValueOldIter);
            TransportSolutionType updateDiff(transValueOldIter);

            while (!converg)
            {
                iter++;
                iterTot++;

                problem_.pressureModel().update();

                //calculate defect of transported quantity
                problem_.transportModel().update(t, dt, updateVec, true);

                updateHelp = updateVec;
                problem_.transportModel().getTransportedQuantity(transportedQuantity);
                transportedQuantity += (updateHelp *= (dt * cFLFactor_));
                transportedQuantity *= omega_;
                transValueOldIter *= (1 - omega_);
                transportedQuantity += transValueOldIter;
                updateDiff = updateVec;
                updateDiff -= updateOldIter;
                transValueOldIter = transportedQuantity;
                updateOldIter = updateVec;
                Dune::dinfo << " defect = " << dt * updateDiff.two_norm() / transportedQuantity.two_norm();
                // break criteria for iteration loop
                if (iterFlag_ == iterToConverged && dt * updateDiff.two_norm() / transportedQuantity.two_norm() <= maxDefect_)
                {
                    converg = true;
                }
                else if (iterFlag_ == iterToNumIter && iter > nIter_)
                {
                    converg = true;
                }

                if (iterFlag_ == iterToConverged && transportedQuantity.infinity_norm() > (1 + maxDefect_))
                {
                    converg = false;
                }
                if (!converg && iter > nIter_)
                {
                    converg = true;
                    std::cout << "Nonlinear loop in IMPET.update exceeded nIter = " << nIter_ << " iterations."
                              << std::endl;
                    std::cout << transportedQuantity.infinity_norm() << std::endl;
                }
            }
            // outputs
            if (iterFlag_ == iterToConverged)
                std::cout << "Iteration steps: " << iterTot << std::endl;
            std::cout.setf(std::ios::scientific, std::ios::floatfield);

            dt *= cFLFactor_;
        }
        else
        {
            DUNE_THROW(Dune::NotImplemented,"IMPET: Iterationtype not implemented!");
        }

    }

    void updateTransport(const Scalar t, Scalar& dt, TransportSolutionType& updateVec)
    {
        problem_.pressureModel().updateVelocity();

        problem_.transportModel().update(t, dt, updateVec, true);
    }

    //! \brief Write data files
    /*!
     *  calls the output methods of both models, pressure and transport.
     *  \param writer the current VTKwriter
     */
    template<class MultiWriter>
    void addOutputVtkFields(MultiWriter &writer)
    {
        problem_.pressureModel().addOutputVtkFields(writer);
        problem_.transportModel().addOutputVtkFields(writer);

        return;
    }

    /*!
     * \brief Mapper for the entities where degrees of freedoms are defined.
     */
    const ElementMapper &dofMapper() const
    {
        return problem_.elementMapper();
    }

    //! Constructs an IMPET object
    /**
     * \param problem Problem
     */
    IMPET(Problem& problem) :
        problem_(problem)
    {
        cFLFactor_ = getParam<Scalar>("Impet.CFLFactor");
        iterFlag_ = getParam<int>("Impet.IterationFlag", 0);
        nIter_ = getParam<int>("Impet.IterationNumber", 2);
        maxDefect_ = getParam<Scalar>("Impet.MaximumDefect", 1e-5);
        omega_ = getParam<Scalar>("Impet.RelaxationFactor", 1.0);
    }

private:
    Problem& problem_; //!< object of type Problem including problem
    Scalar cFLFactor_;
    int iterFlag_; //!< flag to switch the iteration type of the IMPET scheme
    int nIter_; //!< number of iterations if IMPET iterations are enabled by IterationFlag
    Scalar maxDefect_; //!< maximum defect if IMPET iterations are enabled by IterationFlag
    Scalar omega_; //!< 1 = new solution is new solution, 0 = old solution is new solution
};
}
#endif
