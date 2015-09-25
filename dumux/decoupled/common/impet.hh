// $Id$
/*****************************************************************************
 *   Copyright (C) 2007-2009 by Bernd Flemisch                               *
 *   Copyright (C) 2008-2009 by Markus Wolff                                 *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: <givenname>.<name>@iws.uni-stuttgart.de                          *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
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
 * @author Bernd Flemisch, Markus Wolff
 */

#include "impetproperties.hh"

namespace Dumux
{
/**
 * \ingroup IMPET
 * \brief IMplicit Pressure Explicit Transport (IMPET) scheme for the solution of weakly coupled diffusion-transport formulations.
 *
 * The model implements the decoupled equations of two-phase flow.
 * These equations can be derived from the two-phase flow equations shown for the two-phase box model (TwoPBoxModel).
 * The first equation to solve is a pressure equation of elliptic character. The second one is a transport equation (e.g. for saturation, concentration,...),
 * which can be hyperbolic or parabolic.
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
    typedef IMPET<TypeTag> ThisType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;


    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;

    typedef typename GET_PROP(TypeTag, PTAG(SolutionTypes)) SolutionTypes;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(TransportSolutionType)) TransportSolutionType;
    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };

    typedef typename SolutionTypes::ScalarSolution ScalarSolutionType;

    typedef Dune::FieldVector<Scalar, dim> LocalPosition;
    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;

public:
    typedef typename SolutionTypes::ScalarSolution SolutionType;

    //! Set initial solution and initialize parameters
    virtual void initialize()
    {
        //initial values of transported quantity
        problem.transportModel().initialize();
        //call function with true to get a first initialisation of the pressure field
        problem.pressureModel().initialize();
        problem.pressureModel().calculateVelocity();

        //update constitutive functions
        problem.pressureModel().updateMaterialLaws();

        Dune::dinfo << "IMPET: initialization done." << std::endl;

        return;
    }

    //! Calculate the update.
    /*!
     *  \param  t         time
     *  \param dt         time step size
     *  \param updateVec  vector for the update values
     *
     *  Calculates the new pressure and velocity and determines the time step size and the update of the transported quantity for the explicit time step
     */
    void update(const Scalar t, Scalar& dt, TransportSolutionType& updateVec)
    {
        // the method is valid for any transported quantity.
        int transSize = problem.variables().gridSize();
        TransportSolutionType transportedQuantity(problem.variables().transportedQuantity());
        TransportSolutionType transValueOldIter(problem.variables().transportedQuantity());
        TransportSolutionType transValueHelp(transSize);
        TransportSolutionType transValueDiff(transSize);
        TransportSolutionType updateOldIter(transSize);
        TransportSolutionType updateHelp(transSize);
        TransportSolutionType updateDiff(transSize);

        bool converg = false;
        int iter = 0;
        int iterTot = 0;
        updateOldIter = 0;

        while (!converg)
        {
            iter++;
            iterTot++;

            problem.pressureModel().pressure(false);

            //calculate velocities
            problem.pressureModel().calculateVelocity();

            //calculate defect of transported quantity
            problem.transportModel().update(t, dt, updateVec, true);

            if (iterFlag_)
            { // only needed if iteration has to be done
                updateHelp = updateVec;
                transportedQuantity = problem.variables().transportedQuantity();
                transportedQuantity += (updateHelp *= (dt * cFLFactor_));
                transportedQuantity *= omega_;
                transValueHelp = transValueOldIter;
                transValueHelp *= (1 - omega_);
                transportedQuantity += transValueHelp;
                updateDiff = updateVec;
                updateDiff -= updateOldIter;
                transValueOldIter = transportedQuantity;
                updateOldIter = updateVec;
            }
            // break criteria for iteration loop
            if (iterFlag_ == 2 && dt * updateDiff.two_norm() / transportedQuantity.two_norm() <= maxDefect_)
            {
                converg = true;
            }
            else if (iterFlag_ == 1 && iter > nIter_)
            {
                converg = true;
            }
            else if (iterFlag_ == 0)
            {
                converg = true;
            }
            if (iterFlag_ == 2 && transportedQuantity.infinity_norm() > (1 + maxDefect_))
            {
                converg = false;
            }
            if (!converg && iter > nIter_)
            {
                converg = true;
                std::cout << "Nonlinear loop in IMPET.update exceeded nIter = " << nIter_ << " iterations." << std::endl;
                std::cout << transportedQuantity.infinity_norm() << std::endl;
            }
        }
        // outputs
        if (iterFlag_ == 2)
            std::cout << "Iteration steps: " << iterTot << std::endl;
        std::cout.setf(std::ios::scientific, std::ios::floatfield);

        dt *= cFLFactor_;
    }

    //! \brief Write data files
    /*!
     *  calls the output methods of both models, pressure and transport.
     *  \param writer the current VTKwriter
     */
    template<class MultiWriter>
    void addOutputVtkFields(MultiWriter &writer)
    {
        problem.pressureModel().addOutputVtkFields(writer);
        problem.transportModel().addOutputVtkFields(writer);

        return;
    }

    // serialization methods
    //! Function needed for restart option.
    template<class Restarter>
    void serialize(Restarter &res)
    {
        problem.variables().serialize<Restarter> (res);
    }

    //! Function needed for restart option.
    template<class Restarter>
    void deserialize(Restarter &res)
    {
        problem.variables().deserialize<Restarter> (res);
        //update constitutive functions
        problem.pressureModel().updateMaterialLaws();
    }

    //! Constructs an IMPET object
    /**
     * \param prob Problem
     */
    IMPET(Problem& prob) :
        problem(prob)
    {
        cFLFactor_ = GET_PROP_VALUE(TypeTag, PTAG(CFLFactor));
        iterFlag_ = GET_PROP_VALUE(TypeTag, PTAG(IterationFlag));
        nIter_ = GET_PROP_VALUE(TypeTag, PTAG(IterationNumber));
        maxDefect_ = GET_PROP_VALUE(TypeTag, PTAG(MaximumDefect));
        omega_ = GET_PROP_VALUE(TypeTag, PTAG(RelaxationFactor));
    }

protected:

    Problem& problem; //!< object of type Problem including problem
private:
    Scalar cFLFactor_;
    int iterFlag_;
    int nIter_;
    Scalar maxDefect_;
    Scalar omega_;
};
}
#endif
