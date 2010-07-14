// $Id: impes.hh 3784 2010-06-24 13:43:57Z bernd $
/*****************************************************************************
 *   Copyright (C) 2007-2009 by Bernd Flemisch                               *
 *   Copyright (C) 2008-2009 by Markus Wolff                                 *
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
#ifndef DUMUX_IMPES_HH
#define DUMUX_IMPES_HH

/**
 * @file
 * @brief  IMPES scheme
 * @author Bernd Flemisch, Markus Wolff
 */

#include "impesproperties.hh"

namespace Dumux
{
/**
 * \ingroup impes
 * \brief IMplicit Pressure Explicit Saturation (IMPES) scheme for the solution of weakly coupled diffusion/transport problems.
 *
 * The model implements the decoupled equations of two-phase flow of two completely immiscible fluids.
 * These equations can be derived from the two-phase flow equations shown for the two-phase box model (TwoPBoxModel).
 * The first equation to solve is a pressure equation of elliptic character. The second one is a saturation equation,
 * which can be hyperbolic or parabolic.
 *
 * This model allows different combinations of primary variables, which can be \f$p_w\f$-\f$S_w\f$, \f$p_w\f$-\f$S_n\f$, \f$p_n\f$-\f$S_w\f$, \f$p_n\f$-\f$S_n\f$,
 * or \f$p\f$-\f$S_w\f$ and \f$p\f$-\f$S_n\f$, where \f$p\f$ is no phase pressure but a global pressure.
 *
 * As the equations are only weakly coupled they do not have to be solved simultaneously
 * but can be solved sequentially. First the pressure equation is solved implicitly,
 * second the saturation equation can be solved explicitly. This solution procedure is called IMPES algorithm
 * (IMplicit Pressure Explicit Saturation).
 *
 * In comparison to a fully coupled model, different discretization methods can be applied to the different equations.
 * So far, the pressure equation is discretized using a cell centered finite volume scheme (optionally with multi point flux approximation),
 * a mimetic finite difference scheme or a finite element scheme. The saturation equation is discretized using a cell centered finite volume scheme.
 * Default time discretization scheme is an explicit Euler scheme.
 */
template<class TypeTag> class IMPES
{
    typedef IMPES<TypeTag> ThisType;
    typedef typename GET_PROP_TYPE(TypeTag, PTAG(GridView)) GridView;

    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Problem)) Problem;


    typedef typename GET_PROP_TYPE(TypeTag, PTAG(Scalar)) Scalar;

    typedef typename GET_PROP(TypeTag, PTAG(SolutionTypes)) SolutionTypes;
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
    virtual void initial()
    {
        //initial saturations
        problem.saturationModel().initialTransport();
        //call function with true to get a first initialisation of the pressure field
        problem.pressureModel().initial();
        problem.pressureModel().calculateVelocity();

        return;
    }

    //! Calculate the update.
    /*!
     *  \param  t         time
     *  \param dt         time step size
     *  \param updateVec  vector for the update values
     *  \param CLFFac     security factor for the time step criterion (0 < CLFFac <= 1)
     *
     *  Calculates the new pressure and velocity and determines the time step size and the saturation update for the explicit time step
     *  Called from Dumux::Timestep.
     */
    void update(const Scalar t, Scalar& dt, ScalarSolutionType& updateVec)
    {
        int satSize = problem.variables().gridSize();
        ScalarSolutionType saturation(problem.variables().saturation());
        ScalarSolutionType satOldIter(problem.variables().saturation());
        ScalarSolutionType satHelp(satSize);
        ScalarSolutionType satDiff(satSize);
        ScalarSolutionType updateOldIter(satSize);
        ScalarSolutionType updateHelp(satSize);
        ScalarSolutionType updateDiff(satSize);

        //update constitutive functions
        problem.pressureModel().updateMaterialLaws();

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

            //calculate saturation defect
            problem.saturationModel().update(t, dt, updateVec, cFLFactor_, true);

            if (iterFlag_)
            { // only needed if iteration has to be done
                updateHelp = updateVec;
                saturation = problem.variables().saturation();
                saturation += (updateHelp *= (dt * cFLFactor_));
                saturation *= omega_;
                satHelp = satOldIter;
                satHelp *= (1 - omega_);
                saturation += satHelp;
                updateDiff = updateVec;
                updateDiff -= updateOldIter;
                satOldIter = saturation;
                updateOldIter = updateVec;
                //                problem.saturationModel().updateMaterialLaws(saturation, true);
            }
            // break criteria for iteration loop
            if (iterFlag_ == 2 && dt * updateDiff.two_norm() / saturation.two_norm() <= maxDefect_)
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
            if (iterFlag_ == 2 && saturation.infinity_norm() > (1 + maxDefect_))
            {
                converg = false;
            }
            if (!converg && iter > nIter_)
            {
                std::cout << "Nonlinear loop in IMPES.update exceeded nIter = " << nIter_ << " iterations." << std::endl;
                std::cout << saturation.infinity_norm() << std::endl;
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
     *  \param name file name
     *  \param k format parameter
     */
    template<class MultiWriter>
    void addOutputVtkFields(MultiWriter &writer)
    {
        problem.variables().addOutputVtkFields(writer);
        return;
    }

    // serialization methods
    template<class Restarter>
    void serialize(Restarter &res)
    {
        problem.variables().serialize<Restarter> (res);
    }
    template<class Restarter>
    void deserialize(Restarter &res)
    {
        problem.variables().deserialize<Restarter> (res);
    }

    //! Constructs an IMPES object
    /**
     * \param prob Problem
     */
    IMPES(Problem& prob) :
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
