// $Id:$
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
#ifndef DUNE_IMPES_HH
#define DUNE_IMPES_HH

#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include "dumux/fractionalflow/fractionalflow.hh"

/**
 * @file
 * @brief  IMPES scheme
 * @author Bernd Flemisch, Markus Wolff
 */

namespace Dune
{
/**
 * \ingroup impes
 * @brief IMplicit Pressure Explicit Saturation (IMPES) scheme for the solution of
 * weakly coupled diffusion/transport problems.
 * First a pressure equation is solved implicitly to obtain a velocity field, then a saturation equation is solved
 * explicitly to get a new saturation distribution.
 *
 *  Template parameters are:

 - GridView      a DUNE gridview type
 - Diffusion     class defining the diffusion model
 - Transport     class defining the transport model
 - VC            type of a class containing different variables of the model
*/
template<class GridView, class Diffusion, class Transport, class VC> class IMPES: public FractionalFlow<
        GridView, Diffusion, Transport, VC>
{
    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };
    typedef typename GridView::Grid Grid;
    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef Dune::FractionalFlow<GridView, Diffusion, Transport, VC> FractionalFlow;
    typedef typename FractionalFlow::RepresentationType PressType;
    typedef typename Diffusion::ScalarType Scalar;

    typedef Dune::FieldVector<Scalar,dim> LocalPosition;
    typedef Dune::FieldVector<Scalar,dimWorld> GlobalPosition;

public:
    typedef typename Transport::RepresentationType RepresentationType;//!< Data type for a Vector of Scalars

    //! Set initial solution and initialize parameters
    virtual void initial()
    {
        Scalar t = 0;
        //initial saturations
        this->transport.initialTransport();
        //call function with true to get a first initialisation of the pressure field
        this->diffusion.pressure(true,t);
        this->diffusion.calculateVelocity(t);

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
     *  Called from Dune::Timestep.
     */
    virtual int update(const Scalar t, Scalar& dt, RepresentationType& updateVec,
            Scalar cFLFactor = 1)
    {
        int satSize = this->transport.problem().variables().gridSizeTransport();
        RepresentationType saturation(this->transport.problem().variables().saturation());
        RepresentationType satOldIter(this->transport.problem().variables().saturation());
        RepresentationType satHelp(satSize);
        RepresentationType satDiff(satSize);
        RepresentationType updateOldIter(satSize);
        RepresentationType updateHelp(satSize);
        RepresentationType updateDiff(satSize);

        //update constitutive functions
        this->transport.updateMaterialLaws();

        bool converg = false;
        int iter = 0;
        int iterTot = 0;
        updateOldIter = 0;
        while (!converg)
        {
            iter++;
            iterTot++;

            // update pressure: give false as the pressure field is already initialised
            this->diffusion.pressure(false , t);

            //calculate velocities
            this->diffusion.calculateVelocity(t);

            //calculate saturation defect
            this->transport.update(t, dt, updateVec,cFLFactor, true);

            if (iterFlag)
            { // only needed if iteration has to be done
                updateHelp = updateVec;
                saturation = this->transport.problem().variables().saturation();
                saturation += (updateHelp *= (dt*cFLFactor));
                saturation *= omega;
                satHelp = satOldIter;
                satHelp *= (1-omega);
                saturation += satHelp;
                updateDiff = updateVec;
                updateDiff -= updateOldIter;
                satOldIter = saturation;
                updateOldIter = updateVec;
                this->transport.updateMaterialLaws(saturation, true);
            }
            // break criteria for iteration loop
            if (iterFlag == 2 && dt * updateDiff.two_norm() / saturation.two_norm() <= maxDefect )
            {
                converg = true;
            }
            else if (iterFlag == 1 && iter > nIter )
            {
                converg = true;
            }
            else if (iterFlag==0)
            {
                converg = true;
            }
            if (iterFlag == 2 && saturation.infinity_norm() > (1+maxDefect))
            {
                converg = false;
            }
            if (!converg && iter > nIter )
            {
                std::cout << "Nonlinear loop in IMPES.update exceeded nIter = "
                << nIter << " iterations."<< std::endl;
                std::cout<<saturation.infinity_norm()<<std::endl;
                return 1;
            }
        }
        // outputs
        if (iterFlag==2)
        std::cout << "Iteration steps: "<< iterTot << std::endl;
        std::cout.setf(std::ios::scientific, std::ios::floatfield);

        return 0;
    }

    //! \brief Write data files
    /*!
     *  \param name file name
     *  \param k format parameter
     */
    virtual void vtkout(const char* name, int k) const
    {
        this->transport.problem().variables().vtkout(name, k);
        return;
    }

    //! Constructs an IMPES object
    /**
     * \param diffusion an object of type Diffusion
     * \param transport an object of type Transport
     * \param flag iteration flag (0 -> no iterations, 1 -> iterate until max number of iterations is reached, 2 -> iterate until converged or max number of iterations is reached)
     * \param nIt maximum number of iterations
     * \param maxDef maximum defect for convergence criterion
     * \param om under relaxation factor (om <= 1)
     */
    IMPES(Diffusion& diffusion, Transport& transport, int flag = 0, int nIt = 2,
            Scalar maxDef = 1e-5, Scalar om = 1) :
    FractionalFlow(diffusion, transport),
    iterFlag(flag), nIter(nIt), maxDefect(maxDef), omega(om)
    {
    }

protected:
    const int iterFlag;//!<iteration flag
    const int nIter;//!<maximum number of iterations
    const Scalar maxDefect;//!maximum defect for convergence criterion
    const Scalar omega;//!under relaxation factor
};
}
#endif
