// -*- mode: C++; tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 4 -*-
// vi: set et ts=4 sw=4 sts=4:
/*****************************************************************************
 *   Copyright (C) 2010 by Markus Wolff, Benjamin Faigle                     *
 *   Institute for Modelling Hydraulic and Environmental Systems             *
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
#ifndef DUMUX_FVTRANSPORT2P2C_MULTIPHYSICS_HH
#define DUMUX_FVTRANSPORT2P2C_MULTIPHYSICS_HH

#include <dumux/decoupled/2p2c/fvtransport2p2c.hh>

/**
 * @file
 * @brief  Finite Volume discretization of the component transport equation
 * @author Markus Wolff, Jochen Fritz, Benjamin Faigle
 */

namespace Dumux
{
//! Miscible Transport step in a Finite Volume discretization
/*!
 * \ingroup multiphysics
 *  The finite volume model for the solution of the transport equation for compositional
 *  two-phase flow.
 *  \f[
      \frac{\partial C^\kappa}{\partial t} = - \nabla \cdot \left( \sum_{\alpha} X^{\kappa}_{\alpha} \varrho_{alpha} \bf{v}_{\alpha}\right) + q^{\kappa},
 *  \f]
 *  where \f$ \bf{v}_{\alpha} = - \lambda_{\alpha} \bf{K} \left(\nabla p_{\alpha} + \rho_{\alpha} \bf{g} \right) \f$.
 *  \f$ p_{\alpha} \f$ denotes the phase pressure, \f$ \bf{K} \f$ the absolute permeability, \f$ \lambda_{\alpha} \f$ the phase mobility,
 *  \f$ \rho_{\alpha} \f$ the phase density and \f$ \bf{g} \f$ the gravity constant and \f$ C^{\kappa} \f$ the total Component concentration.
 *
 * The model domain is automatically divided
 * in a single-phase and a two-phase domain. The full 2p2c model is only evaluated within the
 * two-phase subdomain, whereas a single-phase transport model is computed in the rest of the
 * domain.
 *
 *  \tparam TypeTag The Type Tag
 */
template<class TypeTag>
class FVTransport2P2CMultiPhysics : public FVTransport2P2C<TypeTag>
{
    typedef typename GET_PROP_TYPE(TypeTag, GridView) GridView;
    typedef typename GET_PROP_TYPE(TypeTag, Scalar) Scalar;
    typedef typename GET_PROP_TYPE(TypeTag, Problem) Problem;

    typedef typename GET_PROP_TYPE(TypeTag, SpatialParameters) SpatialParameters;
    typedef typename SpatialParameters::MaterialLaw MaterialLaw;

    typedef typename GET_PROP_TYPE(TypeTag, Indices) Indices;
    typedef typename GET_PROP_TYPE(TypeTag, FluidSystem) FluidSystem;
    typedef typename GET_PROP_TYPE(TypeTag, FluidState) FluidState;

    typedef typename GET_PROP_TYPE(TypeTag, CellData) CellData;

    typedef typename GET_PROP_TYPE(TypeTag, TransportSolutionType) TransportSolutionType;

    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };
    enum
    {
        pw = Indices::pressureW,
        pn = Indices::pressureNW,
        pglobal = Indices::pressureGlobal,
        vw = Indices::velocityW,
        vn = Indices::velocityNW,
        vt = Indices::velocityTotal,
        Sw = Indices::saturationW,
        Sn = Indices::saturationNW
    };
    enum
    {
        wPhaseIdx = Indices::wPhaseIdx, nPhaseIdx = Indices::nPhaseIdx,
        wCompIdx = Indices::wPhaseIdx, nCompIdx = Indices::nPhaseIdx
    };

    typedef typename GridView::template Codim<0>::Iterator ElementIterator;
    typedef typename GridView::IntersectionIterator IntersectionIterator;

    typedef Dune::FieldVector<Scalar, dimWorld> GlobalPosition;
    typedef Dune::FieldVector<Scalar, 2> PhaseVector;
    typedef typename GET_PROP_TYPE(TypeTag, PrimaryVariables) PrimaryVariables;

    //! Acess function for the current problem
    Problem& problem()
    {return this->problem_;};
public:
    virtual void update(const Scalar t, Scalar& dt, TransportSolutionType& updateVec, bool impes);

    //! Constructs a FVTransport2P2CMultiPhysics object
    /**
     * \param problem a problem class object
     */

    FVTransport2P2CMultiPhysics(Problem& problem) : FVTransport2P2C<TypeTag>(problem)
    {}

    virtual ~FVTransport2P2CMultiPhysics()
    {     }
};

//! \brief Calculate the update vector and determine timestep size
/*!
 *  This method calculates the update vector \f$ u \f$ of the discretized equation
 *  \f[
       C^{\kappa , new} = C^{\kappa , old} + u,
 *  \f]
 *  where \f$ u = \sum_{element faces} \boldsymbol{v}_{\alpha} * \varrho_{\alpha} * X^{\kappa}_{\alpha} * \boldsymbol{n} * A_{element face} \f$,
 *  \f$ \boldsymbol{n} \f$ is the face normal and \f$ A_{element face} \f$ is the face area.
 *
 *  In addition to the \a update vector, the recommended time step size \a dt is calculated
 *  employing a CFL condition.
 *  This method = old concentrationUpdate()
 *
 *  \param t Current simulation time \f$\mathrm{[s]}\f$
 *  \param[out] dt Time step size \f$\mathrm{[s]}\f$
 *  \param[out] updateVec Update vector, or update estimate for secants, resp. Here in \f$\mathrm{[kg/m^3]}\f$
 *  \param impet Flag that determines if it is a real impet step or an update estimate for volume derivatives
 */
template<class TypeTag>
void FVTransport2P2CMultiPhysics<TypeTag>::update(const Scalar t, Scalar& dt, TransportSolutionType& updateVec, bool impet = false)
{
    // initialize dt very large
    dt = 1E100;

    // resize update vector and set to zero
    updateVec.resize(GET_PROP_VALUE(TypeTag, NumComponents));
    updateVec[wCompIdx].resize(problem().gridView().size(0));
    updateVec[nCompIdx].resize(problem().gridView().size(0));
    updateVec[wCompIdx] = 0;
    updateVec[nCompIdx] = 0;

    // Cell which restricts time step size
    int restrictingCell = -1;

    PhaseVector entries(0.), timestepFlux(0.);
    // compute update vector
    ElementIterator eItEnd = problem().gridView().template end<0> ();
    for (ElementIterator eIt = problem().gridView().template begin<0> (); eIt != eItEnd; ++eIt)
    {
        // get cell infos
        int globalIdxI = problem().variables().index(*eIt);
        CellData& cellDataI = problem().variables().cellData(globalIdxI);

        if(impet or cellDataI.subdomain()==2)   // estimate only necessary in subdomain
        {
            // some variables for time step calculation
            double sumfactorin = 0;
            double sumfactorout = 0;

            // run through all intersections with neighbors and boundary
            IntersectionIterator isItEnd = problem().gridView().iend(*eIt);
            for (IntersectionIterator isIt = problem().gridView().ibegin(*eIt); isIt != isItEnd; ++isIt)
            {

                // handle interior face
                if (isIt->neighbor())
                {
                    this->getFlux(entries, timestepFlux, *isIt, cellDataI);
                }

                /******************************************
                 *     Boundary Face
                 ******************************************/
                if (isIt->boundary())
                {
                    this->getFluxOnBoundary(entries, timestepFlux, *isIt, cellDataI);
                }
                if(isnan(entries[0]) or isinf(entries[0]) or isnan(entries[1]) or isinf(entries[1]))
                    std::cout << "args!!!";

                // add to update vector
                updateVec[wCompIdx][globalIdxI] += entries[wCompIdx];
                updateVec[nCompIdx][globalIdxI] += entries[nCompIdx];

                // for time step calculation
                sumfactorin += timestepFlux[0];
                sumfactorout += timestepFlux[1];

            }// end all intersections

            /***********     Handle source term     ***************/
            PrimaryVariables q(NAN);
            problem().source(q, *eIt);
            updateVec[wCompIdx][globalIdxI] += q[Indices::contiWEqIdx];
            updateVec[nCompIdx][globalIdxI] += q[Indices::contiNEqIdx];

            // account for porosity
            sumfactorin = std::max(sumfactorin,sumfactorout)
                            / problem().spatialParameters().porosity(*eIt);

            if ( 1./sumfactorin < dt)
            {
                dt = 1./sumfactorin;
                restrictingCell= globalIdxI;
            }
        }
    } // end grid traversal
    if(impet)
        Dune::dinfo << "Timestep restricted by CellIdx " << restrictingCell << " leads to dt = "<<dt * GET_PARAM(TypeTag, Scalar, CFLFactor)<< std::endl;
    return;
}
}
#endif
