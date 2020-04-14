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
 * \ingroup SequentialTwoPTwoCModel
 * \brief Finite volume discretization of the component transport equation.
 */
#ifndef DUMUX_FVTRANSPORT2P2C_MULTIPHYSICS_HH
#define DUMUX_FVTRANSPORT2P2C_MULTIPHYSICS_HH

#include <dumux/porousmediumflow/2p2c/sequential/fvtransport.hh>
#include <dumux/parallel/vectorcommdatahandle.hh>

namespace Dumux {
/*!
 * \ingroup SequentialTwoPTwoCModel
 * \brief Compositional transport step in a finite volume discretization.
 *
 *  The finite volume model for the solution of the transport equation for compositional
 *  two-phase flow.
 *  \f[
      \frac{\partial C^\kappa}{\partial t} =
      - \nabla \cdot \left( \sum_{\alpha} X^{\kappa}_{\alpha} \varrho_{alpha} \bf{v}_{\alpha}\right) + q^{\kappa},
 *  \f]
 *  where \f$ \bf{v}_{\alpha} = - \lambda_{\alpha} \bf{K} \left(\nabla p_{\alpha} + \rho_{\alpha} \bf{g} \right) \f$.
 *  \f$ p_{\alpha} \f$ denotes the phase pressure, \f$ \bf{K} \f$ the absolute permeability,
 *  \f$ \lambda_{\alpha} \f$ the phase mobility,
 *  \f$ \rho_{\alpha} \f$ the phase density and \f$ \bf{g} \f$ the gravity constant and
 *  \f$ C^{\kappa} \f$ the total Component concentration.
 *
 * The model domain is automatically divided into a single-phase and a two-phase domain. As the flux computation is relatively cheap,
 * the same method is used for the real transport step independently of the subdomain.
 * The pressure equation does not need any derivatives in simple
 * subdomains, therefore in the transport estimate step inter-cell fluxes in the simple subdomain are omitted.
 *
 *  \tparam TypeTag The Type Tag
 */
template<class TypeTag>
class FVTransport2P2CMultiPhysics : public FVTransport2P2C<TypeTag>
{
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;

    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;

    using CellData = GetPropType<TypeTag, Properties::CellData>;

    using TransportSolutionType = GetPropType<TypeTag, Properties::TransportSolutionType>;

    enum
    {
        dim = GridView::dimension
    };
    enum
    {
        wPhaseIdx = Indices::wPhaseIdx, nPhaseIdx = Indices::nPhaseIdx,
        wCompIdx = Indices::wPhaseIdx, nCompIdx = Indices::nPhaseIdx
    };

    using PhaseVector = Dune::FieldVector<Scalar, 2>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;

    //! Acess function for the current problem
    Problem& problem()
    { return this->problem_; }

    using LocalTimesteppingData = typename FVTransport2P2C<TypeTag>::LocalTimesteppingData;

public:
    virtual void update(const Scalar t, Scalar& dt, TransportSolutionType& updateVec, bool impet = false);

    //! Constructs a FVTransport2P2CMultiPhysics object
    /**
     * \param problem a problem class object
     */
    FVTransport2P2CMultiPhysics(Problem& problem) : FVTransport2P2C<TypeTag>(problem)
    {}

    virtual ~FVTransport2P2CMultiPhysics()
    {}
};

/*!
 * \brief Calculate the update vector and determine timestep size
 *
 *  This method calculates the update vector \f$ u \f$ of the discretized equation
 *  \f[
       C^{\kappa , new} = C^{\kappa , old} + u,
 *  \f]
 *  where \f$ u = \sum_{\gamma} \boldsymbol{v}_{\alpha} * \varrho_{\alpha} * X^{\kappa}_{\alpha} * \boldsymbol{n} * A_{\gamma} \f$,
 *  \f$ \boldsymbol{n} \f$ is the face normal and \f$ A_{\gamma} \f$ is the face area of face \f$ \gamma \f$.
 *
 *  In addition to the \a update vector, the recommended time step size \a dt is calculated
 *  employing a CFL condition.
 *
 *  \param t Current simulation time \f$\mathrm{[s]}\f$
 *  \param[out] dt Time step size \f$\mathrm{[s]}\f$
 *  \param[out] updateVec Update vector, or update estimate for secants, resp. Here in \f$\mathrm{[kg/m^3]}\f$
 *  \param impet Flag that determines if it is a real impet step or an update estimate for volume derivatives
 */
template<class TypeTag>
void FVTransport2P2CMultiPhysics<TypeTag>::update(const Scalar t, Scalar& dt, TransportSolutionType& updateVec, bool impet)
{
    // initialize dt very large
    dt = 1E100;

    unsigned int size = problem().gridView().size(0);
    if (this->enableLocalTimeStepping())
    {
        if (this->timeStepData_.size() != size)
            this->timeStepData_.resize(size);
    }
    // store if we do update Estimate for flux functions
    this->impet_ = impet;
    this->averagedFaces_ = 0.;

    // resize update vector and set to zero
    updateVec.resize(getPropValue<TypeTag, Properties::NumComponents>());
    updateVec[wCompIdx].resize(problem().gridView().size(0));
    updateVec[nCompIdx].resize(problem().gridView().size(0));
    updateVec[wCompIdx] = 0;
    updateVec[nCompIdx] = 0;

    // Cell which restricts time step size
    int restrictingCell = -1;

    PhaseVector entries(0.), timestepFlux(0.);
    // compute update vector
    for (const auto& element : elements(problem().gridView()))
    {
        // get cell infos
        int globalIdxI = problem().variables().index(element);
        CellData& cellDataI = problem().variables().cellData(globalIdxI);

        if (impet || cellDataI.subdomain() == 2)   // estimate only necessary in subdomain
        {
            // some variables for time step calculation
            double sumfactorin = 0;
            double sumfactorout = 0;

            // run through all intersections with neighbors and boundary
            for (const auto& intersection : intersections(problem().gridView(), element))
            {
                int indexInInside = intersection.indexInInside();

                /****** interior face   *****************/
                if (intersection.neighbor())
                    this->getFlux(entries, timestepFlux, intersection, cellDataI);

                /******  Boundary Face   *****************/
                if (intersection.boundary())
                    this->getFluxOnBoundary(entries, timestepFlux, intersection, cellDataI);

                if (this->enableLocalTimeStepping())
                {
                    LocalTimesteppingData& localData = this->timeStepData_[globalIdxI];

                    if (localData.faceTargetDt[indexInInside] < this->accumulatedDt_ + this->dtThreshold_)
                    {
                        localData.faceFluxes[indexInInside] = entries;
                    }
                }
                else
                {
                // add to update vector
                updateVec[wCompIdx][globalIdxI] += entries[wCompIdx];
                updateVec[nCompIdx][globalIdxI] += entries[nCompIdx];
                }

                // for time step calculation
                sumfactorin += timestepFlux[0];
                sumfactorout += timestepFlux[1];

            }// end all intersections

            if (this->enableLocalTimeStepping())
            {
                LocalTimesteppingData& localData = this->timeStepData_[globalIdxI];
                for (int i=0; i < 2*dim; i++)
                {
                    updateVec[wCompIdx][globalIdxI] += localData.faceFluxes[i][wCompIdx];
                    updateVec[nCompIdx][globalIdxI] += localData.faceFluxes[i][nCompIdx];
                }
            }

            /***********     Handle source term     ***************/
            PrimaryVariables q(NAN);
            problem().source(q, element);
            updateVec[wCompIdx][globalIdxI] += q[Indices::contiWEqIdx];
            updateVec[nCompIdx][globalIdxI] += q[Indices::contiNEqIdx];

            // account for porosity in fluxes for time-step
            using std::max;
            sumfactorin = max(sumfactorin,sumfactorout)
                            / problem().spatialParams().porosity(element);

            //calculate time step
            if (this->enableLocalTimeStepping())
            {
                this->timeStepData_[globalIdxI].dt = 1./sumfactorin;
                if ( 1./sumfactorin < dt)
                {
                    dt = 1./sumfactorin;
                    restrictingCell= globalIdxI;
                }
            }
            else
            {
                if ( 1./sumfactorin < dt)
                {
                    dt = 1./sumfactorin;
                    restrictingCell= globalIdxI;
                }
            }
        }
    } // end grid traversal

#if HAVE_MPI
    // communicate updated values
    using SolutionTypes = GetProp<TypeTag, Properties::SolutionTypes>;
    using ElementMapper = typename SolutionTypes::ElementMapper;
    using DataHandle = VectorCommDataHandleEqual<ElementMapper, Dune::BlockVector<Dune::FieldVector<Scalar, 1> >, 0/*elementCodim*/>;
    for (int i = 0; i < updateVec.size(); i++)
    {
        DataHandle dataHandle(problem().variables().elementMapper(), updateVec[i]);
        problem().gridView().template communicate<DataHandle>(dataHandle,
                                                            Dune::InteriorBorder_All_Interface,
                                                            Dune::ForwardCommunication);
    }
    dt = problem().gridView().comm().min(dt);
#endif

    if(impet)
    {
        Dune::dinfo << "Timestep restricted by CellIdx " << restrictingCell <<
          " leads to dt = "<<dt * getParam<Scalar>("Impet.CFLFactor")<< std::endl;
        if(this->averagedFaces_ != 0)
            Dune::dinfo  << " Averageing done for " << this->averagedFaces_ << " faces. "<< std::endl;
    }
    return;
}
} // end namespace Dumux
#endif
