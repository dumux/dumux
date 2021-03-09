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
 * \brief Finite volume 2p2c pressure model with multi-physics.
 */
#ifndef DUMUX_FVPRESSURE2P2C_MULTIPHYSICS_HH
#define DUMUX_FVPRESSURE2P2C_MULTIPHYSICS_HH

#include <dune/common/float_cmp.hh>

// dumux environment
#include <dumux/porousmediumflow/2p2c/sequential/fvpressure.hh>
#include <dumux/parallel/vectorcommdatahandle.hh>
#include <dumux/material/constraintsolvers/compositionalflash.hh>

#include <dumux/common/deprecated.hh>

namespace Dumux {
/*!
 * \ingroup SequentialTwoPTwoCModel
 * \brief The finite volume model for the solution of the compositional pressure equation
 *
 * Provides a Finite Volume implementation for the pressure equation of a gas-liquid
 * system with two components. An IMPES-like method is used for the sequential
 * solution of the problem.  Diffusion is neglected, capillarity can be regarded.
 * Isothermal conditions and local thermodynamic
 * equilibrium are assumed.  Gravity is included.
 * \f[
        c_{total}\frac{\partial p}{\partial t} + \sum_{\kappa} \frac{\partial v_{total}}{\partial C^{\kappa}}
        \nabla \cdot \left( \sum_{\alpha} X^{\kappa}_{\alpha} \varrho_{alpha} \bf{v}_{\alpha}\right)
         = \sum_{\kappa} \frac{\partial v_{total}}{\partial C^{\kappa}} q^{\kappa},
 * \f]
 * where \f$\bf{v}_{\alpha} = - \lambda_{\alpha} \bf{K} \left(\nabla p_{\alpha} + \rho_{\alpha} \bf{g} \right) \f$.
 * \f$ c_{total} \f$ represents the total compressibility, for constant porosity this yields
 * \f$ - \frac{\partial V_{total}}{\partial p_{\alpha}} \f$,
 * \f$p_{\alpha} \f$ denotes the phase pressure, \f$ \bf{K} \f$ the absolute permeability,
 * \f$ \lambda_{\alpha} \f$ the phase mobility,
 * \f$ \rho_{\alpha} \f$ the phase density and \f$ \bf{g} \f$ the gravity constant and
 * \f$ C^{\kappa} \f$ the total Component concentration.
 * See paper SPE 99619 or "Analysis of a Compositional Model for Fluid
 * Flow in Porous Media" by Chen, Qin and Ewing for derivation.
 *
 * The partial derivatives of the actual fluid volume \f$ v_{total} \f$ are gained by using a secant method.
 *
 * The model domain is automatically divided
 * in a single-phase and a two-phase domain. The full 2p2c model is only evaluated within the
 * two-phase subdomain, whereas a single-phase transport model is computed in the rest of the
 * domain.
 *
 * \tparam TypeTag The Type Tag
 */
template<class TypeTag>
class FVPressure2P2CMultiPhysics : public FVPressure2P2C<TypeTag>
{
    using ParentType = FVPressure2P2C<TypeTag>;
    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;
    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using SolutionTypes = GetProp<TypeTag, Properties::SolutionTypes>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;

    using SpatialParams = GetPropType<TypeTag, Properties::SpatialParams>;

    using Indices = typename GetPropType<TypeTag, Properties::ModelTraits>::Indices;
    using BoundaryTypes = GetPropType<TypeTag, Properties::SequentialBoundaryTypes>;

    using FluidSystem = GetPropType<TypeTag, Properties::FluidSystem>;
    using FluidState = GetPropType<TypeTag, Properties::FluidState>;

    using CellData = GetPropType<TypeTag, Properties::CellData>;
    enum
    {
        dim = GridView::dimension, dimWorld = GridView::dimensionworld
    };
    enum
    {
        pw = Indices::pressureW
    };
    enum
    {
        wPhaseIdx = Indices::wPhaseIdx, nPhaseIdx = Indices::nPhaseIdx,
        wCompIdx = Indices::wPhaseIdx, nCompIdx = Indices::nPhaseIdx,
        contiWEqIdx = Indices::contiWEqIdx, contiNEqIdx = Indices::contiNEqIdx
    };

    // using declarations to abbreviate several dune classes...
    using Element = typename GridView::Traits::template Codim<0>::Entity;
    using Grid = typename GridView::Grid;
    using Intersection = typename GridView::Intersection;

    // convenience shortcuts for Vectors/Matrices
    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;
    using DimMatrix = Dune::FieldMatrix<Scalar, dim, dim>;
    using PhaseVector = Dune::FieldVector<Scalar, getPropValue<TypeTag, Properties::NumPhases>()>;
    using PrimaryVariables = GetPropType<TypeTag, Properties::PrimaryVariables>;

    //! @copydoc FVPressure::EntryType
    using EntryType = Dune::FieldVector<Scalar, 2>;
//! Access functions to the current problem object
    Problem& problem()
    {    return this->problem_;   }
    const Problem& problem() const
    {    return this->problem_;   }
public:
    //function which assembles the system of equations to be solved
    void assemble(bool first);

    void get1pSource(EntryType& sourceEntry, const Element& elementI, const CellData& cellDataI);

    void get1pStorage(EntryType& storageEntry, const Element& elementI, CellData& cellDataI);

    void get1pFlux(EntryType& entries, const Intersection& intersection,
                   const CellData& cellDataI);

    void get1pFluxOnBoundary(EntryType& entries,
                             const Intersection& intersection,
                             const CellData& cellDataI);

    //initialize multi-physics-specific pressure model constituents
    void initialize(bool solveTwice = false)
    {
        // assign whole domain to most complex subdomain, => 2p
        int size = this->problem().gridView().size(0);
        for (int i = 0; i < size; i++)
        {
            CellData& cellData = this->problem().variables().cellData(i);
            cellData.subdomain() = 2;
        }
        nextSubdomain.resize(size);
        ParentType::initialize();
    }

    /*!
     * \brief  Function for serialization of the pressure field.
     *
     *  Function needed for restart option. Writes the pressure of a grid element to a restart file.
     *
     *  \param outstream Stream into the restart file.
     *  \param element Grid element
     */
    void serializeEntity(std::ostream &outstream, const Element &element)
    {
        ParentType::serializeEntity(outstream,element);
        int eIdxGlobal = problem().variables().index(element);
        CellData& cellData = problem().variables().cellData(eIdxGlobal);
        outstream <<"  "<< cellData.subdomain();
    }

    /*!
     * \brief  Function for deserialization of the pressure field.
     *
     *  Function needed for restart option. Reads the pressure of a grid element from a restart file.
     *
     *  \param instream Stream from the restart file.
     *  \param element Grid element
     */
    void deserializeEntity(std::istream &instream, const Element &element)
    {
        ParentType::deserializeEntity(instream,element);

        int eIdxGlobal = problem().variables().index(element);
        CellData& cellData = problem().variables().cellData(eIdxGlobal);
        int subdomainIdx;
        instream >> subdomainIdx;
        cellData.setSubdomainAndFluidStateType(subdomainIdx);
    }


    //constitutive functions are initialized and stored in the variables object
    void updateMaterialLaws(bool postTimeStep = false);
    //updates singlephase secondary variables for one cell and stores in the variables object
    void update1pMaterialLawsInElement(const Element& elementI, CellData& cellData, bool postTimeStep);

    /*!
     * \brief Write data files
     * \param writer The writer
     */
    template<class MultiWriter>
    void addOutputVtkFields(MultiWriter &writer)
    {
        FVPressureCompositional<TypeTag>::addOutputVtkFields(writer);

        if(problem().vtkOutputLevel()>=1)
        {
            int size = problem().gridView().size(0);
            // add multiphysics stuff
            Dune::BlockVector<Dune::FieldVector<int,1> >* subdomainPtr = writer.template allocateManagedBuffer<int, 1> (size);
            for (int i = 0; i < size; i++)
            {
                CellData& cellData = problem().variables().cellData(i);
                (*subdomainPtr)[i] = cellData.subdomain();
            }
            writer.attachCellData(*subdomainPtr, "subdomain");
        }

        return;
    }

    /*!
     * \brief Constructs a FVPressure2P2CPC object
     * \param problem a problem class object
     */
    FVPressure2P2CMultiPhysics(Problem& problem) : FVPressure2P2C<TypeTag>(problem),
            gravity_(problem.gravity()), timer_(false)
    {}

protected:
    #if HAVE_MPI
        using ElementMapper = typename SolutionTypes::ElementMapper;
        using DataHandle = VectorCommDataHandleEqual<ElementMapper, Dune::BlockVector<Dune::FieldVector<int, 1> >, 0/*elementCodim*/>;
    #endif

    // subdomain map
    Dune::BlockVector<Dune::FieldVector<int,1> > nextSubdomain;  //!< vector holding next subdomain
    const GlobalPosition& gravity_; //!< vector including the gravity constant
    //! gives kind of pressure used (\f$ 0 = p_w \f$, \f$ 1 = p_n \f$, \f$ 2 = p_{global} \f$)
    static constexpr int pressureType = getPropValue<TypeTag, Properties::PressureFormulation>();
    Dune::Timer timer_; //!< A timer for the time spent on the multiphysics framework.

    /*!
     * \brief Indices of matrix and rhs entries
     *
     * During the assembling of the global system of equations get-functions are called (getSource(),
     * getFlux(), etc.), which return global matrix or right hand side entries in a vector.
     * These can be accessed using following indices:
     */
    enum
    {
        rhs = 1,//!<index for the right hand side entry
        matrix = 0//!<index for the global matrix entry

    };
};

/*!
 * \brief  function which assembles the system of equations to be solved
 *
 * for first == true, this function assembles the matrix and right hand side for
 * the solution of the pressure field in the same way as in the class FVPressure2P.
 * for first == false, the approach is changed to \f[-\frac{\partial V}{\partial p}
 * \frac{\partial p}{\partial t}+\sum_{\kappa}\frac{\partial V}{\partial m^{\kappa}}\nabla\cdot
 * \left(\sum_{\alpha}C_{\alpha}^{\kappa}\mathbf{v}_{\alpha}\right)
 * =\sum_{\kappa}\frac{\partial V}{\partial m^{\kappa}}q^{\kappa} \f]. See Paper SPE 99619.
 * This is done to account for the volume effects which appear when gas and liquid are dissolved in each other.
 * \param first Flag if pressure field is unknown
 */
template<class TypeTag>
void FVPressure2P2CMultiPhysics<TypeTag>::assemble(bool first)
{
    if(first)
    {
        ParentType::assemble(true);
        return;
    }
    // initialization: set matrix this->A_ to zero
    this->A_ = 0;
    this->f_ = 0;

    for (const auto& element : elements(problem().gridView()))
    {
        // get the global index of the cell
        int eIdxGlobalI = problem().variables().index(element);

        // assemble interior element contributions
        if (element.partitionType() == Dune::InteriorEntity)
        {
            // get the cell data
            CellData& cellDataI = problem().variables().cellData(eIdxGlobalI);

            Dune::FieldVector<Scalar, 2> entries(0.);

            /*****  source term ***********/
            if(cellDataI.subdomain() != 2)
                problem().pressureModel().get1pSource(entries,element, cellDataI);
            else
                problem().pressureModel().getSource(entries,element, cellDataI, first);

            this->f_[eIdxGlobalI] = entries[rhs];

            /*****  flux term ***********/
            // iterate over all faces of the cell
            for (const auto& intersection : intersections(problem().gridView(), element))
            {
                /************* handle interior face *****************/
                if (intersection.neighbor())
                {
                    int eIdxGlobalJ = problem().variables().index(intersection.outside());

                    if (cellDataI.subdomain() != 2
                            or problem().variables().cellData(eIdxGlobalJ).subdomain() != 2) // cell in the 1p domain
                        get1pFlux(entries, intersection, cellDataI);
                    else
                        problem().pressureModel().getFlux(entries, intersection, cellDataI, first);

                    //set right hand side
                    this->f_[eIdxGlobalI] -= entries[rhs];
                    // set diagonal entry
                    this->A_[eIdxGlobalI][eIdxGlobalI] += entries[matrix];
                    // set off-diagonal entry
                    this->A_[eIdxGlobalI][eIdxGlobalJ] = -entries[matrix];
                }   // end neighbor


                /************* boundary face ************************/
                else
                {
                    if (cellDataI.subdomain() != 2) //the current cell in the 1p domain
                        problem().pressureModel().get1pFluxOnBoundary(entries, intersection, cellDataI);
                    else
                        problem().pressureModel().getFluxOnBoundary(entries, intersection, cellDataI, first);

                    //set right hand side
                    this->f_[eIdxGlobalI] += entries[rhs];
                    // set diagonal entry
                    this->A_[eIdxGlobalI][eIdxGlobalI] += entries[matrix];
                }
            } //end interfaces loop
    //        printmatrix(std::cout, this->A_, "global stiffness matrix", "row", 11, 3);

            /*****  storage term ***********/
            if (cellDataI.subdomain() != 2) //the current cell in the 1p domain
                problem().pressureModel().get1pStorage(entries, element, cellDataI);
            else
                problem().pressureModel().getStorage(entries, element, cellDataI, first);

            this->f_[eIdxGlobalI] += entries[rhs];
            // set diagonal entry
            this->A_[eIdxGlobalI][eIdxGlobalI] += entries[matrix];
        }
        // assemble overlap and ghost element contributions
        else
        {
            this->A_[eIdxGlobalI] = 0.0;
            this->A_[eIdxGlobalI][eIdxGlobalI] = 1.0;
            this->f_[eIdxGlobalI] = this->pressure()[eIdxGlobalI];
        }
    } // end grid traversal
//        printmatrix(std::cout, this->A_, "global stiffness matrix after assempling", "row", 11,3);
//        printvector(std::cout, this->f_, "right hand side", "row", 10);
    return;
}

/*!
 * \brief Assembles the source term
 *
 * The source is translated into a volumentric source term:
 * \f[ V_i \sum_{\kappa} \frac{1}{\varrho} q^{\kappa}_i \; , \f]
 * because under singlephase conditions
 * \f[ \frac{\partial v_{t}}{\partial C^{\kappa}} \approx \frac{1}{\varrho} \f].
 * \param sourceEntry The Matrix and RHS entries
 * \param elementI The element I
 * \param cellDataI Data of cell I
 */
template<class TypeTag>
void FVPressure2P2CMultiPhysics<TypeTag>::get1pSource(Dune::FieldVector<Scalar, 2>& sourceEntry,
        const Element& elementI, const CellData& cellDataI)
{
    sourceEntry=0.;

    // cell volume & perimeter, assume linear map here
    Scalar volume = elementI.geometry().volume();
    int subdomainIdx = cellDataI.subdomain();

    /****************implement source************************/
    PrimaryVariables source(NAN);
    problem().source(source, elementI);
    source[1+subdomainIdx] /= cellDataI.density(subdomainIdx);

    sourceEntry[1] = volume * source[1+subdomainIdx];

    return;
}

/*!
 * \brief  Assembles the storage term for a 1p cell in a multiphysics framework
 *
 * The storage term comprises the (single-phase) compressibility (due to a change in
 * pressure from last timestep):
 *  \f[ V_i c_{i} \frac{p^t_i - p^{t-\Delta t}_i}{\Delta t} \f]
 * and the damped error introduced by the incorrect transport of the last timestep:
 *  \f[ V_i \alpha_r \frac{v_{t} - \phi}{\Delta t} \f].
 * The latter is damped according to Fritz 2011.
 * \param storageEntry The Matrix and RHS entries
 * \param elementI The element I
 * \param cellDataI Data of cell I
 */
template<class TypeTag>
void FVPressure2P2CMultiPhysics<TypeTag>::get1pStorage(Dune::FieldVector<Scalar, 2>& storageEntry,
                                                        const Element& elementI,
                                                        CellData& cellDataI)
{
    storageEntry = 0.;
    // cell index
    int eIdxGlobalI = problem().variables().index(elementI);
    int presentPhaseIdx = cellDataI.subdomain();
    Scalar volume = elementI.geometry().volume();

    // determine maximum error to scale error-term
    Scalar timestep_ = problem().timeManager().timeStepSize();

    // compressibility term: 1p domain, so no dv_dp calculated
    if (true)
    {
        Scalar& incp = this->incp_;

        // numerical derivative of fluid volume with respect to pressure
        PhaseVector p_(incp);
        p_[nPhaseIdx] += cellDataI.pressure(nPhaseIdx);
        p_[wPhaseIdx] += cellDataI.pressure(wPhaseIdx);

        Scalar sumC = (cellDataI.massConcentration(wCompIdx) + cellDataI.massConcentration(nCompIdx));
        Scalar Z0 = cellDataI.massConcentration(wCompIdx) / sumC;
        // initialize simple fluidstate object
        PseudoOnePTwoCFluidState<Scalar, FluidSystem> pseudoFluidState;
        CompositionalFlash<Scalar, FluidSystem> flashSolver;
        flashSolver.concentrationFlash1p2c(pseudoFluidState, Z0, p_, cellDataI.subdomain(),
                cellDataI.temperature(wPhaseIdx));

        Scalar v_ = 1. / pseudoFluidState.density(presentPhaseIdx);
        cellDataI.dv_dp() = (sumC * ( v_ - (1. /cellDataI.density(presentPhaseIdx)))) /incp;

        if (cellDataI.dv_dp()>0)
        {
            // dV_dp > 0 is unphysical: Try inverse increment for secant
            Dune::dinfo << "dv_dp larger 0 at Idx " << eIdxGlobalI << " , try and invert secant"<< std::endl;

            p_ -= 2*incp;
            flashSolver.concentrationFlash1p2c(pseudoFluidState, Z0, p_, cellDataI.subdomain(),
                    cellDataI.temperature(wPhaseIdx));
            v_ = 1. / pseudoFluidState.density(presentPhaseIdx);
            cellDataI.dv_dp() = (sumC * ( v_ - (1. /cellDataI.density(presentPhaseIdx)))) /incp;
            // dV_dp > 0 is unphysical: Try inverse increment for secant
            if (cellDataI.dv_dp()>0)
            {
                Dune::dinfo <<__FILE__<< "dv_dp still larger 0 after inverting secant. regularize"<< std::endl;
                cellDataI.dv_dp() *= -1;
            }
        }


        Scalar compress_term = cellDataI.dv_dp() / timestep_;

        storageEntry[matrix] -= compress_term*volume;
        storageEntry[rhs] -= cellDataI.pressure(pressureType) * compress_term * volume;

        using std::isnan;
        using std::isinf;
        if (isnan(compress_term) ||isinf(compress_term))
            DUNE_THROW(Dune::MathError, "Compressibility term leads to NAN matrix entry at index " << eIdxGlobalI);

        if(!getPropValue<TypeTag, Properties::EnableCompressibility>())
            DUNE_THROW(Dune::NotImplemented, "Compressibility is switched off???");
    }

    // Abort error damping if there will be a possibly tiny timestep compared with last one
    // This might be the case if the episode or simulation comes to an end.
    if( problem().timeManager().episodeWillBeFinished()
            || problem().timeManager().willBeFinished())
    {
        problem().variables().cellData(eIdxGlobalI).errorCorrection() = 0.;
        return;
    }

    // error reduction routine: volumetric error is damped and inserted to right hand side
    // if damping is not done, the solution method gets unstable!
    problem().variables().cellData(eIdxGlobalI).volumeError() /= timestep_;
    Scalar maxError = this->maxError_;
    Scalar erri = fabs(cellDataI.volumeError());
    Scalar x_lo = this->ErrorTermLowerBound_;
    Scalar x_mi = this->ErrorTermUpperBound_;
    Scalar fac  = this->ErrorTermFactor_;
    if (pressureType == pw)
        fac = 0.1*this->ErrorTermFactor_;
    Scalar lofac = 0.;
    Scalar hifac = 0.;

    if ((erri*timestep_ > 5e-5) && (erri > x_lo * maxError) && (!problem().timeManager().willBeFinished()))
    {
        if (erri <= x_mi * maxError)
            storageEntry[rhs] +=
                    problem().variables().cellData(eIdxGlobalI).errorCorrection() =
                            fac* (1-x_mi*(lofac-1)/(x_lo-x_mi) + (lofac-1)/(x_lo-x_mi)*erri/maxError)
                                * cellDataI.volumeError() * volume;
        else
            storageEntry[rhs] +=
                    problem().variables().cellData(eIdxGlobalI).errorCorrection() =
                            fac * (1 + x_mi - hifac*x_mi/(1-x_mi) + (hifac/(1-x_mi)-1)*erri/maxError)
                                * cellDataI.volumeError() * volume;
    }
    else
        problem().variables().cellData(eIdxGlobalI).errorCorrection()=0 ;

    return;
}

/*!
 * \brief The compositional single-phase flux in the multiphysics framework
 *
 * If only single-phase conditions are encountered, the flux expression simplifies to (written for the
 * case where the wetting phase is only present):
   \f[
          A_{\gamma} \mathbf{n}_{\gamma}^T \mathbf{K}
       \lambda_w \mathbf{d}_{ij} \left( \frac{p_{w,j}^t - p^{t}_{w,i}}{\Delta x} + \varrho_{w} \mathbf{g}^T \mathbf{d}_{ij} \right) .
   \f]
 *
 * \param entries The Matrix and RHS entries
 * \param intersection Intersection between cell I and J
 * \param cellDataI Data of cell I
 */
template<class TypeTag>
void FVPressure2P2CMultiPhysics<TypeTag>::get1pFlux(Dune::FieldVector<Scalar, 2>& entries,
        const Intersection& intersection, const CellData& cellDataI)
{
    entries = 0.;
    auto elementI = intersection.inside();

    // get global coordinate of cell center
    const GlobalPosition& globalPos = elementI.geometry().center();

    // get absolute permeability
    DimMatrix permeabilityI(problem().spatialParams().intrinsicPermeability(elementI));

    // get normal vector
    const GlobalPosition& unitOuterNormal = intersection.centerUnitOuterNormal();

    // get face volume
    Scalar faceArea = intersection.geometry().volume();

        // access neighbor
        auto neighbor = intersection.outside();
        int eIdxGlobalJ = problem().variables().index(neighbor);
        CellData& cellDataJ = problem().variables().cellData(eIdxGlobalJ);

        // gemotry info of neighbor
        const GlobalPosition& globalPosNeighbor = neighbor.geometry().center();

        // distance vector between barycenters
        GlobalPosition distVec = globalPosNeighbor - globalPos;

        // compute distance between cell centers
        Scalar dist = distVec.two_norm();

        GlobalPosition unitDistVec(distVec);
        unitDistVec /= dist;

        DimMatrix permeabilityJ
            = problem().spatialParams().intrinsicPermeability(neighbor);

        // compute vectorized permeabilities
        DimMatrix meanPermeability(0);
        harmonicMeanMatrix(meanPermeability, permeabilityI, permeabilityJ);

        Dune::FieldVector<Scalar, dim> permeability(0);
        meanPermeability.mv(unitDistVec, permeability);

        // due to "safety cell" around subdomain, both cells I and J
        // have single-phase conditions, although one is in 2p domain.
        using std::min;
        int phaseIdx = min(cellDataI.subdomain(), cellDataJ.subdomain());

        Scalar rhoMean = 0.5 * (cellDataI.density(phaseIdx) + cellDataJ.density(phaseIdx));

        // 1p => no pc => only 1 pressure, potential
        Scalar potential = (cellDataI.pressure(phaseIdx) - cellDataJ.pressure(phaseIdx)) / dist;

        potential += rhoMean * (unitDistVec * gravity_);

        Scalar lambda;

        if (potential > 0.)
        {
            lambda = cellDataI.mobility(phaseIdx);
            cellDataJ.setUpwindCell(intersection.indexInOutside(), contiWEqIdx, false);  // store in cellJ since cellI is const
            cellDataJ.setUpwindCell(intersection.indexInOutside(), contiNEqIdx, false);  // store in cellJ since cellI is const
        }
        else if (potential < 0.)
        {
            lambda = cellDataJ.mobility(phaseIdx);
            cellDataJ.setUpwindCell(intersection.indexInOutside(), contiWEqIdx, true);
            cellDataJ.setUpwindCell(intersection.indexInOutside(), contiNEqIdx, true);
        }
        else
        {
            lambda = harmonicMean(cellDataI.mobility(phaseIdx) , cellDataJ.mobility(phaseIdx));
            cellDataJ.setUpwindCell(intersection.indexInOutside(), contiWEqIdx, false);
            cellDataJ.setUpwindCell(intersection.indexInOutside(), contiNEqIdx, false);
        }

        entries[0]  = lambda * faceArea * fabs(permeability * unitOuterNormal) / (dist);
        entries[1]  = rhoMean * lambda;
        entries[1] *= faceArea * fabs(permeability * unitOuterNormal) * (unitDistVec * gravity_);

        return;
}

/*!
 * \brief The compositional single-phase flux in the multiphysics framework
 *
 * If only single-phase conditions are encountered, the flux expression simplifies to (written for the
 * case where the wetting phase is only present):
   \f[
          A_{\gamma} \mathbf{n}_{\gamma}^T \mathbf{K}
      \varrho_w \lambda_w \mathbf{d}_{i-Boundary} \left( \frac{p_{w,Boundary}^t - p^{t}_{w,i}}{\Delta x}
      + \varrho_{w} \mathbf{g}^T \mathbf{d}_{i-Boundary} \right) .
   \f]
 *
 * If a Neumann BC is set, the given (mass-)flux is directly multiplied by the volume derivative and inserted.
 *
 * \param entries The Matrix and RHS entries
 * \param intersection Intersection between cell I and J
 * \param cellDataI Data of cell I
 */
template<class TypeTag>
void FVPressure2P2CMultiPhysics<TypeTag>::get1pFluxOnBoundary(Dune::FieldVector<Scalar, 2>& entries,
        const Intersection& intersection, const CellData& cellDataI)
{
    entries = 0.;
    // get global coordinate of cell center
    auto elementI = intersection.inside();
    const GlobalPosition& globalPos = elementI.geometry().center();
//    int eIdxGlobalI = problem().variables().index(elementI);
    int phaseIdx = cellDataI.subdomain();

    // get normal vector
    const GlobalPosition& unitOuterNormal = intersection.centerUnitOuterNormal();
    // get face volume
    Scalar faceArea = intersection.geometry().volume();

                // center of face in global coordinates
                const GlobalPosition& globalPosFace = intersection.geometry().center();

                // geometrical information
                GlobalPosition distVec(globalPosFace - globalPos);
                Scalar dist = distVec.two_norm();
                GlobalPosition unitDistVec(distVec);
                unitDistVec /= dist;

                //get boundary condition for boundary face center
                BoundaryTypes bcType;
                problem().boundaryTypes(bcType, intersection);

                // prepare pressure boundary condition
                PhaseVector pressBC(0.);

                /**********         Dirichlet Boundary        *************/
                if (bcType.isDirichlet(Indices::pressureEqIdx))
                {
                    // get absolute permeability
                    DimMatrix permeabilityI(problem().spatialParams().intrinsicPermeability(elementI));
                    if(this->regulateBoundaryPermeability)
                    {
                        int axis = intersection.indexInInside() / 2;
                        if(permeabilityI[axis][axis] < this->minimalBoundaryPermeability)
                            permeabilityI[axis][axis] = this->minimalBoundaryPermeability;
                    }
                    // get mobilities and fractional flow factors
                    Scalar lambdaI = cellDataI.mobility(phaseIdx);

                    //permeability vector at boundary
                    Dune::FieldVector<Scalar, dim> permeability(0);
                    permeabilityI.mv(unitDistVec, permeability);

                    // create a fluid state for the boundary
                    FluidState BCfluidState;

                    //read boundary values
                    PrimaryVariables primaryVariablesOnBoundary(NAN);
                    problem().dirichlet(primaryVariablesOnBoundary, intersection);

                    {
                        // read boundary values
                        problem().transportModel().evalBoundary(globalPosFace,
                                                                    intersection,
                                                                    BCfluidState,
                                                                    pressBC);

                        // determine fluid properties at the boundary
                        Scalar densityBound =
                            FluidSystem::density(BCfluidState, phaseIdx);
                        Scalar viscosityBound =
                            FluidSystem::viscosity(BCfluidState, phaseIdx);

                        // mobility at the boundary
                        Scalar lambdaBound = 0.;
                        switch (getPropValue<TypeTag, Properties::BoundaryMobility>())
                        {
                            case Indices::satDependent:
                            {
                                lambdaBound = BCfluidState.saturation(phaseIdx) / viscosityBound;
                                break;
                            }
                            case Indices::permDependent:
                            {
                                // old material law interface is deprecated: Replace this by
                                // const auto& fluidMatrixInteraction = problem().spatialParams.fluidMatrixInteractionAtPos(elementI.geometry().center());
                                // after the release of 3.3, when the deprecated interface is no longer supported
                                const auto fluidMatrixInteraction = Deprecated::makePcKrSw(Scalar{}, problem().spatialParams(), elementI);

                                if (phaseIdx == wPhaseIdx)
                                    lambdaBound = fluidMatrixInteraction.krw(BCfluidState.saturation(wPhaseIdx)) / viscosityBound;
                                else
                                    lambdaBound = fluidMatrixInteraction.krn(BCfluidState.saturation(wPhaseIdx)) / viscosityBound;
                                break;
                            }
                        }
                        Scalar rhoMean = 0.5 * (cellDataI.density(phaseIdx) + densityBound);

                        Scalar potential = 0;

                        //calculate potential gradient: pc = 0;
                        potential = (cellDataI.pressure(phaseIdx) - pressBC[phaseIdx]) / dist;

                        potential += rhoMean * (unitDistVec * gravity_);

                        Scalar lambda(0.);

                        if (Dune::FloatCmp::eq<Scalar, Dune::FloatCmp::absolute>(potential, 0.0, 1.0e-30))
                        {
                            lambda = 0.5*(lambdaI + lambdaBound);
                        }
                        else if (potential > 0.)
                        {
                            lambda = lambdaI;
                        }
                        else
                        {
                            lambda = lambdaBound;
                        }

                        //calculate current matrix entry
                        Scalar entry(0.), rightEntry(0.);
                        entry = lambda * (fabs(permeability * unitOuterNormal) / dist) * faceArea;

                        //calculate right hand side
                        rightEntry = lambda * rhoMean * fabs(permeability * unitOuterNormal)
                                * faceArea ;

                        // set diagonal entry and right hand side entry
                        entries[0] += entry;
                        entries[1] += entry * primaryVariablesOnBoundary[Indices::pressureEqIdx];
                        entries[1] -= rightEntry * (gravity_ * unitDistVec);
                    }    //end of if(first) ... else{...
                }   // end dirichlet

                /**********************************
                 * set neumann boundary condition
                 **********************************/
                else if(bcType.isNeumann(Indices::pressureEqIdx))
                {
                    PrimaryVariables J(NAN);
                    problem().neumann(J, intersection);
                    J[1+phaseIdx] /= cellDataI.density(phaseIdx);

                    entries[1] -= J[1+phaseIdx] * faceArea;
                }
                else
                    DUNE_THROW(Dune::NotImplemented, "Boundary Condition neither Dirichlet nor Neumann!");


    return;
}


/*!
 * \brief constitutive functions are updated once if new concentrations are calculated and stored in the variables container
 *
 * In contrast to the standard sequential 2p2c model ( FVPressure2P2C<TypeTag>::updateMaterialLaws() ),
 * this method also holds routines to adapt the subdomain. The subdomain indicates weather we are in 1p domain (value = 1)
 * or in the two phase subdomain (value = 2).
 * Note that the type of flash, i.e. the type of FluidState (FS), present in each cell does not have to
 * coincide with the subdomain. If a cell will be simple and was complex, a complex FS is available, so next time step
 * will use this complex FS, but updateMaterialLaw afterwards will finally transform that to simple FS.
 * \param postTimeStep Flag indicating method is called from Problem::postTimeStep()
 */
template<class TypeTag>
void FVPressure2P2CMultiPhysics<TypeTag>::updateMaterialLaws(bool postTimeStep)
{
    //get timestep for error term
    Scalar maxError = 0.;

    // next subdomain map
    if (Dune::FloatCmp::eq<Scalar, Dune::FloatCmp::absolute>(problem().timeManager().time(), 0.0, 1.0e-30))
        nextSubdomain = 2;  // start with complicated sub in initialization
    else
        nextSubdomain = -1;  // reduce complexity after first TS

    // Loop A) through leaf grid
    for (const auto& element : elements(problem().gridView()))
    {
        // get global coordinate of cell center
        int eIdxGlobal = problem().variables().index(element);
        CellData& cellData = problem().variables().cellData(eIdxGlobal);

        if(cellData.subdomain() == 2)    // complex
        {
            this->updateMaterialLawsInElement(element, postTimeStep);

            // check subdomain consistency
            timer_.start();
            // enshure we are not at source
            // get sources from problem
            PrimaryVariables source(NAN);
            problem().source(source, element);

            if ((cellData.saturation(wPhaseIdx) > 0.0 && cellData.saturation(wPhaseIdx) < 1.0)
                || Dune::FloatCmp::ne<Scalar, Dune::FloatCmp::absolute>(source.one_norm(), 0.0, 1.0e-30)) // cell still 2p
            {
                // mark this element
                nextSubdomain[eIdxGlobal] = 2;

                // mark neighbors
                for (const auto& intersection : intersections(problem().gridView(), element))
                {
                    if (intersection.neighbor())
                    {
                        int eIdxGlobalJ = problem().variables().index(intersection.outside());
                        // mark neighbor Element
                        nextSubdomain[eIdxGlobalJ] = 2;
                    }
                }
            }
            else if(nextSubdomain[eIdxGlobal] != 2)// update next subdomain if possible
            {
                if(Dune::FloatCmp::ne<Scalar, Dune::FloatCmp::absolute>(cellData.saturation(wPhaseIdx), 0.0, 1.0e-30))
                    nextSubdomain[eIdxGlobal] = wPhaseIdx;
                else if (Dune::FloatCmp::ne<Scalar, Dune::FloatCmp::absolute>(cellData.saturation(nPhaseIdx), 0.0, 1.0e-30))
                    nextSubdomain[eIdxGlobal] = nPhaseIdx;
            }
            timer_.stop();
            // end subdomain check
        }// end complex domain
        else if (nextSubdomain[eIdxGlobal] != 2) //check if cell remains in simple subdomain
            nextSubdomain[eIdxGlobal] = cellData.subdomain();

    } //end define complex area of next subdomain

    timer_.start();
    //communicate next subdomain if parallel
    #if HAVE_MPI
    // communicate updated values
    DataHandle dataHandle(problem().variables().elementMapper(), nextSubdomain);
    problem().gridView().template communicate<DataHandle>(dataHandle,
                                                        Dune::InteriorBorder_All_Interface,
                                                        Dune::ForwardCommunication);
    #endif

    // Loop B) thorugh leaf grid
    // investigate cells that were "simple" in current TS
    for (const auto& element : elements(problem().gridView()))
    {
        int eIdxGlobal = problem().variables().index(element);
        CellData& cellData = problem().variables().cellData(eIdxGlobal);

        // store old subdomain information and assign new info
        int oldSubdomainI = cellData.subdomain();
        cellData.subdomain() = nextSubdomain[eIdxGlobal];

        //first check if simple will become complicated
        if(oldSubdomainI != 2
                    && nextSubdomain[eIdxGlobal] == 2)
        {
            // use complex update of the fluidstate
            timer_.stop();
            this->updateMaterialLawsInElement(element, postTimeStep);
            timer_.start();
        }
        else if(oldSubdomainI != 2
                    && nextSubdomain[eIdxGlobal] != 2)    // will be simple and was simple
        {
            // perform simple update
            this->update1pMaterialLawsInElement(element, cellData, postTimeStep);
        }
        //else
        // a) will remain complex -> everything already done in loop A
        // b) will be simple and was complex: complex FS available, so next TS
        //         will use comlex FS, next updateMaterialLaw will transform to simple

        using std::max;
        maxError = max(maxError, fabs(cellData.volumeError()));
    }// end grid traversal
    this->maxError_ = maxError/problem().timeManager().timeStepSize();

    timer_.stop();

    if(problem().timeManager().willBeFinished() or problem().timeManager().episodeWillBeFinished())
        Dune::dinfo << "Subdomain routines took " << timer_.elapsed() << " seconds" << std::endl;

    return;
}

/*!
 * \brief updates secondary variables of one single phase cell
 *
 * For each element, the secondary variables are updated according to the
 * primary variables. Only a simple flash calulation has to be carried out,
 * as phase distribution is already known: single-phase.
 * \param elementI The element
 * \param cellData The cell data of the current element
 * \param postTimeStep Flag indicating if we have just completed a time step
 */
template<class TypeTag>
void FVPressure2P2CMultiPhysics<TypeTag>::update1pMaterialLawsInElement(const Element& elementI, CellData& cellData, bool postTimeStep)
{
    // get global coordinate of cell center
    GlobalPosition globalPos = elementI.geometry().center();
    int eIdxGlobal = problem().variables().index(elementI);

    // determine which phase should be present
    int presentPhaseIdx = cellData.subdomain(); // this is already =nextSubomainIdx

    // reset to calculate new timeStep if we are at the end of a time step
    if(postTimeStep)
        cellData.reset();

    // acess the simple fluid state and prepare for manipulation
    auto& pseudoFluidState = cellData.manipulateSimpleFluidState();

    // old material law interface is deprecated: Replace this by
    // const auto& fluidMatrixInteraction = problem().spatialParams.fluidMatrixInteractionAtPos(elementI.geometry().center());
    // after the release of 3.3, when the deprecated interface is no longer supported
    const auto fluidMatrixInteraction = Deprecated::makePcKrSw(Scalar{}, problem().spatialParams(), elementI);

    // prepare phase pressure for fluid state
    // both phase pressures are necessary for the case 1p domain is assigned for
    // the next 2p subdomain
    PhaseVector pressure(0.);
    Scalar pc = 0;
    if(getPropValue<TypeTag, Properties::EnableCapillarity>())
        pc = fluidMatrixInteraction.pc(((presentPhaseIdx == wPhaseIdx) ? 1. : 0.)); // assign sw = 1 if wPhase present, else 0
    if(pressureType == wPhaseIdx)
    {
        pressure[wPhaseIdx] = this->pressure(eIdxGlobal);
        pressure[nPhaseIdx] = this->pressure(eIdxGlobal)+pc;
    }
    else
    {
        pressure[wPhaseIdx] = this->pressure(eIdxGlobal)-pc;
        pressure[nPhaseIdx] = this->pressure(eIdxGlobal);
    }

    // get the overall mass of first component:  Z0 = C^0 / (C^0+C^1) [-]
    Scalar sumConc = cellData.massConcentration(wCompIdx)
            + cellData.massConcentration(nCompIdx);
    Scalar Z0 = cellData.massConcentration(wCompIdx)/ sumConc;

    CompositionalFlash<Scalar, FluidSystem> flashSolver;
    flashSolver.concentrationFlash1p2c(pseudoFluidState, Z0, pressure, presentPhaseIdx, problem().temperatureAtPos(globalPos));

    // write stuff in fluidstate
    assert(presentPhaseIdx == pseudoFluidState.presentPhaseIdx());

//    cellData.setSimpleFluidState(pseudoFluidState);

    // initialize viscosities
    cellData.setViscosity(presentPhaseIdx, FluidSystem::viscosity(pseudoFluidState, presentPhaseIdx));

    // initialize mobilities
    if(presentPhaseIdx == wPhaseIdx)
    {
        cellData.setMobility(wPhaseIdx,
            fluidMatrixInteraction.krw(pseudoFluidState.saturation(wPhaseIdx)) / cellData.viscosity(wPhaseIdx));
        cellData.setMobility(nPhaseIdx, 0.);
    }
    else
    {
        cellData.setMobility(nPhaseIdx,
            fluidMatrixInteraction.krn(pseudoFluidState.saturation(wPhaseIdx)) / cellData.viscosity(nPhaseIdx));
        cellData.setMobility(wPhaseIdx, 0.);
    }

    // error term handling
    Scalar vol(0.);
    vol = sumConc / pseudoFluidState.density(presentPhaseIdx);

    if (Dune::FloatCmp::ne<Scalar, Dune::FloatCmp::absolute>(problem().timeManager().timeStepSize(), 0.0, 1.0e-30))
        cellData.volumeError() = (vol - problem().spatialParams().porosity(elementI));
    return;
}

}//end namespace Dumux
#endif
