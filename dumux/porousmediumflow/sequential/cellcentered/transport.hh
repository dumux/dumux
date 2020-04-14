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
#ifndef DUMUX_FVTRANSPORT_HH
#define DUMUX_FVTRANSPORT_HH

#include <dune/grid/common/gridenums.hh>
#include <dumux/porousmediumflow/sequential/transportproperties.hh>
#include <dumux/porousmediumflow/sequential/properties.hh>
#include <dumux/parallel/vectorcommdatahandle.hh>
#include <unordered_map>

/**
 * @file
 * @brief  Finite Volume discretization of a  transport equation
 */

namespace Dumux
{
//! \ingroup IMPET
/*!\brief The finite volume discretization of a transport equation
 *
 *  Base class for finite volume (FV) implementations of an explicitly treated transport equation.
 *  The class provides a method to calculate the explicit update to get a new solution of the transported quantity:
 *  \f[
 *      u_{new} = u_{old} + \Delta t \Delta u_{update}
 *  \f]
 *  A certain transport equation defined in a implementation of this base class must be splitted
 *  into a flux term and a source term.
 *  Corresponding functions (<tt>getSource()</tt>, <tt>getFlux()</tt> and <tt>getFluxOnBoundary()</tt>)
 *  have to be defined in the implementation.
 *
 * \tparam TypeTag The Type Tag
 */
template<class TypeTag>
class FVTransport
{
    using Implementation = GetPropType<TypeTag, Properties::TransportModel>;

    using GridView = typename GetPropType<TypeTag, Properties::GridGeometry>::GridView;

    enum
        {
            dim = GridView::dimension, dimWorld = GridView::dimensionworld
        };

    using Scalar = GetPropType<TypeTag, Properties::Scalar>;
    using Problem = GetPropType<TypeTag, Properties::Problem>;

    using TransportSolutionType = GetPropType<TypeTag, Properties::TransportSolutionType>;
    using CellData = GetPropType<TypeTag, Properties::CellData>;

    using EvalCflFluxFunction = GetPropType<TypeTag, Properties::EvalCflFluxFunction>;

    using Element = typename GridView::Traits::template Codim<0>::Entity;
    using Intersection = typename GridView::Intersection;

    using GlobalPosition = Dune::FieldVector<Scalar, dimWorld>;

    struct LocalTimesteppingData
    {
        Dune::FieldVector<Scalar, 2*dim> faceFluxes;
        Dune::FieldVector<Scalar, 2*dim> faceTargetDt;
        Scalar dt;
        LocalTimesteppingData():faceFluxes(0.0), faceTargetDt(0.0), dt(0)
        {}
    };

protected:
    //! \cond \private
    EvalCflFluxFunction& evalCflFluxFunction()
    {
        return *evalCflFluxFunction_;
    }

    const EvalCflFluxFunction& evalCflFluxFunction() const
    {
        return *evalCflFluxFunction_;
    }
    //! \endcond

    void innerUpdate(TransportSolutionType& updateVec);

public:

    // Calculate the update vector.
    void update(const Scalar t, Scalar& dt, TransportSolutionType& updateVec, bool impet = false);

    void updateTransport(const Scalar t, Scalar& dt, TransportSolutionType& updateVec)
    {
        asImp_().updateMaterialLaws();
        asImp_().update(t, dt, updateVec);
    }

    /*! \brief Function which calculates the flux update
     *
     * Function computes the inter-cell flux term and adds it to the update.
     *
     * \param update The cell update
     * \param intersection Intersection of two grid elements
     * \param cellDataI Object containing all model relevant cell data
     */
    void getFlux(Scalar& update, const Intersection& intersection, CellData& cellDataI);

    /*! \brief Function which calculates the boundary flux update
     *
     * Function computes the boundary-flux term and  adds it to the update.
     *
     * \param update The cell update
     * \param intersection Intersection of two grid elements
     * \param cellDataI Object containing all model relevant cell data
     */
    void getFluxOnBoundary(Scalar& update, const Intersection& intersection, CellData& cellDataI);

    /*! \brief Function which calculates the source update
     *
     * Function computes the source term and adds it to the update.
     *
     * \param update The cell update
     * \param element Grid element
     * \param cellDataI Object containing all model relevant cell data
     */
    void getSource(Scalar& update, const Element& element, CellData& cellDataI);

    //! Sets the initial solution \f$ S_0 \f$.
    void initialize()
    {
        evalCflFluxFunction_->initialize();
    }

    /*! \brief Updates constitutive relations and stores them in the variable class*/
    void updateMaterialLaws();

    /*! \brief Writes the current values of the primary transport variable into the
     *  <tt>transportedQuantity</tt>-vector (comes as function argument)
     *
     * \param transportedQuantity Vector of the size of global numbers of degrees of freedom
     *  of the primary transport variable.
     */
    void getTransportedQuantity(TransportSolutionType& transportedQuantity);

    template<class DataEntry>
    bool inPhysicalRange(DataEntry& entry);

    /*! \brief Writes the current values of the primary transport variable into the variable container
     *
     * \param transportedQuantity Vector of the size of global numbers of degrees of freedom of the primary transport variable.
     */
    void setTransportedQuantity(TransportSolutionType& transportedQuantity);

    /*! \brief Updates the primary transport variable.
     *
     * \param updateVec Vector containing the global update.
     */
    void updateTransportedQuantity(TransportSolutionType& updateVec);

    /*! \brief Adds transport output to the output file
     *
     * \tparam MultiWriter Class defining the output writer
     * \param writer The output writer (usually a <tt>VTKMultiWriter</tt> object)
     *
     */
    template<class MultiWriter>
    void addOutputVtkFields(MultiWriter &writer)
    {}

    /*! \brief  Function for serialization of the primary transport variable.
     *
     *  Function needed for restart option. Writes the primary transport variable of a grid element to a restart file.
     *
     *  \param outstream Stream into the restart file.
     *  \param element Grid element
     */
    void serializeEntity(std::ostream &outstream, const Element &element)
    {}

    /*! \brief  Function for deserialization of the primary transport variable.
     *
     *  Function needed for restart option. Reads the the primary transport variable of a grid element from a restart file.
     *
     *  \param instream Stream from the restart file.
     *  \param element Grid element
     */
    void deserializeEntity(std::istream &instream, const Element &element)
    {}

    bool enableLocalTimeStepping()
    {
        return localTimeStepping_;
    }


    //! Constructs a FVTransport object
    /**

     * \param problem A problem class object
     */
    FVTransport(Problem& problem) :
        problem_(problem), switchNormals_(getParam<bool>("Impet.SwitchNormals")),
        subCFLFactor_(1.0), accumulatedDt_(0), dtThreshold_(1e-6)
    {
        evalCflFluxFunction_ = std::make_shared<EvalCflFluxFunction>(problem);

        Scalar cFLFactor = getParam<Scalar>("Impet.CFLFactor");
        using std::min;
        subCFLFactor_ = min(getParam<Scalar>("Impet.SubCFLFactor"), cFLFactor);
        verbosity_ = getParam<int>("TimeManager.SubTimestepVerbosity");

        localTimeStepping_ = subCFLFactor_/cFLFactor < 1.0 - dtThreshold_;

        if (localTimeStepping_)
            std::cout<<"max CFL-Number of "<<cFLFactor<<", max Sub-CFL-Number of "<<subCFLFactor_<<": Enable local time-stepping!\n";
    }

private:
    //! Returns the implementation of the problem (i.e. static polymorphism)
    Implementation &asImp_()
    { return *static_cast<Implementation *>(this); }

    //! \copydoc IMPETProblem::asImp_()
    const Implementation &asImp_() const
    { return *static_cast<const Implementation *>(this); }

    void updatedTargetDt_(Scalar &dt);

    void resetTimeStepData_()
    {
        timeStepData_.clear();
        accumulatedDt_ = 0;
    }

    Problem& problem_;
    bool switchNormals_;

    std::shared_ptr<EvalCflFluxFunction> evalCflFluxFunction_;
    std::vector<LocalTimesteppingData> timeStepData_;
    bool localTimeStepping_;
    Scalar subCFLFactor_;
    Scalar accumulatedDt_;
    const Scalar dtThreshold_;
    int verbosity_;
};


/*! \brief Calculate the update vector.
 *  \param t         current time
 *  \param dt        time step size
 *  \param updateVec  vector containing the update values
 *  \param impet      variable should be true if an impet algorithm is used and false if the transport part is solved independently
 *
 *  Additionally to the \a update vector, the recommended time step size \a dt is calculated
 *  employing a CFL condition.
 */
template<class TypeTag>
void FVTransport<TypeTag>::update(const Scalar t, Scalar& dt, TransportSolutionType& updateVec, bool impet)
{
    if (!impet)
    {
        asImp_().updateMaterialLaws();
    }

    unsigned int size = problem_.gridView().size(0);
    if (localTimeStepping_)
    {
        if (timeStepData_.size() != size)
            timeStepData_.resize(size);
    }
    // initialize dt very large
    dt = std::numeric_limits<Scalar>::max();

    // resize update vector and set to zero
    updateVec.resize(size);
    updateVec = 0.0;

    // compute update vector
    for (const auto& element : elements(problem_.gridView()))
    {
#if HAVE_MPI
        if (element.partitionType() != Dune::InteriorEntity)
        {
            continue;
        }
#endif

        // cell index
        int globalIdxI = problem_.variables().index(element);

        CellData& cellDataI = problem_.variables().cellData(globalIdxI);

        Scalar update = 0;
        evalCflFluxFunction().reset();

        if (localTimeStepping_)
        {
            LocalTimesteppingData& localData = timeStepData_[globalIdxI];
            for (int i = 0; i < 2*dim; i++)
            {
                if (localData.faceTargetDt[i] < accumulatedDt_ + dtThreshold_)
                {
                    localData.faceFluxes[i] = 0.0;
                }
            }
        }

        // run through all intersections with neighbors and boundary
        for (const auto& intersection : intersections(problem_.gridView(), element))
        {
            GlobalPosition unitOuterNormal = intersection.centerUnitOuterNormal();
            if (switchNormals_)
                unitOuterNormal *= -1.0;

            int indexInInside = intersection.indexInInside();

            // handle interior face
            if (intersection.neighbor())
            {
                if (localTimeStepping_)
                {
                    LocalTimesteppingData& localData = timeStepData_[globalIdxI];

                    if (localData.faceTargetDt[indexInInside] < accumulatedDt_ + dtThreshold_)
                    {
                        asImp_().getFlux(localData.faceFluxes[indexInInside], intersection, cellDataI);
                    }
                    else
                    {
                        asImp_().getFlux(update, intersection, cellDataI);//only for time-stepping
                    }
                }
                else
                {
                    //add flux to update
                    asImp_().getFlux(update, intersection, cellDataI);
                }
            } //end intersection with neighbor element
            // handle boundary face
            else if (intersection.boundary())
            {
                if (localTimeStepping_)
                {
                    LocalTimesteppingData& localData = timeStepData_[globalIdxI];
                    if (localData.faceTargetDt[indexInInside] < accumulatedDt_ + dtThreshold_)
                    {
                        asImp_().getFluxOnBoundary(localData.faceFluxes[indexInInside], intersection, cellDataI);
                    }
                    else
                    {
                        asImp_().getFluxOnBoundary(update, intersection, cellDataI);//only for time-stepping
                    }
                }
                else
                {
                    //add boundary flux to update
                    asImp_().getFluxOnBoundary(update, intersection, cellDataI);
                }
            } //end boundary
        } // end all intersections

        if (localTimeStepping_)
        {
            LocalTimesteppingData& localData = timeStepData_[globalIdxI];
            for (int i=0; i < 2*dim; i++)
            {
                updateVec[globalIdxI] += localData.faceFluxes[i];
            }
        }
        else
        {
            //add flux update to global update vector
            updateVec[globalIdxI] += update;
        }

        //        std::cout<<"updateVec at "<<element.geometry().center()<<" : "<<updateVec[globalIdxI]<<"\n";

        //add source to global update vector
        Scalar source = 0.;
        asImp_().getSource(source,element, cellDataI);
        updateVec[globalIdxI] += source;

        //calculate time step
        using std::min;
        if (localTimeStepping_)
        {
            Scalar dtCfl = evalCflFluxFunction().getDt(element);

            timeStepData_[globalIdxI].dt = dtCfl;
            dt = min(dt, dtCfl);
        }
        else
        {
            //calculate time step
            dt = min(dt, evalCflFluxFunction().getDt(element));
        }

        //store update
        cellDataI.setUpdate(updateVec[globalIdxI]);
    } // end grid traversal


#if HAVE_MPI
    // communicate updated values
    using SolutionTypes = GetProp<TypeTag, Properties::SolutionTypes>;
    using ElementMapper = typename SolutionTypes::ElementMapper;
    using DataHandle = VectorCommDataHandleEqual<ElementMapper, Dune::BlockVector<Dune::FieldVector<Scalar, 1> >, 0/*elementCodim*/>;
    DataHandle dataHandle(problem_.elementMapper(), updateVec);
    problem_.gridView().template communicate<DataHandle>(dataHandle,
                                                         Dune::InteriorBorder_All_Interface,
                                                         Dune::ForwardCommunication);

    if (localTimeStepping_)
    {
    using TimeDataHandle = VectorCommDataHandleEqual<ElementMapper, std::vector<LocalTimesteppingData>, 0/*elementCodim*/>;

    TimeDataHandle timeDataHandle(problem_.elementMapper(), timeStepData_);
    problem_.gridView().template communicate<TimeDataHandle>(timeDataHandle,
                                                         Dune::InteriorBorder_All_Interface,
                                                         Dune::ForwardCommunication);
    }

    dt = problem_.gridView().comm().min(dt);
#endif
}

template<class TypeTag>
void FVTransport<TypeTag>::updatedTargetDt_(Scalar &dt)
{
    dt = std::numeric_limits<Scalar>::max();

    // update target time-step-sizes
    for (const auto& element : elements(problem_.gridView()))
    {
#if HAVE_MPI
        if (element.partitionType() != Dune::InteriorEntity)
        {
            continue;
        }
#endif

        // cell index
        int globalIdxI = problem_.variables().index(element);

        LocalTimesteppingData& localDataI = timeStepData_[globalIdxI];


        using FaceDt = std::unordered_map<int, Scalar>;
        FaceDt faceDt;

        // run through all intersections with neighbors and boundary
        for (const auto& intersection : intersections(problem_.gridView(), element))
        {
            int indexInInside = intersection.indexInInside();
            using std::min;
            if (intersection.neighbor())
            {
                auto neighbor = intersection.outside();
                int globalIdxJ = problem_.variables().index(neighbor);

                int levelI = element.level();
                int levelJ = neighbor.level();

                if (globalIdxI < globalIdxJ && levelI <= levelJ)
                {
                    LocalTimesteppingData& localDataJ = timeStepData_[globalIdxJ];

                    int indexInOutside = intersection.indexInOutside();

                    if (localDataI.faceTargetDt[indexInInside] < accumulatedDt_ + dtThreshold_
                        || localDataJ.faceTargetDt[indexInOutside] < accumulatedDt_ + dtThreshold_)
                    {
                        Scalar timeStep  = min(localDataI.dt, localDataJ.dt);

                        if (levelI < levelJ)
                        {
                            typename FaceDt::iterator it = faceDt.find(indexInInside);
                            if (it != faceDt.end())
                            {
                                it->second = min(it->second, timeStep);
                            }
                            else
                            {
                                faceDt.insert(std::make_pair(indexInInside, timeStep));
                            }
                        }
                        else
                        {
                            localDataI.faceTargetDt[indexInInside] += subCFLFactor_ * timeStep;
                            localDataJ.faceTargetDt[indexInOutside] += subCFLFactor_ * timeStep;
                        }

                        dt = min(dt, timeStep);
                    }
                }
            }
            else if (intersection.boundary())
            {
                if (localDataI.faceTargetDt[indexInInside] < accumulatedDt_ + dtThreshold_)
                {
                    localDataI.faceTargetDt[indexInInside] += subCFLFactor_ * localDataI.dt;
                    dt =min(dt, subCFLFactor_ * localDataI.dt);
                }
            }
        }
        if (faceDt.size() > 0)
        {
            typename FaceDt::iterator it = faceDt.begin();
            for (;it != faceDt.end();++it)
            {
                localDataI.faceTargetDt[it->first] += subCFLFactor_ * it->second;
            }

            for (const auto& intersection : intersections(problem_.gridView(), element))
            {
                if (intersection.neighbor())
                {
                    int indexInInside = intersection.indexInInside();

                    it = faceDt.find(indexInInside);
                    if (it != faceDt.end())
                    {
                        int globalIdxJ = problem_.variables().index(intersection.outside());

                        LocalTimesteppingData& localDataJ = timeStepData_[globalIdxJ];

                        int indexInOutside = intersection.indexInOutside();

                        localDataJ.faceTargetDt[indexInOutside] += subCFLFactor_ * it->second;
                    }
                }
            }
        }
    }

#if HAVE_MPI
    // communicate updated values
    using SolutionTypes = GetProp<TypeTag, Properties::SolutionTypes>;
    using ElementMapper = typename SolutionTypes::ElementMapper;
    using TimeDataHandle = VectorCommDataHandleEqual<ElementMapper, std::vector<LocalTimesteppingData>, 0/*elementCodim*/>;

    TimeDataHandle timeDataHandle(problem_.elementMapper(), timeStepData_);
    problem_.gridView().template communicate<TimeDataHandle>(timeDataHandle,
                                                         Dune::InteriorBorder_All_Interface,
                                                         Dune::ForwardCommunication);

    dt = problem_.gridView().comm().min(dt);
#endif
}

template<class TypeTag>
void FVTransport<TypeTag>::innerUpdate(TransportSolutionType& updateVec)
{
    if (localTimeStepping_)
    {
        Scalar realDt = problem_.timeManager().timeStepSize();

        Scalar subDt = realDt;

        updatedTargetDt_(subDt);

        Scalar accumulatedDtOld = accumulatedDt_;
        accumulatedDt_ += subDt;

        Scalar t = problem_.timeManager().time();

        if (accumulatedDt_ < realDt)
        {
            using std::min;
            while(true)
            {
                Scalar dtCorrection = min(0.0, realDt - accumulatedDt_);
                subDt += dtCorrection;

                if (verbosity_ > 0)
                    std::cout<<"    Sub-time-step size: "<<subDt<<"\n";

                TransportSolutionType transportedQuantity;
                asImp_().getTransportedQuantity(transportedQuantity);

                bool stopTimeStep = false;
                int size = transportedQuantity.size();
                for (int i = 0; i < size; i++)
                {
                    Scalar newVal = transportedQuantity[i] += updateVec[i] * subDt;
                    if (!asImp_().inPhysicalRange(newVal))
                    {
                        stopTimeStep = true;

                        break;
                    }
                }

#if HAVE_MPI
                int rank = 0;
                if (stopTimeStep)
                    rank = problem_.gridView().comm().rank();

                rank = problem_.gridView().comm().max(rank);
                problem_.gridView().comm().broadcast(&stopTimeStep,1,rank);
#endif


                if (stopTimeStep && accumulatedDtOld > dtThreshold_)
                {
                    problem_.timeManager().setTimeStepSize(accumulatedDtOld);
                    break;
                }
                else
                {
                    asImp_().setTransportedQuantity(transportedQuantity);
                }


                if (accumulatedDt_ >= realDt)
                {
                    break;
                }

                problem_.model().updateTransport(t, subDt, updateVec);

                updatedTargetDt_(subDt);

                accumulatedDtOld = accumulatedDt_;
                accumulatedDt_ += subDt;
            }
        }
        else
        {
            asImp_().updateTransportedQuantity(updateVec, realDt);
        }

        resetTimeStepData_();
    }
}
}
#endif
