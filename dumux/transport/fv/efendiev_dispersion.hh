// $Id$

#ifndef DUNE_EFENDIEV_DISPERSION_HH
#define DUNE_EFENDIEV_DISPERSION_HH

#include<dune/grid/utility/hierarchicsearch.hh>
#include<dumux/transport/fv/diffusivepart.hh>
#include<dumux/diffusion/diffusionproblem.hh>

namespace Dune
{

/**        \ingroup diffPart
 * @brief Class for dispersive 2 Phase flow based on subscale heterogeneities
 *
 * Class for dispersive 2 Phase flow based on papers "Modeling of subgrid effects in
 * coarse-scale simulations of transport in heterogeneous porous media" by Efendiev and Durlofsky
 * in Water Resources Research 2000 and "Numerical modeling of subgrid heterogeneity in two phase flow
 * simulations" by Efendiev and Durlovsky in Water Resources Research 2002. A brief description can also
 * be seen in "Multiscale Modeling of Multi-Phase-Multi-Component Processes in Heterogeneous Media" by
 * Niessner, PhD Thesis 2006 see:
 * http://elib.uni-stuttgart.de/opus/volltexte/2006/2769/pdf/opus_niessner_heft151.pdf
 */

template<class G, class RT, class TransportProblem, class DiffusionProblem>
class EfendievDurlovskiDispersion : public DiffusivePart<G,RT>
{
private:
    enum{dim = G::dimension};
    typedef typename G::Traits::template Codim<0>::Entity Entity;
    typedef typename IntersectionIteratorGetter<G,LevelTag>::IntersectionIterator IntersectionIterator;
    typedef typename G::template Codim<1>::EntityPointer EntityPointer;
    typedef BlockVector< BlockVector<FieldVector<RT, dim> > > VelType;
    typedef FieldVector<RT, dim> FieldVector;

    G& grid;
    TransportProblem & tranProb;
    DiffusionProblem & diffProb;
    VelType stddev;
    static const double sigma = 1.5;
    const typename G::Traits::LevelIndexSet& isetC;
    const typename G::Traits::LevelIndexSet& isetFine;
    bool linear;
    const int finelevel_, coarselevel_;

public:
    // constructor
    /**
     * @param transport object must provide a member 'transport.problem.vTotal(Entity, int)'
     *  with Entity being a codim 0 element and a member 'transport.problem.materialLaw.fractionalW(RT)'
     * @param diffprob if the transport object t does not contain a member of type Diffusion or derived, you should
     * provide another DiffusionProblem here.
     * @param lin linear switch: if transport of a component is considered
     */
    EfendievDurlovskiDispersion(TransportProblem& t, DiffusionProblem& d, G& g, int level_coarse, int level_fine, bool lin=false):
        grid(g), isetC(g.levelIndexSet(level_coarse)), isetFine(g.levelIndexSet(level_fine)),
        diffProb(d), tranProb(t), stddev(), finelevel_(level_fine), coarselevel_(level_coarse)
    {
        PermStddev(d);
        stddev.resize(isetC.size(0));
        int endi = stddev.size();
        for (int i = 0; i<endi; i++) stddev[i].resize(2*dim);
    }

    //! Traces back Streamlines from any point in domain
    Dune::FieldVector<RT,G::dimension>
    StreamlineBackTrack( FieldVector faceglobal, const RT& time) const;

    virtual Dune::FieldVector<RT,G::dimension> operator() (const Entity& entity, int numberInSelf,
                                                           const RT satIntersection, const FieldVector& satGradient, const RT time) const;

private:
    Dune::FieldVector<RT,G::dimension> PermStddev(DiffusionProblem& problem);
};

/*! @brief Provides a dispersive flux at the given intersection
 *
 * Based on the papers by Efendiev and Durlovsky, this function calculates a dispersive flux for a given Intersection.
 * Thereto not only the given parameters but furthermore the fractional flow function from the member transport.problem and
 * the permeability function is needed. The latter is called by the constructer of this class from the member transport.diffusion.diffusionproblem
 * @param entity a two dimensional codim-0 entity. Note: this function works only for rectangles!
 * @param is a LevelIntersectionIterator pointing to the face of interest
 * @param satIntersection the saturation at the intersection
 * @param satGradient the Gradient of saturation at the intersection
 * @param time the actual time. This is needed to determine the streamline length
 * @return dispersive flux
 */
template<class G, class RT, class TransportProblem, class DiffusionProblem>
Dune::FieldVector<RT, G::dimension > EfendievDurlovskiDispersion<G, RT, TransportProblem, DiffusionProblem> :: operator() (const Entity& entity, const int numberInSelf,
                                                                                                                           const RT satIntersection, const FieldVector& satGradient, const RT time) const
{
    int index = isetC.index(entity);
    // get value and derivatives of fractional flow function for satIntersection
    double Si = std::max(satIntersection,1e-5);
    Si = std::min(Si,1-1e-5);
    double fw_i = tranProb.materialLaw.fractionalW(Si - 1e-5);
    double fw_j = tranProb.materialLaw.fractionalW(Si);
    double fw_k = tranProb.materialLaw.fractionalW(Si + 1e-5);
    // central derivative
    double fw_s = (fw_k-fw_i)/2e-5;
    // second derivative
    double fw_ss = (fw_i+fw_k-2*fw_j) / (1e-5*1e-5);

    //get position of face center
    EntityPointer is = entity.template entity<1>(numberInSelf);
    Dune::GeometryType gtf = (*is).type();
    const Dune::FieldVector<RT,dim-1>&
        facelocal = Dune::ReferenceElements<RT,dim-1>::general(gtf).position(0,0);
    const Dune::FieldVector<RT,dim> faceglobal = (*is).geometry().global(facelocal);

    // streamline length
    Dune::FieldVector<RT,dim> length(StreamlineBackTrack(faceglobal, time));

    Dune::FieldVector<RT,dim> dispFlux;
    double satGrad = satGradient *length;
    Dune::FieldVector<RT,dim> fluct(stddev[index][numberInSelf]);
    for (int d = 0; d < dim; d++) fluct[d] = stddev[index][numberInSelf][d] * tranProb.vTotal(entity,numberInSelf)[d];
    if (fabs(fw_ss) < 0.0001 || linear)
    {
        if (linear) fw_j = satIntersection;
        double D = fabs(length*fluct);
        D *= pow(0.5*sigma,4);
        dispFlux = satGradient;
        dispFlux *= -D * fw_j;
    }
    else
    {
        double sigmafac = 0.2;
        dispFlux[0] = fw_j * fabs(pow(sigmafac * sigma, 4) * fluct[0]
                                  * fw_s * fw_s * (1 - exp(-fabs( (satGradient[0] * length[0]) * fw_ss) ) ) / fw_ss);
        dispFlux[1] = fw_j * fabs(pow(sigmafac * sigma, 4) * fluct[1]
                                  * fw_s * fw_s * (1 - exp(-fabs( (satGradient[1] * length[1]) * fw_ss) ) ) / fw_ss);
    }

    return dispFlux;
}

/*!
  StreamlineBackTrack evaluates the streamline starting from any point within the grid-domain.
  This function is written for rectangular, axis-parallel, conform grids in 2D!
  Since this is Backtracking, the input point global is taken to be the end of the Streamline at t=time (input).
  Then the Streamline is traced back to the point where it started at t=0.
  @param global  position to start backtracking from in global coordinates in a Dune::FieldVector
  @param time    total simulation time. Determines the maximum length of the streamlines
  @return vector increment between the point global (input) and the start point of the streamline
*/
template<class G, class RT, class TransportProblem, class DiffusionProblem>
FieldVector<RT, G::dimension > EfendievDurlovskiDispersion<G, RT, TransportProblem, DiffusionProblem> :: StreamlineBackTrack(FieldVector global, const RT& time) const
{
    // intersection iterator type
    typedef typename G::template Codim<0>::LevelIntersectionIterator IntersectionIterator;

    Dune::FieldVector<RT, G::dimension > length;
    Dune::GeometryType gtf;
    // exit point of the streamline
    Dune::FieldVector<RT,dim> pos_out_glob = global;
    Dune::FieldVector<RT,dim> pos_in_glob;
    // repetedely used coordinate-vectors
    Dune::FieldVector<RT,dim>  faceglobal, faceInElem;
    Dune::FieldVector<RT,dim-1> facelocal;
    // Hierarchic Search object for searching first grid cell
    Dune::HierarchicSearch<G,typename G::template Codim<0>::LevelIndexSet> child( grid, isetC);
    // a pointer to the first grid cell to start backtracking from
    typename G::template Codim<0>::EntityPointer nextElem = child.findEntity(global);


    bool stop = false;
    bool flag = false; // when flag == true, streamline info is printed to StreamLineTrace.dat this is for debugging purposes.
    double delta_t=0.;
    double Lx = 0.;
    double Ly = 0.;

    // stuff needed to write streamline info to a file...
    static bool first_time = true;
    // delete data in traceFile when functions is run the first time
    if (first_time)
    {
        std::fstream traceFile("StreamLineTrace.dat", std::ios::out|std::ios::trunc);
        traceFile.close();
    }
    first_time = false;
    std::fstream traceFile;
    traceFile.open("StreamLineTrace.dat", std::ios::out|std::ios::app);
    if (flag) traceFile << "start streamline backtracking from point [ " << global[0] << " , " << global[1] <<" ] " << " in Element " << isetC.index(*nextElem)<<"  at time "<< time<<std::endl;
    int countloops = 0;

    if (time == 0.)
    {
        length[0] = 0.;
        length[1] = 0.;
        if(flag) traceFile<< "time is zero, abort backtracking" <<std::endl;
        traceFile.close();
        return length;
    }

    // start streamline backtracking (this is gonna be nasty!)
    // since this is backtracking, we search the entry point starting from the exit point.
    while(!stop)
    {
        countloops++;
        if (flag)
        {
            if (countloops>1) traceFile << "proceeding through element " << isetC.index(*nextElem)<<std::endl;
        }
        Dune::FieldVector<RT,dim> pos_out = nextElem->geometry().local(pos_out_glob); // exit point
        Dune::FieldVector<RT,dim> pos_in; // entry point
        double x0, x1, y0, y1; // global positions of the faces of the cell
        double vx0, vx1, vy0, vy1; // velocities at the faces
        const int index = isetC.index(*nextElem);

        // first Intersection traversal for determining the current cell's velocity field:
        IntersectionIterator endis = nextElem->ilevelend();
        for (IntersectionIterator is = nextElem->ilevelbegin(); is != endis; ++is)
        {

            gtf = is.intersectionSelfLocal().type();
            facelocal = Dune::ReferenceElements<RT,dim-1>::general(gtf).position(0,0);
            faceglobal = is.intersectionGlobal().global(facelocal);
            faceInElem = nextElem->geometry().local(faceglobal);
            int numberInSelf = is.numberInSelf();

            // calculate total velocities at faces of current cell
            if (faceInElem[0] < 1e-5) // left face
            {
                x0 = is.intersectionGlobal().global(facelocal)[0];
                vx0 = tranProb.vTotal(*nextElem, numberInSelf)[0];
            }
            if (faceInElem[0] > 1. - 1e-5) // right face
            {
                x1 = is.intersectionGlobal().global(facelocal)[0];
                vx1 = tranProb.vTotal(*nextElem, numberInSelf)[0] ;
            }
            if (faceInElem[1] < 1e-5) // lower face
            {
                y0 = is.intersectionGlobal().global(facelocal)[1];
                vy0 = tranProb.vTotal(*nextElem, numberInSelf)[1] ;
            }
            if (faceInElem[1] > 1. - 1e-5) // upper face
            {
                y1 = is.intersectionGlobal().global(facelocal)[1];
                vy1 = tranProb.vTotal(*nextElem, numberInSelf)[1] ;
            }
        }
        // further computations are performed for a unit-cell --> normalise velocities!
        vx0 /= fabs(x0-x1);
        vx1 /= fabs(x0-x1);
        vy0 /= fabs(y0-y1);
        vy1 /= fabs(y0-y1);

        // variables characterising the cell's velocity field:
        const double A = vx1 - vx0 ;
        const double B = vy1 - vy0 ;

        // vector of the velocity at the exit point
        FieldVector v_out;
        v_out[0] = vx0   + pos_out[0] * A;
        v_out[1] = vy0   + pos_out[1] * B;

        IntersectionIterator inface = endis; // endis is only a dummy instantiation!
        double TOF = 1e99; // instatiate time of flight very high

        std::cout.precision(8);
        std::cout.setf (std::ios::scientific, std::ios::floatfield);

        // second Intersection traversal for determining the time of flight in the current cell
        for (IntersectionIterator is = nextElem->ilevelbegin(); is != endis; ++is)
        {
            //get time of flight (TOF)
            double tof;
            gtf = is.intersectionSelfLocal().type();
            facelocal = Dune::ReferenceElements<RT,dim-1>::general(gtf).position(0,0);
            faceglobal = is.intersectionGlobal().global(facelocal);
            faceInElem = nextElem->geometry().local(faceglobal);
            if (faceInElem[0] < 1e-5)
            {
                if ( fabs((vx0-vx1)/(vx0+vx1)) < 1e-5 || A==0. )
                {
                    tof = 2*(pos_out[0] - faceInElem[0]) / (vx0 + vx1);
                }
                else
                {
                    tof = log(v_out[0] / vx0) / A;
                }
                if (tof < TOF && vx0 > 0.)
                {
                    TOF = fabs(tof);
                    inface = is;
                }
            }
            if (faceInElem[0] > 1. - 1e-5)
            {
                if ( fabs((vx0-vx1)/(vx0+vx1)) < 1e-5 || A==0. )
                {
                    tof = 2*(pos_out[0] - faceInElem[0]) / (vx0 + vx1);
                }
                else
                {
                    tof = log(v_out[0] / vx1) / A;
                }
                if (tof < TOF && vx1 < 0.)
                {
                    TOF = fabs(tof);
                    inface = is;
                }
            }
            if (faceInElem[1] < 1e-5)
            {
                if ( fabs((vy0-vy1)/(vy0+vy1)) < 1e-5 || B==0. )
                {
                    tof = 2*(pos_out[1] - faceInElem[1]) / (vy0 + vy1);
                }
                else
                {
                    tof = log(v_out[1] / vy0) / B;
                }
                if (tof < TOF && vy0 > 0.)
                {
                    TOF = fabs(tof);
                    inface = is;
                }
            }
            if (faceInElem[1] > 1. - 1e-5)
            {
                if ( fabs((vy0-vy1)/(vy0+vy1)) < 1e-5 || B==0. )
                {
                    tof = 2*(pos_out[1] - faceInElem[1]) / (vy0 + vy1);
                }
                else
                {
                    tof = log(v_out[1] / vy1) / B;
                }
                if (tof < TOF && vy1 < 0.)
                {
                    TOF = fabs(tof);
                    inface = is;
                }
            }
        }

        // determining inflow position
        if ( fabs((vx0-vx1)/(vx0+vx1)) < 1e-5 || A==0. )
        {
            pos_in[0] = pos_out[0] - 0.5 * TOF *  (vx0 + vx1);
        }
        else
        {
            pos_in[0] = pos_out[0] - v_out[0] * (exp(A * TOF) - 1.) / A;
        }

        if ( fabs((vy0-vy1)/(vy0+vy1)) < 1e-5 || B==0. )
        {
            pos_in[1] = pos_out[1] - 0.5 * TOF *  (vy0 + vy1);
        }
        else
        {
            pos_in[1] = pos_out[1] - v_out[1] * (exp(B * TOF) - 1.) / B;
        }
        // regularisation to make sure that entry point lies on inflow face
        gtf = inface.intersectionSelfLocal().type();
        facelocal = Dune::ReferenceElements<RT,dim-1>::general(gtf).position(0,0);
        if (facelocal[0] == 0.) pos_in = 0.;
        if (facelocal[0] == 1.) pos_in = 1.;
        if (facelocal[1] == 0.) pos_in = 0.;
        if (facelocal[1] == 1.) pos_in = 1.;
        // regularisation to make sure that entry point lies in cell
        if (pos_in[0] < 0.) pos_in[0] = 0.;
        if (pos_in[0] > 1.) pos_in[0] = 1.;
        if (pos_in[1] < 0.) pos_in[1] = 0.;
        if (pos_in[1] > 1.) pos_in[1] = 1.;

        // Update of streamline length
        pos_in_glob = nextElem -> geometry().global(pos_in);
        pos_out_glob = nextElem -> geometry().global(pos_out);
        double factor = 1.;
        if (delta_t + TOF >= time)
        {
            factor = (time - delta_t) / TOF;
            stop = true; // stop backtracking if current simulation time is reached
        }
        Lx += factor*(pos_out_glob[0]-pos_in_glob[0]);
        Ly += factor*(pos_out_glob[1]-pos_in_glob[1]);
        delta_t += TOF;

        // set next element to continue backtracking in or stop backtracking if boundary is reached.
        if (inface.neighbor())
        {
            nextElem = inface.outside();
        }
        else stop = true;

        // set outflow position of next element to inflow position of current element.
        pos_out_glob = pos_in_glob;
        if (flag)
        {
            traceFile<<"entry point [ "<<pos_in_glob[0]<<" , "<<pos_in_glob[1]<<" ]"<<
                "  streamline length: [ "<< Lx<<" , "<<Ly <<" ] "<<" TOF: "<<delta_t<< std::endl;
        }
    }
    if (flag) traceFile<< std::endl;
    traceFile.close();
    length[0] = Lx;
    length[1] = Ly;
    return length;
}

// this class computes the standard deviations of permeability for every coarse scale element-face.
// these standard deviations are needed for the calculation of the dispersive flux.
template<class G, class RT, class TransportProblem, class DiffusionProblem>
Dune::FieldVector<RT,G::dimension> EfendievDurlovskiDispersion<G, RT, TransportProblem, DiffusionProblem> :: PermStddev(DiffusionProblem& problem)
{
    int maxlevel = grid.maxLevel();
    const int blocksize = 2*dim;
    const int nFine = (int)pow(2, ((finelevel_ - coarselevel_)*(dim-1)) ); //fine element edges per coarse element edge
    int hitcount[blocksize];
    Dune::BlockVector<BlockVector<FieldVector> > fineOnCoarse(blocksize);
    for (int i = 0; i < blocksize; i++) fineOnCoarse[i].resize(nFine);
    VelType deviation(isetC.size(0));
    int endi = deviation.size();
    for (int i = 0; i<endi; i++) deviation[i].resize(2*dim);

    typedef typename G::template Codim<0>::LevelIterator ElementLevelIterator;
    typedef typename G::template Codim<0>::LevelIntersectionIterator LevelIntersectionIterator;
    typedef typename G::template Codim<0>::LeafIntersectionIterator IntersectionIterator;
    typedef typename G::template Codim<0>::HierarchicIterator HierarchicIterator;
    typedef typename G::template Codim<0>::EntityPointer ElementEntityPointer;

    ElementLevelIterator endcit = grid.template lend<0>(coarselevel_);
    for (ElementLevelIterator cit = grid.template lbegin<0>(coarselevel_); cit!=endcit; ++cit)
    {
        for (int i=0; i<blocksize ; i++)  hitcount[i] = 0 ;

        for (HierarchicIterator it = cit->hbegin(maxlevel); it != cit-> hend(finelevel_); ++it)
        {
            if (isetFine.contains(*it))
            {

                // get some cell properties
                Dune::GeometryType gt = it->geometry().type();
                const Dune::FieldVector<RT,dim>&
                    local = Dune::ReferenceElements<RT,dim>::general(gt).position(0,0);
                Dune::FieldVector<RT,dim> global = it->geometry().global(local); // cell center in global coordinates
                double volume = it->geometry().integrationElement(local)
                    *Dune::ReferenceElements<RT,dim>::general(gt).volume();
                int indexi = isetFine.index(*it);

                // run through all intersections with neighbors and boundary
                IntersectionIterator isend = it->ilevelend();
                for (IntersectionIterator is = it->ilevelbegin(); is!=isend; ++is)
                {
                    // get some face properties
                    Dune::GeometryType gtf = is.intersectionSelfLocal().type();
                    int numberInSelf = is.numberInSelf();
                    const Dune::FieldVector<RT,dim-1>&
                        facelocal = Dune::ReferenceElements<RT,dim-1>::general(gtf).position(0,0);
                    FieldVector unitOuterNormal = is.unitOuterNormal(facelocal);
                    Dune::FieldVector<RT,dim>
                        faceglobal = is.intersectionGlobal().global(facelocal); // center of face in global coordinates

                    // check if face lies on one of father's faces
                    for (LevelIntersectionIterator cis = cit->ilevelbegin(); cis != cit->ilevelend();  ++cis)
                    {
                        // checking if leafelement face *is lies on level-zero element face *cis
                        Dune::GeometryType cgtf = cis.intersectionSelfLocal().type();
                        Dune::FieldVector<RT,dim-1> cfacelocal = Dune::ReferenceElements<RT,dim-1>::general(cgtf).position(0,0);
                        bool isInCis;
                        if ( cis.neighbor() )
                        {
                            /* if face *is lies in face *cis, *is must be defined both in the inside and outside
                               element of *cis! for the inside element this is true anyway, because of the hierarchic
                               Iteration. */
                            isInCis = cis.outside()->geometry().checkInside( cis.outside()->geometry().local(faceglobal));
                            if (isInCis)
                            {
                                // get permeability
                                FieldMatrix<RT,dim,dim> Ki = problem.K(global, *it, local);

                                // access neighbor
                                ElementEntityPointer outside = is.outside();
                                int indexj = isetFine.index(*outside);
                                GeometryType nbgt = outside->geometry().type();
                                const FieldVector& nblocal = ReferenceElements<RT,dim>::general(nbgt).position(0,0);
                                FieldVector nbglobal = outside->geometry().global(nblocal); // global coordinate of neighbor's cell center

                                // get neighbor permeability
                                FieldMatrix<RT,dim,dim> Kj = problem.K(global, *it, local);

                                // compute vectorized permeabilities
                                FieldVector Kni(0);
                                FieldVector Knj(0);
                                Ki.umv(unitOuterNormal, Kni);
                                Kj.umv(unitOuterNormal, Knj);
                                // compute permeability normal to intersection and take harmonic mean
                                double K_n_i = Kni * unitOuterNormal;
                                double K_n_j = Knj * unitOuterNormal;
                                double Kn    = 2 * K_n_i * K_n_j / (K_n_i + K_n_j);
                                // compute permeability tangential to intersection and take arithmetic mean
                                FieldVector uON = unitOuterNormal;
                                FieldVector K_t_i = Kni - (uON *= K_n_i);
                                uON = unitOuterNormal;
                                FieldVector K_t_j = Knj - (uON *= K_n_j);
                                FieldVector Kt = (K_t_i += K_t_j);
                                Kt *= 0.5;
                                // Build vectorized averaged permeability
                                uON = unitOuterNormal;
                                FieldVector K = (Kt += (uON *=Kn));

                                int coarseNumberInSelf = cis.numberInSelf();
                                int hits =  hitcount[coarseNumberInSelf];
                                fineOnCoarse[coarseNumberInSelf][hits] = K;
                                hitcount[coarseNumberInSelf] += 1;
                            }
                        }
                        else if (cis.boundary() && is.boundary())
                        {
                            /* if face *is lies in face *cis and
                             *cis lies on a boundary, *is must lie on a boundry too!
                             to be sure, that it is the same boundary, check if the normals are the same. */
                            isInCis = ( cis.unitOuterNormal(cfacelocal) * is.unitOuterNormal(facelocal) ) > 0.999;
                            if (isInCis)
                            {
                                // get permeability
                                FieldMatrix<RT,dim,dim> Ki = problem.K(global, *it, local);
                                FieldVector Kni(0);
                                Ki.umv(unitOuterNormal, Kni);
                                int coarseNumberInSelf = cis.numberInSelf();
                                int hits =  hitcount[coarseNumberInSelf];
                                fineOnCoarse[coarseNumberInSelf][hits] = Kni;
                                hitcount[coarseNumberInSelf] += 1;
                            }
                        }
                    }
                } // end intersection traversal
            } // end the if( it2->isLeaf) environment
        }// end hierarchic iteration

        for (int i=0;i<blocksize;i++)
        {
            int h = hitcount[i];
            h++;
        }
        // evaluate mean permeabilities and standard deviations at coarse element edges
        // and write them to the velocity struct.
        int coarseIndex = isetC.index(*cit);
        for (int i = 0; i<blocksize; i++)
        {
            FieldVector mean(0);
            FieldVector stddev_(0);
            // the mean
            for (int j = 0; j<nFine; j++)
            {
                mean += fineOnCoarse[i][j];
            }
            mean/= nFine;

            // standard deviation
            for (int j = 0; j<nFine; j++)
            {
                for (int d = 0; d<dim; d++) stddev_[d] += pow((mean[d]-fineOnCoarse[i][j][d]),2);
            }
            for (int d = 0; d<dim; d++) if (mean[d] != 0) stddev_[d] = sqrt(stddev_[d]/nFine)/mean[d];
            // write to velocity struct
            deviation[coarseIndex][i] = stddev_;
        }
    } // end grid traversal
    stddev = deviation;
}

}
#endif
