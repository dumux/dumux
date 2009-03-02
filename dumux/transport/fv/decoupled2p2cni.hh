//$Id$

#ifndef DECOUPLED2P2CNI_HH
#define DECOUPLED2P2CNI_HH

// commons:
#include <float.h>
#include <cmath>
#include <string>
#include <fstream>
#include <dune/common/helpertemplates.hh>
#include <dune/common/typetraits.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include <dune/grid/utility/intersectiongetter.hh>
#include <dune/istl/bvector.hh>
#include <dune/common/fvector.hh>
#include <dune/common/stdstreams.hh>

// transport:
#include "dumux/transport/fv/numericalflux.hh"
#include "dumux/transport/fv/diffusivepart.hh"
#include "dumux/transport/transportproblem2p2cni.hh"

// pressure:
#include <dune/disc/operators/boundaryconditions.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>
#include "dumux/pardiso/pardiso.hh"


//! author: Jochen Fritz
// last change: 15.10.2008

namespace Dune
{

/*####################################################*
 *                                                                                                         *
 *     CLASS DECLARATION                                                            *
 *                                                                                                         *
 *####################################################*/

//! Implementation of a decoupled formulation of a non-isothermal two phase two component process
/** \author Jochen Fritz
 *
 */
template<class G, class RT>
class Decoupled2p2cni
{
    template<int dim>
    struct ElementLayout
    {
        bool contains (GeometryType gt)
        {
            return gt.dim() == dim;
        }
    };

    enum{dim = G::dimension};
    enum{dimworld = G::dimensionworld};

    // typedef mania...
    typedef typename G::LevelGridView GV;
    typedef typename G::Traits::template Codim<0>::Entity Entity;
    typedef typename GV::IndexSet IS;
    typedef typename GV::template Codim<0>::Iterator Iterator;
    typedef typename G::template Codim<0>::HierarchicIterator HierarchicIterator;
    typedef MultipleCodimMultipleGeomTypeMapper<G,IS,ElementLayout> EM;
    typedef typename G::template Codim<0>::EntityPointer EntityPointer;
    typedef typename IntersectionIteratorGetter<G,LevelTag>::IntersectionIterator IntersectionIterator;
    typedef typename G::ctype ct;
    typedef FieldMatrix<double,1,1> MB;
    typedef BCRSMatrix<MB> MatrixType;
    typedef FieldVector<double, 1> VB;
    typedef BlockVector<VB> Vector;
    typedef FieldVector<double,dim> R2;
    typedef BlockVector<R2> R3;
    typedef BlockVector<R3> LocVelType;

public:
    typedef BlockVector< FieldVector<RT,1> > RepresentationType;
    typedef BlockVector< FieldVector<FieldVector<RT, G::dimension>, 2*G::dimension> > VelocityType;

    // initialization procedure
    void initial()
    {
        dinfo << "start initialization procedure..." << std::endl;
        upd = 0;
        timestep = 0;
        problem.variables.pressure = 1e5;
        initializeMatrix();            dinfo << "matrix initialization" << std::endl;;
        initialguess();                    dinfo << "first saturation guess" << std::endl;
        pressure(true, 0);            dinfo << "rough pressure calculation" << std::endl;
        transportInitial();            dinfo << "first transport initialization" << std::endl;
        pressure(false,0);            dinfo << "exact pressure calculation" << std::endl;
        transportInitial();            dinfo << "final transport initialization" << std::endl;
        totalVelocity(0);
        dinfo << "initialization procedure done" << std::endl;
    }

    // main part of the timestep
    void update(double t, double& dt, RepresentationType& updateVec)
    {
        dinfo << "start update..." << std::endl;
        upd = 0;
        dinfo << "transport prediction..."<< std::endl;
        concentrationUpdate(t, dt, upd);    dinfo << "transport prediction done"<<std::endl;
        upd *= dt;
        timestep = dt;
        dinfo << "pressure calculation..." << std::endl;
        pressure(false, t);                          dinfo << "pressure calculation done" << std::endl;    dinfo << "velocity calculation..." << std::endl;
        totalVelocity(t);                            dinfo << "velocity calculation done" << std::endl; dinfo << "transport calculation..."<<std::endl;
        concentrationUpdate(t, dt, updateVec);       dinfo << "transport calculation done" << std::endl; dinfo << "update done" << std::endl;
    }

    int level()
    {
        return level_;
    }

    // returns transport primary variables
    RepresentationType& operator*()
    {
        return problem.variables.totalConcentration;
    }

    // calculates the pressure field
    void pressure(bool first, const RT t=0)
    {
        assemble(first, t);
        solve();
        return;
    }

    void totalVelocity(const RT t);

    int concentrationUpdate(const RT t, RT& dt, RepresentationType& updateVec);

    // must be called at end of the timestep
    void postupdate(double t, double dt);

    // graphical output
    void vtkout(const char* name, int k)
    {
        problem.variables.vtkout(name, k);
    }

private:
    // pressure functions
    //~~~~~~~~~~~~~~~~~~~
    // must be called due to the functionality of the BCRSMAtrix class
    void initializeMatrix();

    // assembles the matrix and right hand side vector for the pressure calculation
    void assemble(bool first, const RT t);

    void solve();

    // transport functions
    //~~~~~~~~~~~~~~~~~~~~
    // makes a first guess on the saturation for a first pressure calculation
    void initialguess();

    // initializes all variables as soon as the pressure is known
    void transportInitial();

    // calculates phase saturations, temperature and phase composition
    void up_flash(double C1, double C2, double p, double h, double cp_matrix, double& temp, double& nuw, double& Xw1, double& Xn1, double& hgas, double& hliquid);

    // inner function of up_flash
    inline void flashCalculation(double Z1, double p, double temp, double& nuw, double& Xw1, double& Xn1);

    // calculates total concentration and phase composition from saturation
    void satFlash(double sat, double p, double temp, double poro, double& C1, double& C2, double& Xw1, double& Xn1);

    // variables
    // common variables
    G& grid;
    int level_;
    const IS& indexset;
    EM elementmapper;
    double timestep;

    // variables for transport equation:
    TransportProblem2p2cni<G, RT>& problem;
    RepresentationType upd;

    bool reconstruct;
    const NumericalFlux<RT>& numFlux;
    const DiffusivePart<G, RT>& diffusivePart;
    double alphamax;

    // variables for pressure equation:
    MatrixType A;
    RepresentationType f;
    std::string solverName_;
    std::string preconditionerName_;

public:
    //! constructor
    /**
     * \param g the grid
     * \param prob the problem to be solved
     * \param lev the grid level on which the calculation is to be run
     * \param diffPart a diffusion class
     * \param rec if set to true, a linear reconstruction shall be done (not implemented yet)
     * \param amax control parameter for linear reconstruction
     * \param solverName "BiCSTAB", "CG" or "Loop"
     * \param preconditionerName "SeqILU0" (combine only with CG or BiCGSTAB solvers!) or "SeqPardiso" (combine only with Loop solver!)
     */
    Decoupled2p2cni(
                    G& g,
                    TransportProblem2p2cni<G, RT>& prob,
                    int lev = 0,
                    DiffusivePart<G, RT>& diffPart = *(new DiffusivePart<G, RT>),
                    bool rec = false,
                    double amax = 0.8,
                    const NumericalFlux<RT>& numFl = *(new Upwind<RT>),
                    const std::string solverName = "BiCGSTAB",
                    const std::string preconditionerName = "SeqILU0" )
        :    grid(g), level_(lev), indexset(g.levelView(lev).indexSet()), elementmapper(g, g.levelView(lev).indexSet()),
             problem(prob), reconstruct(rec), numFlux(numFl), diffusivePart(diffPart), alphamax(amax),
             A(g.size(lev, 0),g.size(lev, 0), (2*dim+1)*g.size(lev, 0), BCRSMatrix<MB>::random), f(g.size(lev, 0)),
             solverName_(solverName), preconditionerName_(preconditionerName)
    {
        problem.variables.volErr = 0;
        upd.resize(3*elementmapper.size());
    };

}; //end class declaration



/*####################################################*
 *                                                                                                         *
 *     FUNCTION DEFINITIONS 1: PRESSURE EQUATION            *
 *                                                                                                         *
 *####################################################*/

template<class G, class RT>
void Decoupled2p2cni<G, RT>::initializeMatrix()
{
    // determine matrix row sizes
    Iterator eendit = grid.template lend<0>(level_);
    for (Iterator it = grid.template lbegin<0>(level_); it != eendit; ++it)
        {
            // cell index
            int indexi = elementmapper.map(*it);

            // initialize row size
            int rowSize = 1;

            // run through all intersections with neighbors
            IntersectionIterator endit = IntersectionIteratorGetter<G,LevelTag>::end(*it);
            for (IntersectionIterator is = IntersectionIteratorGetter<G,LevelTag>::begin(*it); is!=endit; ++is)
                if (is.neighbor()) rowSize++;
            A.setrowsize(indexi, rowSize);
        }
    A.endrowsizes();

    // determine position of matrix entries
    for (Iterator it = grid.template lbegin<0>(level_); it != eendit; ++it)
        {
            // cell index
            int indexi = elementmapper.map(*it);

            // add diagonal index
            A.addindex(indexi, indexi);

            // run through all intersections with neighbors
            IntersectionIterator endit = IntersectionIteratorGetter<G,LevelTag>::end(*it);
            for (IntersectionIterator is = IntersectionIteratorGetter<G,LevelTag>::begin(*it); is!=endit; ++is)
                if (is.neighbor())
                    {
                        // access neighbor
                        EntityPointer outside = is.outside();
                        int indexj = elementmapper.map(*outside);

                        // add off diagonal index
                        A.addindex(indexi, indexj);
                    }
        }
    A.endindices();

    return;
} // end function initialize matrix


/** for first == true, this function assembles the matrix and right hand side for
 * the solution of the pressure field in the same way as in the class FVDiffusion.
 * for first == false, the approach is changed to \f$ \[-\frac{\partial V}{\partial p}
 * \frac{\partial p}{\partial t}+\sum_{\kappa}\frac{\partial V}{\partial m^{\kappa}}\nabla\cdot
 * \left(\sum_{\alpha}C_{\alpha}^{\kappa}\mathbf{v}_{\alpha}\right)
 * =\sum_{\kappa}\frac{\partial V}{\partial m^{\kappa}}q^{\kappa}\] \f$. See Paper SPE 99619.
 * This is done to account for the volume effects which appear when gas and liquid are dissolved in each other.
 * \param first control switch for solution method
 * \param t simulation time
 */
template<class G, class RT>
void Decoupled2p2cni<G, RT>::assemble(bool first, const RT t=0)
{
    // initialization: set matrix A to zero
    A = 0;

    // iterate over all cells
    Iterator eendit = grid.template lend<0>(level_);
    for (Iterator it = grid.template lbegin<0>(level_); it != eendit; ++it)
        {
            // get geometry infos about the cell...
            GeometryType gt = it->geometry().type(); // cell geometry type
            const FieldVector<ct,dim>&
                local = ReferenceElements<ct,dim>::general(gt).position(0,0); // cell center in reference element
            FieldVector<ct,dim> global = it->geometry().global(local);             // global coordinate of cell center
            double volume = it->geometry().integrationElement(local)
                *ReferenceElements<ct,dim>::general(gt).volume(); // cell volume

            // cell index
            int indexi = elementmapper.map(*it);

            // get absolute permeability
            FieldMatrix<ct,dim,dim> Ki(this->problem.soil.K(global,*it,local));

            // get the cell's saturation
            double sati = problem.variables.saturation[indexi];

            // relative permeabilities
            std::vector<double> kr(problem.materialLaw.kr(sati, global, *it, local, problem.variables.temperature[indexi]));

            // phase viscosities
            double viscosityL, viscosityG;

            // total mobility and fractional flow factors
            double lambdaI, fw_I, fn_I;

            // derivatives of the fluid volume with respect to mass of compnents and pressure
            double dV_dm1 = 0;
            double dV_dm2 = 0;
            double dV_dh = 0;
            double dV_dp = 0;

            // specific volume of the phases
            double Vg, Vw;

            if (first)
                {
                    // total mobility and fractional flow factors
                    viscosityL = problem.liquidPhase.viscosity(problem.variables.temperature[indexi], problem.variables.pressure[indexi], 0.);
                    viscosityG = problem.gasPhase.viscosity(problem.variables.temperature[indexi], problem.variables.pressure[indexi], 0.);
                    lambdaI = kr[0] / viscosityL + kr[1] / viscosityG;
                    fw_I = kr[0] / viscosityL / lambdaI;
                    fn_I = kr[1] / viscosityG / lambdaI;

                    // specific volume of the phases
                    double Vg = 1. / problem.gasPhase.density(problem.variables.temperature[indexi][0], 1e5, 0.);
                    double Vw = 1. / problem.liquidPhase.density(problem.variables.temperature[indexi][0], 1e5, 0.);

                    FieldVector<RT,2> q = problem.q(global,*it,local);
                    f[indexi] = volume * (q[0] * Vw + q[1] * Vg);
                }
            else
                {
                    // get current temperature and pressure in cell
                    double T_I = problem.variables.temperature[indexi];
                    double p = problem.variables.pressure[indexi];

                    // total mobility and fractional flow factors
                    viscosityL = problem.liquidPhase.viscosity(T_I, p, 1. - problem.variables.wet_X1[indexi]);
                    viscosityG = problem.gasPhase.viscosity(T_I, p, problem.variables.nonwet_X1[indexi]);
                    lambdaI = kr[0] / viscosityL + kr[1] / viscosityG;
                    fw_I = kr[0] / viscosityL / lambdaI;
                    fn_I = kr[1] / viscosityG / lambdaI;

                    // specific volume of the phases
                    Vg = 1. / problem.variables.density_nonwet[indexi];
                    Vw = 1. / problem.variables.density_wet[indexi];

                    // matrix porosity
                    double poro = problem.soil.porosity(global, *it, local);

                    // mass of components inside the cell
                    double m1 = problem.variables.totalConcentration[indexi]*volume;
                    double m2 = problem.variables.totalConcentration[indexi+elementmapper.size()]*volume;
                    double h = problem.variables.totalConcentration[indexi+2 * elementmapper.size()] /*+ p * poro * (1 - problem.variables.saturation[indexi])*/;
                    // mass fraction of wetting phase
                    double nuw1 = sati / Vw / (sati/Vw + (1-sati)/Vg);
                    // actual fluid volume
                    double volalt = (m1+m2) * (nuw1 * Vw + (1-nuw1) * Vg);

                    // increments for numerical derivatives
                    double inc1 = (fabs(upd[indexi][0])*poro > 1e-8 /Vw) ?  upd[indexi][0]*poro : 1e-8/Vw;
                    double inc2 =(fabs(upd[indexi+elementmapper.size()][0])*poro > 1e-8 / Vg) ?  upd[indexi+elementmapper.size()][0]*poro : 1e-8 / Vg;
                    inc1 *= volume;
                    inc2 *= volume;
                    if (m1 + inc1 < 0 || m1 + inc1 > poro * volume / Vw) inc1 = 1e-8 / Vw  * volume;
                    if (m2 + inc2 < 0 || m2 + inc2 > poro * volume / Vg) inc2 = 1e-8 / Vg  * volume;

                    // numerical derivative of fluid volume with respect to mass of component 1
                    m1 +=  inc1;
                    double Z1 = m1 / (m1 + m2);
                    double dummy1, dummy2, nuw;
                    flashCalculation(Z1, p, T_I, nuw, dummy1, dummy2);
                    dV_dm1 = ((m1+m2) * (nuw * Vw + (1-nuw) * Vg) - volalt) /inc1;
                    m1 -= inc1;

                    // numerical derivative of fluid volume with respect to mass of component 2
                    m2 += inc2;
                    Z1 = m1 / (m1 + m2);
                    flashCalculation(Z1, p, T_I, nuw, dummy1, dummy2);
                    dV_dm2 = ((m1+m2) * (nuw * Vw + (1-nuw) * Vg) - volalt)/ inc2;
                    m2 -= inc2;

                    // numerical derivative of fluid volume with respect to pressure
                    double incp = 1e-5;
                    double p_ = p + incp;
                    double Vg_ = 1. / problem.gasPhase.density(T_I, p_, problem.variables.nonwet_X1[indexi]);
                    double Vw_ = 1. / problem.liquidPhase.density(T_I, p_, 1. - problem.variables.wet_X1[indexi]);
                    dV_dp = ((m1+m2) * (nuw1 * Vw_ + (1-nuw1) * Vg_) - volalt) /incp;

                    //numerical derivative of fluid volume with respect to enthalpy
                    Z1 = m1 / (m1 + m2);
                    inc1 /= volume;
                    inc2 /= volume;
                    double T_, Xw2, Xn1;
                    T_ = T_I +1;
                    double inch = problem.liquidPhase.intEnergy(T_I, p, 1. - problem.variables.wet_X1[indexi]) * poro * problem.variables.saturation[indexi] * problem.liquidPhase.density(T_I + 1, p, 1. - problem.variables.wet_X1[indexi])
                        + problem.gasPhase.intEnergy(T_I + 1, p, problem.variables.nonwet_X1[indexi]) * poro * (1 - problem.variables.saturation[indexi]) * problem.gasPhase.density(T_I + 1, p, problem.variables.nonwet_X1[indexi])
                        + problem.soil.heatCap(global, *it, local) * (T_I - 273.15) - problem.variables.totalConcentration[indexi+2 * elementmapper.size()] - h;
                    flashCalculation(Z1, p, T_, nuw, dummy1, dummy2);
                    dV_dh = ((m1+m2) * (nuw * Vw + (1-nuw) * Vg) - volalt)/ inch;

                    if (fabs(dV_dm1-0.001) > 0.001)
                        std::cout<< "cell " << indexi << " has dV_dm1 "<<dV_dm1 <<std::endl;
                    if (fabs(dV_dm2-1) > 1)
                        std::cout<< "cell " << indexi << " has dV_dm2 "<<dV_dm2 <<std::endl;

                    // right hand side entry: sources
                    FieldVector<RT,2> q = problem.q(global,*it,local);
                    double qh = problem.qh(global,*it,local);
                    f[indexi] = volume * (dV_dm1 * q[0] + dV_dm2 * q[1] + dV_dh * qh);
                }

            // iterate over all faces of the cell
            IntersectionIterator endit = IntersectionIteratorGetter<G,LevelTag>::end(*it);
            for (IntersectionIterator is = IntersectionIteratorGetter<G,LevelTag>::begin(*it);
                 is!=endit; ++is)
                {
                    // some geometry informations of the face
                    GeometryType gtf = is.intersectionSelfLocal().type(); // get geometry type of face
                    const FieldVector<ct,dim-1>&
                        facelocal = ReferenceElements<ct,dim-1>::general(gtf).position(0,0); // center in face's reference element
                    const FieldVector<ct,dim>&
                        facelocalDim = ReferenceElements<ct,dim>::general(gtf).position(is.numberInSelf(),1); // center of face inside volume reference element
                    FieldVector<ct,dimworld> unitOuterNormal = is.unitOuterNormal(facelocal);// get normal vector
                    FieldVector<ct,dimworld> integrationOuterNormal = is.integrationOuterNormal(facelocal);
                    integrationOuterNormal *= ReferenceElements<ct,dim-1>::general(gtf).volume(); //normal vector scaled with volume
                    double faceVol = is.intersectionGlobal().volume(); // get face volume

                    // compute directed permeability vector Ki.n
                    FieldVector<ct,dim> Kni(0);
                    Ki.umv(unitOuterNormal, Kni);

                    // handle interior face
                    if (is.neighbor())
                        {
                            // acces neighbor
                            EntityPointer outside = is.outside();
                            int indexj = elementmapper.map(*outside);

                            // some geometry infos of the neighbor
                            GeometryType nbgt = outside->geometry().type();
                            const FieldVector<ct,dim>& nblocal = ReferenceElements<ct,dim>::general(nbgt).position(0,0);
                            FieldVector<ct,dimworld>
                                nbglobal = outside->geometry().global(nblocal); // neighbor cell center in global coordinates

                            // distance vector between barycenters
                            FieldVector<ct,dimworld> distVec = global - nbglobal;

                            // compute distance between cell centers
                            double dist = distVec.two_norm();

                            // get absolute permeability
                            FieldMatrix<ct,dim,dim> Kj(problem.soil.K(nbglobal, *outside, nblocal));

                            // compute vectorized permeabilities
                            FieldVector<ct,dim> Knj(0);
                            Kj.umv(unitOuterNormal, Knj);
                            double K_n_i = Kni * unitOuterNormal;
                            double K_n_j = Knj * unitOuterNormal;
                            double Kn    = 2 * K_n_i * K_n_j / (K_n_i + K_n_j);
                            // compute permeability tangential to intersection and take arithmetic mean
                            FieldVector<ct,dim> uON = unitOuterNormal;
                            FieldVector<ct,dim> K_t_i = Kni - (uON *= K_n_i);
                            uON = unitOuterNormal;
                            FieldVector<ct,dim> K_t_j = Knj - (uON *= K_n_j);
                            FieldVector<ct,dim> Kt = (K_t_i += K_t_j);
                            Kt *= 0.5;
                            // Build vectorized averaged permeability
                            uON = unitOuterNormal;
                            FieldVector<ct,dim> K = (Kt += (uON *=Kn));

                            //compute mobilities
                            double fw_J, fn_J, lambdaJ;
                            double satj = problem.variables.saturation[indexj];
                            kr = problem.materialLaw.kr(satj, nbglobal, *outside, nblocal, problem.variables.temperature[indexj][0]);
                            if (first)
                                {
                                    viscosityL = problem.liquidPhase.viscosity(problem.variables.temperature[indexj][0], problem.variables.pressure[indexj], 0.);
                                    viscosityG = problem.gasPhase.viscosity(problem.variables.temperature[indexj][0], problem.variables.pressure[indexj], 0.);
                                    lambdaJ = kr[0] / viscosityL + kr[1] / viscosityG;
                                }
                            else
                                {
                                    viscosityL = problem.liquidPhase.viscosity(problem.variables.temperature[indexj][0], problem.variables.pressure[indexj], 1. - problem.variables.wet_X1[indexj]);
                                    viscosityG = problem.gasPhase.viscosity(problem.variables.temperature[indexj][0], problem.variables.pressure[indexj], problem.variables.nonwet_X1[indexj]);
                                    lambdaJ = kr[0] / viscosityL + kr[1] / viscosityG;
                                    fw_J = kr[0] / viscosityL / lambdaJ;
                                    fn_J = kr[1] / viscosityG / lambdaJ;
                                }

                            // compute averaged total mobility
                            // CAREFUL: Harmonic weightig can generate zero matrix entries,
                            // use arithmetic weighting instead:
                            double lambda;
                            lambda = 0.5*(lambdaI + lambdaJ);

                            // update diagonal entry
                            double entry;
                            if (first)
                                entry = fabs(lambda*faceVol*(K*distVec)/(dist*dist));
                            else
                                {
                                    // phase densities in cell in neighbor
                                    double rho_w_I = 1 / Vw;
                                    double rho_n_I = 1 / Vg;
                                    double rho_w_J = problem.variables.density_wet[indexj];
                                    double rho_n_J = problem.variables.density_nonwet[indexj];
                                    if (problem.variables.pressure[indexi] > problem.variables.pressure[indexj])
                                        {
                                            entry = fabs(
                                                         dV_dm1 * ( rho_w_I * problem.variables.wet_X1[indexi] * fw_I + rho_n_I * problem.variables.nonwet_X1[indexi] * fn_I)
                                                         + dV_dm2 * ( rho_w_I * (1. - problem.variables.wet_X1[indexi]) * fw_I + rho_n_I * (1.- problem.variables.nonwet_X1[indexi]) * fn_I)
                                                         + dV_dh  * ( rho_w_I * problem.variables.enthalpy_l[indexi] * fw_I + rho_n_I * problem.variables.enthalpy_g[indexi] * fn_I)
                                                         );
                                        }
                                    else
                                        {
                                            entry = fabs(
                                                         dV_dm1 * ( rho_w_J * problem.variables.wet_X1[indexj] * fw_J + rho_n_J * problem.variables.nonwet_X1[indexj] * fn_J)
                                                         + dV_dm2 * ( rho_w_J * (1. - problem.variables.wet_X1[indexj]) * fw_J + rho_n_J * (1. - problem.variables.nonwet_X1[indexj]) * fn_J)
                                                         + dV_dh  * ( rho_w_J * problem.variables.enthalpy_l[indexj] * fw_J + rho_n_J * problem.variables.enthalpy_g[indexj] * fn_J)
                                                         );
                                        }
                                    entry *= lambda * fabs(faceVol*(K*distVec)/(dist*dist));
                                }
                            // set diagonal entry
                            A[indexi][indexi] += entry;

                            // set off-diagonal entry
                            A[indexi][indexj] = -entry;
                        }

                    // boundary face
                    else
                        {
                            // center of face in global coordinates
                            FieldVector<ct,dimworld>
                                faceglobal = is.intersectionGlobal().global(facelocal);

                            // compute total mobility
                            double lambda = lambdaI;

                            //get boundary condition for boundary face center
                            BoundaryConditions::Flags bctype = problem.pbctype(faceglobal, *it, facelocalDim);

                            //dirichlet boundary
                            if (bctype == BoundaryConditions::dirichlet)
                                {
                                    // distance vector to boundary face
                                    FieldVector<ct,dimworld> distVec(global - faceglobal);
                                    double dist = distVec.two_norm();
                                    if (first)
                                        {
                                            A[indexi][indexi] -= lambda * faceVol * (Kni * distVec) / (dist * dist);
                                            double pressBC = problem.gPress(faceglobal, *it, facelocalDim);
                                            f[indexi] -= lambda * faceVol * pressBC * (Kni * distVec) / (dist * dist);
                                        }
                                    else
                                        {
                                            double pressBC = problem.gPress(faceglobal, *it, facelocalDim);

                                            double satBound, C1Bound, C2Bound, Xw1Bound, Xn1Bound;
                                            double TBound = problem.gT(faceglobal, *it, facelocalDim);

                                            //get boundary condition type for compositional transport
                                            BoundaryConditions2p2c::Flags bctype = problem.cbctype(faceglobal, *it, facelocalDim);
                                            if (bctype == BoundaryConditions2p2c::saturation) // saturation given
                                                {
                                                    satBound = problem.gS(faceglobal, *it, facelocalDim);
                                                    satFlash(satBound, pressBC, TBound, problem.soil.porosity(global, *it, local), C1Bound, C2Bound, Xw1Bound, Xn1Bound);
                                                }
                                            if (bctype == BoundaryConditions2p2c::concentration) // mass fraction given
                                                {
                                                    double Z1Bound = problem.gZ(faceglobal, *it, facelocalDim);
                                                    flashCalculation(Z1Bound, pressBC, TBound, satBound, Xw1Bound, Xn1Bound);
                                                }

                                            // phase densities in cell and on boundary
                                            double rho_w_I = 1 / Vw;
                                            double rho_n_I = 1 / Vg;
                                            double rho_w_J = problem.liquidPhase.density(TBound, pressBC, 1. - Xw1Bound);
                                            double rho_n_J = problem.gasPhase.density(TBound, pressBC, Xn1Bound);

                                            // phase enthalpies on boundary
                                            double enthalpy_w_Bound = problem.liquidPhase.enthalpy(TBound, pressBC, 1. - Xw1Bound);
                                            double enthalpy_n_Bound = problem.gasPhase.enthalpy(TBound, pressBC, Xn1Bound);

                                            double entry;
                                            if (problem.variables.pressure[indexi] > pressBC)
                                                entry = fabs(
                                                             dV_dm1 * ( rho_w_I * problem.variables.wet_X1[indexi] * fw_I + rho_n_I * problem.variables.nonwet_X1[indexi] * fn_I)
                                                             + dV_dm2 * ( rho_w_I * (1. - problem.variables.wet_X1[indexi]) * fw_I + rho_n_I * (1. - problem.variables.nonwet_X1[indexi]) * fn_I)
                                                             + dV_dh  * ( rho_w_I * problem.variables.enthalpy_l[indexi] * fw_I + rho_n_I * problem.variables.enthalpy_g[indexi] * fn_I)
                                                             );
                                            else
                                                entry = fabs(
                                                             dV_dm1 * ( rho_w_J * Xw1Bound * fw_I + rho_n_J * Xn1Bound * fn_I)
                                                             + dV_dm2 * ( rho_w_J * (1. - Xw1Bound) * fw_I + rho_n_J * (1. - Xn1Bound) * fn_I)
                                                             + dV_dh  * ( rho_w_J * enthalpy_w_Bound * fw_I + rho_n_J * enthalpy_n_Bound * fn_I)
                                                             );

                                            entry *= - lambda * faceVol*(Kni*distVec)/(dist*dist);

                                            // set diagonal entry and right hand side entry
                                            A[indexi][indexi] += entry;
                                            f[indexi] += entry * pressBC;
                                        }
                                }
                            else
                                {
                                    FieldVector<RT,3> J = problem.J(faceglobal, *it, facelocalDim);
                                    if (first)
                                        f[indexi] += faceVol * (J[0] * Vw + J[1] * Vg);
                                    else
                                        {
                                            f[indexi] += faceVol* (J[0] * dV_dm1 + J[1] * dV_dm2 + J[2] * dV_dh);
                                        }
                                }
                        }
                } // end all intersections

            // compressibility term
            if (!first && timestep != 0)
                {
                    A[indexi][indexi] -= dV_dp / timestep;
                    f[indexi] -= problem.variables.pressure[indexi] *dV_dp / timestep;
                }

            // error reduction routine: volumetric error is damped and inserted to right hand side
            // if damping is not done, the solution method gets unstable!
            double maxErr = fabs(problem.variables.volErr.infinity_norm());
            double erri = fabs(problem.variables.volErr[indexi]);
            double x_lo = 0.1;
            double x_mi = 0.99;
            double fac  = 0.9;
            double lofac = 0.;
            double hifac = 0.;
            hifac /= fac;

            //            if (erri > 5e-4)
            if (erri > x_lo * maxErr)
                {
                    if (erri <= x_mi * maxErr)
                        f[indexi] += fac* (1-x_mi*(lofac/fac-1)/(x_lo-x_mi) + (lofac/fac-1)/(x_lo-x_mi)*erri/maxErr)* problem.variables.volErr[indexi] / timestep * volume;
                    else
                        f[indexi] += fac * (1 + x_mi - hifac*x_mi/(1-x_mi) + (hifac/(1-x_mi)-1)*erri/maxErr) * problem.variables.volErr[indexi] / timestep * volume;
                }

        } // end grid traversal

    //        printmatrix(std::cout,A,"stiffnesmatrix","row");
    //        printvector(std::cout,f,"right hand side","row");
    return;
} // end function assemble


//! solves the system of equations that is made in the function assemble
/**
 *  which solver is used, is controlled with the constructor.
 */
template<class G, class RT>
void Decoupled2p2cni<G, RT>::solve()
{
    MatrixAdapter<MatrixType,Vector,Vector> op(A);
    InverseOperatorResult r;

    if (preconditionerName_ == "SeqILU0") {
        SeqILU0<MatrixType,Vector,Vector> preconditioner(A, 1.0);
        if (solverName_ == "CG") {
            CGSolver<Vector> solver(op, preconditioner, 1E-14, 10000, 0);
            solver.apply(problem.variables.pressure, f, r);
        }
        else if (solverName_ == "BiCGSTAB") {
            BiCGSTABSolver<Vector> solver(op, preconditioner, 1E-14, 10000, 0);
            solver.apply(problem.variables.pressure, f, r);
        }
        else
            DUNE_THROW(NotImplemented, "FVDiffusion :: solve : combination " << preconditionerName_
                       << " and " << solverName_ << ".");
    }
    else if (preconditionerName_ == "SeqPardiso") {
        SeqPardiso<MatrixType,Vector,Vector> preconditioner(A);
        if (solverName_ == "Loop") {
            LoopSolver<Vector> solver(op, preconditioner, 1E-14, 10000, 0);
            solver.apply(problem.variables.pressure, f, r);
        }
        else
            DUNE_THROW(NotImplemented, "FVDiffusion :: solve : combination " << preconditionerName_
                       << " and " << solverName_ << ".");
    }
    else
        DUNE_THROW(NotImplemented, "FVDiffusion :: solve : preconditioner " << preconditionerName_ << ".");

    return;
}

//! computes the total velocity from the pressure field.
template<class G, class RT>
void Decoupled2p2cni<G, RT>::totalVelocity(const RT t=0)
{
    // find out whether gravity effects are relevant
    bool hasGravity = false;
    const FieldVector<ct,dim>& gravity = problem.gravity();
    for (int k = 0; k < dim; k++)
        if (gravity[k] != 0)
            hasGravity = true;

    Iterator eendit = grid.template lend<0>(level_);
    for (Iterator it = grid.template lbegin<0>(level_); it != eendit; ++it)
        {
            // some geometry infos about the cell
            GeometryType gt = it->geometry().type();  // cell geometry type
            const FieldVector<ct,dim>& local = ReferenceElements<ct,dim>::general(gt).position(0,0); // cell center in reference element
            FieldVector<ct,dimworld> global = it->geometry().global(local);  // cell center in global coordinates

            // cell index
            int indexi = elementmapper.map(*it);

            // get pressure  in element
            double pressi = this->problem.variables.pressure[indexi];

            // get absolute permeability
            FieldMatrix<ct,dim,dim> Ki(problem.soil.K(global,*it,local));


            // total mobility and fractional flow factors
            double sati = problem.variables.saturation[indexi];
            std::vector<double> kr(problem.materialLaw.kr(sati, global, *it, local, problem.variables.temperature[indexi]));
            double viscosityL = problem.liquidPhase.viscosity(problem.variables.temperature[indexi], problem.variables.pressure[indexi], (1. - problem.variables.wet_X1[indexi]));
            double viscosityG = problem.gasPhase.viscosity(problem.variables.temperature[indexi], problem.variables.pressure[indexi], problem.variables.nonwet_X1[indexi]);
            double lambdaI = kr[0] / viscosityL + kr[1] / viscosityG;
            double fractionalWI = kr[0] / viscosityL / lambdaI;

            double faceVol[2*dim];

            // run through all intersections with neighbors and boundary
            IntersectionIterator endit = IntersectionIteratorGetter<G,LevelTag>::end(*it);
            for (IntersectionIterator is = IntersectionIteratorGetter<G,LevelTag>::begin(*it); is!=endit; ++is)
                {
                    // get geometry type of face
                    GeometryType gtf = is.intersectionSelfLocal().type();

                    //Geometry dg = is.intersectionSelfLocal();
                    // local number of facet
                    int numberInSelf = is.numberInSelf();

                    faceVol[numberInSelf] = is.intersectionGlobal().volume();

                    // center in face's reference element
                    const FieldVector<ct,dim-1>& facelocal = ReferenceElements<ct,dim-1>::general(gtf).position(0,0);

                    // center of face inside volume reference element
                    const FieldVector<ct,dim>& facelocalDim = ReferenceElements<ct,dim>::general(gtf).position(numberInSelf,1);

                    // get normal vector
                    FieldVector<ct,dimworld> unitOuterNormal = is.unitOuterNormal(facelocal);

                    // center of face in global coordinates
                    FieldVector<ct,dimworld> faceglobal = is.intersectionGlobal().global(facelocal);

                    // handle interior face
                    if (is.neighbor())
                        {
                            // access neighbor
                            EntityPointer outside = is.outside();
                            int indexj = elementmapper.map(*outside);

                            // get neighbor pressure and permeability
                            double pressj = this->problem.variables.pressure[indexj];

                            // geometry informations of neighbor
                            GeometryType nbgt = outside->geometry().type(); //geometry type
                            const FieldVector<ct,dim>& nblocal = ReferenceElements<ct,dim>::general(nbgt).position(0,0); // cell center in local coordinates
                            FieldVector<ct,dimworld> nbglobal = outside->geometry().global(nblocal); // neighbor cell center in global coordinates

                            // distance vector between barycenters
                            FieldVector<ct,dimworld> distVec = global - nbglobal;

                            // compute distance between cell centers
                            double dist = distVec.two_norm();

                            // get absolute permeability
                            FieldMatrix<ct,dim,dim> Kj(problem.soil.K(nbglobal, *outside, nblocal));

                            // compute vectorized permeabilities
                            FieldVector<ct,dim> Kni(0);
                            FieldVector<ct,dim> Knj(0);
                            Ki.umv(unitOuterNormal, Kni);
                            Kj.umv(unitOuterNormal, Knj);
                            // compute permeability normal to intersection and take harmonic mean
                            double K_n_i = Kni * unitOuterNormal;
                            double K_n_j = Knj * unitOuterNormal;
                            double Kn    = 2 * K_n_i * K_n_j / (K_n_i + K_n_j);
                            // compute permeability tangential to intersection and take arithmetic mean
                            FieldVector<ct,dim> uON = unitOuterNormal;
                            FieldVector<ct,dim> K_t_i = Kni - (uON *= K_n_i);
                            uON = unitOuterNormal;
                            FieldVector<ct,dim> K_t_j = Knj - (uON *= K_n_j);
                            FieldVector<ct,dim> Kt = (K_t_i += K_t_j);
                            Kt *= 0.5;
                            // Build vectorized averaged permeability
                            uON = unitOuterNormal;
                            FieldVector<ct,dim> K = (Kt += (uON *=Kn));


                            // total mobility and fractional flow factors
                            double satj = problem.variables.saturation[indexj];
                            kr = problem.materialLaw.kr(satj, nbglobal, *outside, nblocal, problem.variables.temperature[indexj]);
                            viscosityL = problem.liquidPhase.viscosity(problem.variables.temperature[indexj], problem.variables.pressure[indexj], 1. - problem.variables.wet_X1[indexj]);
                            viscosityG = problem.gasPhase.viscosity(problem.variables.temperature[indexj], problem.variables.pressure[indexj], problem.variables.nonwet_X1[indexj]);
                            double lambdaJ = kr[0] / viscosityL + kr[1] / viscosityG;
                            double fractionalWJ = kr[0] / viscosityL / lambdaJ;

                            // compute averaged total mobility
                            // CAREFUL: Harmonic weightig can generate zero matrix entries,
                            // use arithmetic weighting instead:
                            double lambda = 1;
                            double fractionalW;
                            lambda = 0.5 * (lambdaI + lambdaJ);
                            if (hasGravity)
                                fractionalW = 0.5 * (fractionalWI + fractionalWJ);

                            FieldVector<ct,dimworld> vTotal(K);
                            vTotal *= lambda * (pressi - pressj) / dist;
                            problem.variables.velocity[indexi][numberInSelf] = vTotal;
                        }
                    // boundary face
                    else
                        {
                            //get boundary condition for boundary face center
                            BoundaryConditions::Flags bctype = problem.pbctype(faceglobal, *it, facelocalDim);
                            if (bctype == BoundaryConditions::dirichlet)
                                {
                                    // uniform direction vector between barycenters
                                    FieldVector<ct,dimworld> distVec = global - faceglobal;
                                    double dist = distVec.two_norm();
                                    distVec /= dist;

                                    // compute directed permeability vector Ki.n
                                    FieldVector<ct,dim> Kni(0);
                                    Ki.umv(distVec, Kni);

                                    double pressBC = problem.gPress(faceglobal, *it, facelocalDim);

                                    FieldVector<ct,dim> vTotal(Kni);
                                    vTotal *= lambdaI * (pressBC - pressi) / dist;
                                    problem.variables.velocity[indexi][numberInSelf] = vTotal;
                                }
                            else
                                {
                                    // for Neumann boundary, only mass inflow but no volumetric flow is known.
                                    // To avoid the use of wrong values, set velocity to nonsens value.
                                    problem.variables.velocity[indexi][numberInSelf] = DBL_MAX; // highest number in double precission.
                                }
                        }
                } // end all intersections
        } // end grid traversal

    return;
} // end function totalVelocity


/*####################################################*
 *                                                                                                         *
 *     FUNCTION DEFINITIONS 2: TRANSPORT EQUATION            *
 *                                                                                                         *
 *####################################################*/

template<class G, class RT>
void Decoupled2p2cni<G,RT>::initialguess()
{
    // iterate through leaf grid an evaluate c0 at cell center
    Iterator eendit = grid.template lend<0>(level_);
    for (Iterator it = grid.template lbegin<0>(level_); it != eendit; ++it)
        {
            int indexi = elementmapper.map(*it);

            // get geometry information of cell
            GeometryType gt = it->geometry().type();
            const FieldVector<ct,dim>&
                local = ReferenceElements<ct,dim>::general(gt).position(0,0);
            FieldVector<ct,dimworld> global = it->geometry().global(local);

            // initial conditions
            double sat_0 = 0;
            double Z1_0 = 0;
            double T_0 = problem.T_0(global, *it, local); // get tempertaure initial condition
            BoundaryConditions2p2c::Flags ictype = problem.ictype(global, *it, local); // get type of initial condition

            if (ictype == BoundaryConditions2p2c::saturation)// saturation initial condition
                sat_0 = problem.S0(global, *it, local);
            else if(ictype == BoundaryConditions2p2c::concentration)            // saturation initial condition
                {
                    Z1_0 = problem.Z1_0(global, *it, local);
                    double rho_l = problem.liquidPhase.density(T_0, 1e5, 0.);
                    sat_0 = Z1_0 / rho_l;
                    sat_0 /= Z1_0 / rho_l + (1 - Z1_0) * problem.materialLaw.nonwettingPhase.density(T_0, 1e5, 0.);
                }
            else
                {
                    DUNE_THROW(Dune::NotImplemented, "Boundary condition " << ictype);
                }

            // initialize cell saturation
            this->problem.variables.saturation[indexi][0] = sat_0;
            //initialize temperature
            this->problem.variables.temperature[indexi][0] = T_0;
        }
    return;
}//end function initialguess



template<class G, class RT>
void Decoupled2p2cni<G,RT>::transportInitial()
{
    // iterate through leaf grid an evaluate c0 at cell center
    Iterator eendit = grid.template lend<0>(level_);
    for (Iterator it = grid.template lbegin<0>(level_); it != eendit; ++it)
        {
            int indexi = elementmapper.map(*it);

            // get geometry information of cell
            GeometryType gt = it->geometry().type();
            const FieldVector<ct,dim>&
                local = Dune::ReferenceElements<ct,dim>::general(gt).position(0,0);
            FieldVector<ct,dimworld> global = it->geometry().global(local);

            // initial conditions
            double sat_0, C1_0, C2_0;
            double p_0 = problem.variables.pressure[indexi];
            double T_0 = problem.variables.temperature[indexi][0];
            double poro = problem.soil.porosity(global, *it, local);
            double rho_w, rho_n;
            Dune::BoundaryConditions2p2c::Flags ictype = problem.ictype(global, *it, local);            // get type of initial condition

            if (ictype == Dune::BoundaryConditions2p2c::saturation)  // saturation initial condition
                {
                    sat_0 = problem.S0(global, *it, local);
                    satFlash(sat_0, p_0, T_0, poro, C1_0, C2_0, problem.variables.wet_X1[indexi][0], problem.variables.nonwet_X1[indexi][0]);
                    rho_w = problem.liquidPhase.density(T_0, p_0, 1. - problem.variables.wet_X1[indexi][0]);
                    rho_n = problem.gasPhase.density(T_0, p_0, problem.variables.nonwet_X1[indexi][0]);
                }
            else if (ictype == Dune::BoundaryConditions2p2c::concentration) // concentration initial condition
                {
                    double Z1_0 = problem.Z1_0(global, *it, local);
                    double nuw;
                    flashCalculation(Z1_0, p_0, T_0, nuw, problem.variables.wet_X1[indexi][0], problem.variables.nonwet_X1[indexi][0]);
                    rho_w = problem.liquidPhase.density(T_0, p_0, 1. - problem.variables.wet_X1[indexi][0]);
                    rho_n = problem.gasPhase.density(T_0, p_0, problem.variables.nonwet_X1[indexi][0]);
                    sat_0 = (nuw/rho_w) / (nuw/rho_w + (1-nuw)/rho_n);
                    C1_0 = poro* (sat_0 * problem.variables.wet_X1[indexi] * rho_w + (1-sat_0) * problem.variables.nonwet_X1[indexi] * rho_n);
                    C2_0 = poro* (sat_0 * (1. - problem.variables.wet_X1[indexi]) * rho_w + (1-sat_0) * (1. - problem.variables.nonwet_X1[indexi]) * rho_n);
                }

            // initialize cell concentration
            problem.variables.totalConcentration[indexi] = C1_0;
            problem.variables.totalConcentration[indexi + elementmapper.size()] = C2_0;
            problem.variables.saturation[indexi][0] = sat_0;
            problem.variables.density_wet[indexi] = rho_w;
            problem.variables.density_nonwet[indexi] = rho_n;

            // initialize phase enthalpies
            problem.variables.enthalpy_l[indexi][0] = problem.liquidPhase.enthalpy(T_0, p_0, 1. - problem.variables.wet_X1[indexi]);
            problem.variables.enthalpy_g[indexi][0] = problem.gasPhase.enthalpy(T_0, p_0, problem.variables.nonwet_X1[indexi]);
            problem.variables.totalConcentration[indexi + 2*elementmapper.size()] =
                problem.liquidPhase.intEnergy(T_0, p_0, 1. - problem.variables.wet_X1[indexi]) * poro * sat_0 * rho_w
                + problem.gasPhase.intEnergy(T_0, p_0, problem.variables.nonwet_X1[indexi]) * poro * (1 - sat_0) * rho_n
                + problem.soil.heatCap(global, *it, local) * (T_0 - 273.15);
        }
    return;
} //end function transportInitial


//! calculates the temporal derivatives of total concentrations and internal energy
template<class G, class RT>
int Decoupled2p2cni<G,RT>::concentrationUpdate(const RT t, RT& dt, RepresentationType& updateVec)
{
    // initialize timestep dt very large
    dt = 1E100;

    // set update vector to zero
    updateVec = 0;

    // stores number of cell which restricts timestep
    int which;

    // compute update vector
    Iterator eendit = grid.template lend<0>(level_);
    for (Iterator it = grid.template lbegin<0>(level_); it != eendit; ++it)
        {
            // get cell geometry informations
            GeometryType gt = it->geometry().type(); //geometry type
            const FieldVector<ct,dim>& local = ReferenceElements<ct,dim>::general(gt).position(0,0); // cell center in reference element
            FieldVector<ct,dimworld> global = it->geometry().global(local); // cell center in global coordinates
            double volume = it->geometry().integrationElement(local) * ReferenceElements<ct,dim>::general(gt).volume(); // cell volume, assume linear map here

            // cell index
            int indexi = elementmapper.map(*it);

            // get saturation and concentration value at cell center
            double satI = problem.variables.saturation[indexi];
            double Xw1_I = problem.variables.wet_X1[indexi];
            double Xn1_I = problem.variables.nonwet_X1[indexi];
            double Xw2_I = 1. - problem.variables.wet_X1[indexi];
            double Xn2_I = 1. - problem.variables.nonwet_X1[indexi];
            double T_I = problem.variables.temperature[indexi];

            // phase densities in cell
            double rho_w_I = problem.variables.density_wet[indexi];
            double rho_n_I = problem.variables.density_nonwet[indexi];

            // total mobility and fractional flow factors
            std::vector<double> kr(problem.materialLaw.kr(satI, global, *it, local, T_I));
            double viscosityL = problem.liquidPhase.viscosity(T_I, problem.variables.pressure[indexi], Xw2_I);
            double viscosityG = problem.gasPhase.viscosity(T_I, problem.variables.pressure[indexi], Xn1_I);
            double lambda = kr[0] / viscosityL + kr[1] / viscosityG;
            double fwI = kr[0] / viscosityL / lambda;
            double fnI = kr[1] / viscosityG / lambda;

            // some variables for time step calculation
            double sumfactor = 0;
            double sumfactor2 = 0;
            double sumDiff = 0;
            double sumDiff2 = 0;

            // run through all intersections with neighbors and boundary
            IntersectionIterator endit = IntersectionIteratorGetter<G,LevelTag>::end(*it);
            for (IntersectionIterator is = IntersectionIteratorGetter<G,LevelTag>::begin(*it);
                 is!=endit; ++is)
                {
                    // local number of facet
                    int numberInSelf = is.numberInSelf();

                    // get geometry informations of face
                    GeometryType gtf = is.intersectionSelfLocal().type(); //geometry type
                    const FieldVector<ct,dim-1>& facelocal = ReferenceElements<ct,dim-1>::general(gtf).position(0,0); // center in face's reference element
                    const FieldVector<ct,dim>& facelocalDim = ReferenceElements<ct,dim>::general(gtf).position(is.numberInSelf(),1); // center of face inside volume reference element

                    // get normal vector scaled with volume of face
                    FieldVector<ct,dimworld> integrationOuterNormal = is.integrationOuterNormal(facelocal);
                    integrationOuterNormal *= ReferenceElements<ct,dim-1>::general(gtf).volume();

                    // standardizes velocity
                    double velocityIJ = std::max(problem.variables.velocity[indexi][numberInSelf] * integrationOuterNormal / (volume), 0.0);

                    // variable for timestep calculation
                    double factor;

                    // variables which describe the change of the update vector per interface
                    double factorC1, factorC2, factorH;


                    if (is.neighbor()) // handle interior face
                        {
                            // access neighbor
                            EntityPointer outside = is.outside();
                            int indexj = elementmapper.map(*outside);

                            // neighbor geometry informations
                            GeometryType nbgt = outside->geometry().type();
                            const FieldVector<ct,dim>& nblocal = ReferenceElements<ct,dim>::general(nbgt).position(0,0);
                            FieldVector<ct,dimworld> nbglobal = outside->geometry().global(nblocal); // neighbor cell center in global coordinates

                            // standardized velocity
                            double velocityJI = std::max(-(problem.variables.velocity[indexi][numberInSelf] * integrationOuterNormal / volume), 0.0);

                            // distance vector between barycenters
                            FieldVector<ct,dimworld> distVec = global - nbglobal;

                            // distance between barycenters
                            double dist = distVec.two_norm();

                            // get saturation and concentration value at neighbor cell center
                            double satJ = this->problem.variables.saturation[indexj];
                            double Xw1_J = problem.variables.wet_X1[indexj];
                            double Xn1_J = problem.variables.nonwet_X1[indexj];
                            double Xw2_J = 1. - problem.variables.wet_X1[indexj];
                            double Xn2_J = 1. - problem.variables.nonwet_X1[indexj];
                            double T_J = problem.variables.temperature[indexj];

                            // phase densities in neighbor
                            double rho_w_J = problem.variables.density_wet[indexj];
                            double rho_n_J = problem.variables.density_nonwet[indexj];

                            // neighbor cell
                            viscosityL = problem.liquidPhase.viscosity(T_J, problem.variables.pressure[indexj], Xw2_J);
                            viscosityG = problem.gasPhase.viscosity(T_J, problem.variables.pressure[indexj], Xn1_J);
                            kr = problem.materialLaw.kr(satJ, nbglobal, *outside, nblocal, T_J);
                            lambda = kr[0] / viscosityL + kr[1] / viscosityG;
                            double fwJ = kr[0] / viscosityL / lambda;
                            double fnJ = kr[1] / viscosityG / lambda;

                            // for timestep control
                            {
                                double coeffW = fwI / satI;
                                if(isinf(coeffW) || isnan(coeffW)) coeffW = 0;
                                double coeffN = fnI / (1-satI);
                                if(isinf(coeffN) || isnan(coeffN)) coeffN = 0;
                                factor = (velocityJI - velocityIJ) * std::max(std::max(coeffW, coeffN),1.);
                            }

                            //  diffFactor =diffPart / volume; TODO include diffusion into timestep control
                            factorC1 =
                                velocityJI * Xw1_J * rho_w_J * numFlux(satJ, satI, fwJ, fwI)
                                - velocityIJ * Xw1_I * rho_w_I * numFlux(satI, satJ, fwI, fwJ)
                                + velocityJI * Xn1_J * rho_n_J * numFlux(1.0-satJ, 1.0-satI, fnJ, fnI)
                                - velocityIJ * Xn1_I * rho_n_I * numFlux(1.0-satI, 1.0-satJ, fnI, fnJ);
                            factorC2 =
                                velocityJI * Xw2_J * rho_w_J * numFlux(satJ, satI, fwJ, fwI)
                                - velocityIJ * Xw2_I * rho_w_I * numFlux(satI, satJ, fwI, fwJ)
                                + velocityJI * Xn2_J * rho_n_J * numFlux(1.0-satJ, 1.0-satI, fnJ, fnI)
                                - velocityIJ * Xn2_I * rho_n_I * numFlux(1.0-satI, 1.0-satJ, fnI, fnJ);
                            factorH =
                                velocityJI * problem.variables.enthalpy_l[indexj] * rho_w_J * numFlux(satJ, satI, fwJ, fwI)
                                - velocityIJ * problem.variables.enthalpy_l[indexi] * rho_w_I * numFlux(satI, satJ, fwI, fwJ)
                                + velocityJI * problem.variables.enthalpy_g[indexj] * rho_n_J * numFlux(1.0-satJ, 1.0-satI, fnJ, fnI)
                                - velocityIJ * problem.variables.enthalpy_g[indexi] * rho_n_I * numFlux(1.0-satI, 1.0-satJ, fnI, fnJ);

                            FieldVector<ct,dimworld> faceglobal = is.intersectionGlobal().global(facelocal);

                        }

                    else // handle boundary face
                        {
                            // cell center in global coordinates
                            FieldVector<ct,dimworld> global = it->geometry().global(local);

                            // face center in globel coordinates
                            FieldVector<ct,dim> faceglobal = is.intersectionGlobal().global(facelocal);

                            // distance vector between cell and face center
                            FieldVector<ct,dimworld> distVec = global - faceglobal;
                            double dist = distVec.two_norm();

                            // get saturation value at cell center
                            double satI = this->problem.variables.saturation[indexi];

                            // standardized velocity
                            double velocityJI = std::max(-(problem.variables.velocity[indexi][numberInSelf] * integrationOuterNormal / volume), 0.0);

                            //get boundary conditions
                            double TBound = problem.gT(faceglobal, *it, facelocalDim); // temperature BC

                            BoundaryConditions::Flags pressBCtype = problem.pbctype(faceglobal, *it, facelocalDim); // BC type
                            if (pressBCtype == BoundaryConditions::dirichlet)
                                {
                                    double pBound = problem.gPress(faceglobal, *it, facelocalDim); //pressure on Boundary

                                    double satBound, C1Bound, C2Bound, Xw1Bound, Xn1Bound;
                                    double rho_w_Bound, rho_n_Bound;
                                    BoundaryConditions2p2c::Flags bctype = problem.cbctype(faceglobal, *it, facelocalDim);
                                    if (bctype == BoundaryConditions2p2c::saturation)
                                        {
                                            satBound = problem.gS(faceglobal, *it, facelocalDim);
                                            satFlash(satBound, pBound, TBound, problem.soil.porosity(global, *it, local), C1Bound, C2Bound, Xw1Bound, Xn1Bound);
                                            rho_w_Bound = problem.liquidPhase.density(TBound, pBound, 1. - Xw1Bound);
                                            rho_n_Bound = problem.gasPhase.density(TBound, pBound, Xn1Bound);
                                        }
                                    if (bctype == BoundaryConditions2p2c::concentration)
                                        {
                                            double Z1Bound = problem.gZ(faceglobal, *it, facelocalDim);
                                            double nuw;
                                            flashCalculation(Z1Bound, pBound, TBound, nuw, Xw1Bound, Xn1Bound);
                                            rho_w_Bound = problem.liquidPhase.density(TBound, pBound, (1. - Xw1Bound));
                                            rho_n_Bound = problem.gasPhase.density(TBound, pBound, Xn1Bound);
                                            satBound = nuw / rho_w_Bound / (nuw / rho_w_Bound + (1-nuw) / rho_n_Bound);
                                        }

                                    // phase enthalpies on boundary
                                    double HwBound = problem.liquidPhase.enthalpy(TBound, pBound, (1. - Xw1Bound));
                                    double HnBound = problem.gasPhase.enthalpy(TBound, pBound, Xn1Bound);

                                    // on boundary
                                    kr = problem.materialLaw.kr(satBound, global, *it, local, TBound);
                                    lambda = kr[0] / viscosityL + kr[1] / viscosityG;
                                    double fwBound = kr[0] / viscosityL / lambda;
                                    double fnBound = kr[1] / viscosityG / lambda;

                                    // for timestep control
                                    {
                                        double coeffW = fwI / satI;
                                        if(isinf(coeffW) || isnan(coeffW)) coeffW = 0;
                                        double coeffN = fnI / (1-satI);
                                        if(isinf(coeffN) || isnan(coeffN)) coeffN = 0;
                                        factor = (velocityJI - velocityIJ) * std::max(std::max(coeffW, coeffN),1.);
                                    }

                                    factorC1 =
                                        + velocityJI * Xw1Bound * rho_w_Bound * numFlux(satBound, satI, fwBound, fwI)
                                        - velocityIJ * Xw1_I * rho_w_I * numFlux(satI, satBound, fwI, fwBound)
                                        + velocityJI * Xn1Bound * rho_n_Bound * numFlux(1.0-satBound, 1.0-satI, fnBound, fnI)
                                        - velocityIJ * Xn1_I * rho_n_I * numFlux(1.0-satI, 1.0-satBound, fnI, fnBound);
                                    factorC2 =
                                        velocityJI * (1. - Xw1Bound) * rho_w_Bound * numFlux(satBound, satI, fwBound, fwI)
                                        - velocityIJ * Xw2_I * rho_w_I * numFlux(satI, satBound, fwI, fwBound)
                                        + velocityJI * (1.- Xn1Bound) * rho_n_Bound * numFlux(1.0-satBound, 1.0-satI, fnBound, fnI)
                                        - velocityIJ * Xn2_I * rho_n_I * numFlux(1.0-satI, 1.0-satBound, fnI, fnBound);
                                    factorH =
                                        velocityJI * HwBound * rho_w_Bound * numFlux(satBound, satI, fwBound, fwI)
                                        - velocityIJ * problem.variables.enthalpy_l[indexi] * rho_w_I * numFlux(satI, satBound, fwI, fwBound)
                                        + velocityJI * HnBound * rho_n_Bound * numFlux(1.0-satBound, 1.0-satI, fnBound, fnI)
                                        - velocityIJ * problem.variables.enthalpy_g[indexi] * rho_n_I * numFlux(1.0-satI, 1.0-satBound, fnI, fnBound);
                                }
                            else if (pressBCtype == BoundaryConditions::neumann)
                                {
                                    FieldVector<RT,3> J = problem.J(faceglobal, *it, facelocalDim);
                                    double faceVol = integrationOuterNormal.two_norm();
                                    factorC1 = J[0] * faceVol / volume;
                                    factorC2 = J[1] * faceVol / volume;
                                    factorH = J[2] * faceVol / volume;

                                    // for timestep control
                                    double coeffW = fwI / satI;
                                    if(isinf(coeffW) || isnan(coeffW)) coeffW = 0;
                                    double coeffN = fnI / (1-satI);
                                    if(isinf(coeffN) || isnan(coeffN)) coeffN = 0;
                                    factor = fabs(J[0] * problem.liquidPhase.density(T_I, problem.variables.pressure[indexi][0], 0.) + J[1] * problem.gasPhase.density(T_I, problem.variables.pressure[indexi][0], 0.));
                                    factor *= std::max(std::max(coeffW, coeffN),1.);
                                }

                            else DUNE_THROW(NotImplemented, "there is no process boundary condition implemented");
                        }

                    // add to update vector
                    updateVec[indexi] += factorC1;
                    updateVec[elementmapper.size() + indexi] += factorC2;
                    updateVec[2*elementmapper.size() + indexi] += factorH;

                    // for time step calculation
                    if (factor>=0)
                        sumfactor += factor;
                    else
                        sumfactor2 += (-factor);
                } // end all intersections

            // account for porosity
            double poro = problem.soil.porosity(global, *it, local);

            // get source term
            updateVec[indexi] += problem.q(global, *it, local)[0];
            updateVec[indexi + elementmapper.size()] += problem.q(global, *it, local)[1];
            updateVec[indexi + 2*elementmapper.size()] += problem.qh(global, *it, local);

            sumfactor = std::max(sumfactor,sumfactor2) / poro;
            sumDiff = std::max(sumDiff,sumDiff2);
            sumfactor = std::max(sumfactor,100*sumDiff);
            if (1./sumfactor < dt)
                {
                    dt = 1./sumfactor;
                    which = indexi;
                }

        } // end grid traversal

    dverb << "timestep restricting cell: " << which << std::endl;
    //        printvector(std::cout,updateVec,"update","row");
    return 0;
} // end function "update"


//! isoenergetic flash
/**
 * \return C1 total concentration of component 1
 * \return C2 total concentration of component 2
 * \param p pressure
 * \param h enthalpy in control volume
 * \param cp_matrix matrix heat capacity
 * \param[out] temp temperature
 * \param[out] nuw mass fraction of water phase
 * \param[out] Xw1 mass fraction of component 1 in wetting phase
 * \param[out] Xw2 mass fraction of component 2 in wetting phase
 * \param[out] Xn1 mass fraction of component 1 in nonwetting phase
 * \param[out] Xn2 mass fraction of component 2 in nonwetting phase
 * \param[out] hgas enthalpy of the gas phase
 * \param[out] hliquid enthalpy of the liquid phase
 */
template<class G, class RT>
void Decoupled2p2cni<G,RT>::up_flash(double C1, double C2, double p, double h, double cp_matrix, double& temp, double& nuw, double& Xw1, double& Xn1, double& hliquid, double& hgas)
{
    dverb << "up_flash run for C1 = "<< C1 << ", C2 = "<< C2 <<", p = "<< p <<", h = "<< h << std::endl;
    // how this function works is explained in Agarwal, Li, Nghiem, Coombe: "Multiphase multicomponent isenthalpic flash calculation", Reservoir Engineering 1991, vol. 30, No. 3
    double Z1 = C1 / (C1 + C2);
    double mass = C1 + C2;
    double Th = problem.liquidPhase.T_vap(p);
    double sol = problem.liquidPhase.Xa_Max(Th, p);
    double rho_w = problem.liquidPhase.density(Th, p, sol);
    double rho_n = problem.gasPhase.density(Th, p, Z1);
    double gh = mass * problem.gasPhase.intEnergy(Th, p, Z1) + cp_matrix * (Th - 273.15) - h;
    double gi = mass * problem.liquidPhase.intEnergy(Th, p,sol) + cp_matrix * (Th - 273.15) - h;
    int count = 0;
    if (Z1 > 1-sol && gh > 0 && gi < 0) //for high concentrations of the liquid phase main component, the normal iteration would take too long.
        {
            double ggas = gh;
            Xw1 = 1;
            while (fabs(gh) > 1e-2)
                {
                    count ++;
                    nuw = ggas / (ggas - gi);
                    Xn1 = 1 - (1-Z1) / (1-nuw);
                    ggas = mass * problem.gasPhase.intEnergy(Th, p, Xn1) + cp_matrix * (Th - 273.15) - h;
                    gh = nuw * gi + (1-nuw) * ggas;
                }
            temp = Th;
            dverb << "took "<< count <<" iteration steps with pure liquid iteration sceme"<< std::endl;
        }
    else
        {
            double Ti, Tnew;
            double ugas, uliquid;
            Ti = gh>0 ? Th-1 : Th+1;

            double Tb[2] = {-1e20,1e20}; // boundarys, T has to lie between. Initialize far away
            double gb[2] = {-1e20,1e20}; // boundarys of target function at temperatures Tb
            bool b_[2] = {false,false};  // flags determining if boundarys already exist

            if (gh > 0) // Th can be used as first boundary
                {
                    Tb[1] = Th;
                    gb[1] = gh;
                    b_[1] = true;
                }
            else
                {
                    Tb[0] = Th;
                    gb[0] = gh;
                    b_[0] = true;
                }

            for (;;)
                {
                    flashCalculation(Z1, p, Ti, nuw, Xw1, Xn1);
                    uliquid = problem.liquidPhase.intEnergy(Ti, p, 1. - Xw1);
                    ugas = problem.gasPhase.intEnergy(Ti, p, Xn1);
                    gi = mass * nuw *  uliquid + mass * (1-nuw) * ugas + cp_matrix * (Ti - 273.15) - h;
                    count++;

                    if (fabs(gi) < 1e-2) break;

                    if (gi < 0 && Ti > Tb[0]) // can Ti be used as new lower boundary?
                        {
                            gb[0] = gi;
                            Tb[0] = Ti;
                            b_[0] = true;
                        }
                    else if (gi > 0 && Ti < Tb[1]) // can Ti be used as new upper boundary?
                        {
                            gb[1] = gi;
                            Tb[1] = Ti;
                            b_[1] = true;
                        }

                    Tnew = std::max(Ti - gi * (Ti-Th) / (gi-gh),1.);

                    if (b_[0] && b_[1] && (Tnew > Tb[1] || Tnew < Tb[0])) // if both boundaries exist and Ti lies outside, do alternative iteration step
                        Tnew = Tb[0] - gb[0]* (Tb[1]-Tb[0])/(gb[1]-gb[0]);

                    if (isnan(Tnew))
                        DUNE_THROW(RangeError, "Decoupled2p2cni::up_flash NANs appear!");
                    if (isinf(Tnew))
                        DUNE_THROW(RangeError, "Decoupled2p2cni::up_flash INFs appear!");

                    Th = Ti;
                    gh = gi;
                    Ti = Tnew;
                }// end loop
            temp = Ti;
            dverb << "took "<< count <<" iterations with standard iteration sceme"<<std::endl;
        }

    hliquid = problem.liquidPhase.enthalpy(temp, p, 1. - Xw1);
    hgas = problem.gasPhase.enthalpy(temp, p, Xn1);

    if (isnan(nuw + temp + Xw1 + Xn1))
        dgrave << "NAN in up_flash" << std::endl;
}

template<class G, class RT>
inline void Decoupled2p2cni<G,RT>::flashCalculation(double Z1, double p, double temp, double& nuw, double& Xw1, double& Xn1)
{
    // mole fraction equilibrium ratios
    double K1 = problem.liquidPhase.p_vap(temp) / p;
    double K2 = 1. / (p * problem.liquidPhase.henry(temp));
    // equilibrium mole fractions
    double xw1 = (1. - K2) / (K1 -K2);
    double xn1 = xw1 * K1;
    // regularization
    xw1 = std::min(1., std::max(xw1, 0.));
    xn1 = std::min(1., std::max(xn1, 0.));
    // Calculate mass fractions from mole fractions
    Xw1 = xw1 * problem.liquidPhase.molarMass_w()
        / ( xw1 * problem.liquidPhase.molarMass_w() + (1.-xw1) * problem.liquidPhase.molarMass_a() );
    Xn1 = xn1 * problem.liquidPhase.molarMass_w()
        / ( xn1 * problem.liquidPhase.molarMass_w() + (1.-xn1) * problem.liquidPhase.molarMass_a() );
    // mass fraction equilibrium ratios
    K1 = Xn1 / Xw1;
    K2 = (1.-Xn1) / (1.-Xw1);

    double nu = 0; // mass fraction of the gas phase

    if (Z1 > Xn1 && Z1 < Xw1)  // two phases appear
        nu = -((K1-1)*Z1 + (K2-1)*(1-Z1)) / (K1-1) / (K2-1);

    else if (Z1 < Xn1)  // only a gaseous phase appears
        {
            nu = 1;
            Xn1 = Z1;
        }

    else if (Z1 > Xw1)  // onlaa liquid phase appears
        {
            nu = 0;
            Xw1 = Z1;
        }

    nuw = 1 - nu;

    if (isnan(nuw + Xw1 + Xn1))
        DUNE_THROW(RangeError, "Decoupled2p2cni::flashCalculation NANs appear!");
} // end function flashCalculation

template<class G, class RT>
void Decoupled2p2cni<G,RT>::satFlash(double sat, double p, double temp, double poro, double& C1, double& C2, double& Xw1, double& Xn1)
{
    if (sat <= 0 || sat >= 1)
        DUNE_THROW(RangeError,
                   "Decoupled2p2cni :: satFlash cannot handle saturations that equal zero or one!");
    double K1 = problem.liquidPhase.p_vap(temp) / p;
    double K2 = 1. / (p * problem.liquidPhase.henry(temp));
    double xw1 = (1. - K2) / (K1 -K2);
    double xn1 = xw1 * K1;
    Xw1 = xw1 * problem.liquidPhase.molarMass_w()
        / ( xw1 * problem.liquidPhase.molarMass_w() + (1.-xw1) * problem.liquidPhase.molarMass_a() );
    Xn1 = xn1 * problem.liquidPhase.molarMass_w()
        / ( xn1 * problem.liquidPhase.molarMass_w() + (1.-xn1) * problem.liquidPhase.molarMass_a() );
    K1 = Xn1 / Xw1;
    K2 = (1.-Xn1) / (1.-Xw1);

    double rho_w = problem.liquidPhase.density(temp, p, 1.-Xw1);
    double rho_n = problem.gasPhase.density(temp, p, Xn1);

    C1  = poro* (sat * Xw1 * rho_w + (1-sat) * Xn1 * rho_n);
    C2  = poro* (sat * (1. - Xw1) * rho_w + (1-sat) * (1. - Xn1) * rho_n);
}

//! Carries out flash calculation and evaluates volumetric error in each cell.
template<class G, class RT>
void Decoupled2p2cni<G,RT>::postupdate(double t, double dt)
{
    int size = elementmapper.size();
    // iterate through leaf grid an evaluate c0 at cell center
    Iterator eendit = grid.template lend<0>(level_);
    for (Iterator it = grid.template lbegin<0>(level_); it != eendit; ++it)
        {
            int indexi = elementmapper.map(*it);
            // get cell geometry informations
            GeometryType gt = it->geometry().type(); //geometry type
            const FieldVector<ct,dim>& local = ReferenceElements<ct,dim>::general(gt).position(0,0); // cell center in reference element
            FieldVector<ct,dimworld> global = it->geometry().global(local); // cell center in global coordinates

            double poro = problem.soil.porosity(global, *it, local);
            double cp_matrix = problem.soil.heatCap(global, *it, local);
            double C1 = problem.variables.totalConcentration[indexi];
            double C2 = problem.variables.totalConcentration[indexi + size];

            double p = problem.variables.pressure[indexi];
            double nuw;
            up_flash(C1,C2, p, problem.variables.totalConcentration[indexi+2*size][0], cp_matrix, problem.variables.temperature[indexi][0], nuw, problem.variables.wet_X1[indexi][0], problem.variables.nonwet_X1[indexi][0], problem.variables.enthalpy_l[indexi][0], problem.variables.enthalpy_g[indexi][0]);

            double rho_l = problem.liquidPhase.density(problem.variables.temperature[indexi], p, (1. - problem.variables.wet_X1[indexi]));
            double rho_g = problem.gasPhase.density(problem.variables.temperature[indexi], p, problem.variables.nonwet_X1[indexi]);
            problem.variables.saturation[indexi][0] = (nuw / rho_l) / (nuw / rho_l + (1-nuw) / rho_g);
            problem.variables.density_wet[indexi] = rho_l;
            problem.variables.density_nonwet[indexi] = rho_g;

            double massw = (C1 + C2) * nuw;
            double massn = (C1 + C2) * (1-nuw);
            double vol = massw / rho_l + massn / rho_g;
            problem.variables.volErr[indexi] = vol - poro;
        }
    timestep = dt;
}

}//end namespace Dune

#endif /*DECOUPLED2P2C_HH_*/
