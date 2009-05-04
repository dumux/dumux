//$Id$
#ifndef MULTIPHYSICS2P2C_HH
#define MULTIPHYSICS2P2C_HH

// standard library:
#include <float.h>
#include <cmath>
#include <string>
#include <fstream>

// dune environment:
#include <dune/common/helpertemplates.hh>
#include <dune/common/typetraits.hh>
#include <dune/grid/common/mcmgmapper.hh>
#include <dune/grid/io/file/vtk/vtkwriter.hh>

#include <dune/istl/bvector.hh>
#include <dune/common/fvector.hh>
#include <dune/common/stdstreams.hh>
#include <dune/disc/operators/boundaryconditions.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include <dune/istl/preconditioners.hh>

// dumux enviroment:
#include "dumux/transport/fv/numericalflux.hh"
#include "dumux/transport/fv/diffusivepart.hh"
#include "dumux/transport/transportproblem2p2c.hh"
#include "dumux/pardiso/pardiso.hh"

//! @author: Jochen Fritz
// last change: 09.03.2009

namespace Dune
{

/*####################################################*
 *                                                    *
 *     CLASS DECLARATION                              *
 *                                                    *
 *####################################################*/

//! Multiphysics model for two-phase two-component model and single-phase transport model
/**
 * This class implements a compressible two-phase two-component model just as in class Decoupled2p2c connected with
 * an incompressible one-phase transport model.
 * The 2p2c model is applied everywhere, where two phases are present,the 1p transport model everywhere else. The choice
 * of the corresponding subdomains is done automatically.
 */
template<class Grid, class Scalar>
class Multiphysics2p2c
{
    template<int dim>
    struct ElementLayout
    {
        bool contains (GeometryType gt)
        {
            return gt.dim() == dim;
        }
    };

    enum{dim = Grid::dimension};
    enum{dimworld = Grid::dimensionworld};

    // grid typedefs
    typedef typename Grid::LevelGridView GV;
    typedef typename GV::IndexSet IS;
    typedef typename GV::template Codim<0>::Iterator ElementIterator;
    typedef typename Grid::template Codim<0>::EntityPointer EntityPointer;
    typedef typename Grid::LevelGridView::IntersectionIterator IntersectionIterator;
    typedef typename Grid::ctype ct;

    // data typedefs
    typedef FieldMatrix<double,1,1> MB;
    typedef BCRSMatrix<MB> MatrixType;

public:
    typedef BlockVector< FieldVector<Scalar,1> > RepresentationType;

    //! initializes the model and computes initial state
    void initial()
    {
        dinfo << "initialization of the model's variables ..."<<std::endl;
        subdomain = true;
        upd = 0;
        timestep = 0;
        problem.variables.pressure = 3e7;
        initializeMatrix();        dinfo << "matrix initialization"<<std::endl;
        initialguess();                dinfo << "first saturation guess"<<std::endl;
        pressure(true, 0);        dinfo << "first pressure guess"<<std::endl;
        transportInitial();        dinfo << "first guess for mass fractions"<<std::endl;
#if DUNE_MINIMAL_DEBUG_LEVEL >= 3
        vtkout("firstguessmp2p2c",0); // output of first pressure guess and first initial state estimation
#endif
        double dt;
        concentrationUpdate(0, dt, upd);
        upd *= dt;
        pressure(false,0);        dinfo << "final pressure initialization"<<std::endl;
        transportInitial();        dinfo << "final initializiation of saturations and mass fractions"<<std::endl;
        dinfo << "initialization done."<<std::endl;
    }

    //! evaluates the pressure equation and performs one transport time step
    /**
     *  The calculated pressure is stored in the variableclass2p2c object in the problem.
     *  The derivative of the change of the total concentration is given back in the vector updateVec.
     *  @param[in] t current time
     *  @param[out] maximum time step size based on CFL-criterion
     *  @param[out] updateVec derivative of the change of the total concentration
     */
    void update(double t, double& dt, RepresentationType& updateVec)
    {
        upd = 0;
        concentrationUpdate(t, dt, upd);
        upd *= dt;
        timestep = dt;

        pressure(false, t);
        int which = concentrationUpdate(t, dt, updateVec);
        dinfo << "timestep restricting cell: " << which << std::endl;
    }

    //! Post processing step must be called at the end of every timestep
    void postProcessUpdate(double t, double dt);

    //! grid level on which the simulation works
    int level()
    {
        return level_;
    }

    //! gives acces to the totalConcentration vector
    RepresentationType& operator*()
    {
        return problem.variables.totalConcentration;
    }

    // pressure equation functions:
    void initializeMatrix();

    // graphical output
    /**
     * @param name file name of output
     * @param k number of output file
     */
    void vtkout(const char* name, int k)
    {
        problem.variables.volErr *= timestep;

        char fname[128];
        sprintf(fname, "%s-%05d", name, k);

        int globalsize = indexset.size(0);
        RepresentationType C1(globalsize);
        RepresentationType C2(globalsize);
        for (int i = 0; i < globalsize; i++)
        {
            C2[i] = problem.variables.totalConcentration[i + globalsize];
            if (i < globalsize)
                C1[i] = problem.variables.totalConcentration[i];
        }

        Dune::VTKWriter<GV> vtkwriter(grid_.levelView(0));
        vtkwriter.addCellData(problem.variables.saturation, "saturation [-]");
        vtkwriter.addCellData(problem.variables.pressure, "pressure[Pa]");
        vtkwriter.addCellData(C1, "total concentration 1 [kg/m^3]");
        vtkwriter.addCellData(C2, "total concentration 2 [kg/m^3]");
        vtkwriter.addCellData(problem.variables.volErr, "volumetric error [-]");
        vtkwriter.addCellData(problem.variables.wet_X1, "mass fraction 1 in wetting phase [-]");
        vtkwriter.addCellData(problem.variables.nonwet_X1, "mass fraction 1 in nonwetting phase [-]");
        vtkwriter.addCellData(problem.variables.density_wet, "wetting phase density [kg/m^3]");
        vtkwriter.addCellData(problem.variables.density_nonwet, "non-wetting phase density [kg/m^3]");
        vtkwriter.addCellData(problem.variables.mobility_wet, "wetting phase mobility [m*s/kg]");
        vtkwriter.addCellData(problem.variables.mobility_nonwet, "non-wetting phase mobility [m*s/kg]");
        vtkwriter.addCellData(subdomain, "subdomain");
        vtkwriter.write(fname, Dune::VTKOptions::ascii);

        problem.variables.volErr /= timestep;
    }

    //! access to grid object
    const Grid& grid() const
    {
        return grid_;
    }

private:
    //********************************
    // private function declarations
    //********************************
    void initialguess();

    void transportInitial();

    void pressure(bool first, const Scalar t)
    {
        assemble(first, t);
        solve();
        return;
    }

    void assemble(bool first, const Scalar t);

    void volumeDerivatives(FieldVector<ct,dim> globalPos, EntityPointer ep, FieldVector<ct,dim> localPos, double T, double& dV_dm1, double& dV_dm2, double& dV_dp);

    void solve();

    int concentrationUpdate(const Scalar t, Scalar& dt, RepresentationType& updateVec);

    void flashCalculation(double Z1, double p, double temp, double poro, double& sat, double& C1, double& C2, double& Xw1, double& Xn1);

    void satFlash(double sat, double p, double temp, double poro, double& C1, double& C2, double& Xw1, double& Xn1);

    //********************************
    // private variables declarations
    //********************************

    // common variables
    Grid& grid_;
    int level_;
    const IS& indexset;
    double timestep;
    const FieldVector<Scalar,dim> gravity;

    // subdomain map
    BlockVector<FieldVector<bool,1> > subdomain;

    // variables for transport equation:
    TransportProblem2p2c<Grid, Scalar>& problem;
    RepresentationType upd;

    // variables for pressure equation:
    MatrixType A;
    RepresentationType f;
    std::string solverName_;
    std::string preconditionerName_;

public:
    // constructor
    Multiphysics2p2c(
                     Grid& g,
                     TransportProblem2p2c<Grid, Scalar>& prob,
                     int lev = 0,
                     const std::string solverName = "BiCGSTAB",
                     const std::string preconditionerName = "SeqILU0" )
        :grid_(g), level_(lev), indexset(g.levelView(lev).indexSet()), gravity(prob.gravity()), problem(prob),
         A(g.size(lev, 0),g.size(lev, 0), (2*dim+1)*g.size(lev, 0), BCRSMatrix<MB>::random), f(g.size(lev, 0)),
         solverName_(solverName), preconditionerName_(preconditionerName)
    {
        problem.variables.volErr = 0;
        upd.resize(2 * indexset.size(0));
        subdomain.resize(indexset.size(0));
        subdomain = false;
    };

}; //end class declaration



/*####################################################*
 *                                                    *
 *     FUNCTION DEFINITIONS 1: PRESSURE EQUATION      *
 *                                                    *
 *####################################################*/

template<class Grid, class Scalar>
void Multiphysics2p2c<Grid, Scalar>::initializeMatrix()
{
    // determine matrix row sizes
    ElementIterator eItEnd = grid_.template lend<0>(level_);
    for (ElementIterator eIt = grid_.template lbegin<0>(level_); eIt != eItEnd; ++eIt)
    {
        // cell index
        int globalIdxi = indexset.index(*eIt);

        // Initialize row size
        int rowSize = 1;

        // run through all intersections with neighbors
        IntersectionIterator endit = eIt->ilevelend();
        for (IntersectionIterator isIt = eIt->ilevelbegin(); isIt!=endit; ++isIt)
            if (isIt->neighbor()) rowSize++;
        A.setrowsize(globalIdxi, rowSize);
    }
    A.endrowsizes();

    // determine position of matrix entries
    for (ElementIterator eIt = grid_.template lbegin<0>(level_); eIt != eItEnd; ++eIt)
    {
        // cell index
        int globalIdxi = indexset.index(*eIt);

        // add diagonal index
        A.addindex(globalIdxi, globalIdxi);

        // run through all intersections with neighbors
        IntersectionIterator endit = eIt->ilevelend();
        for (IntersectionIterator isIt = eIt->ilevelbegin(); isIt!=endit; ++isIt)
            if (isIt->neighbor())
            {
                // access neighbor
                EntityPointer outside = isIt->outside();
                int globalIdxj = indexset.index(*outside);

                // add off diagonal index
                A.addindex(globalIdxi, globalIdxj);
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
 * This is done to account for the volume effects which appear when gas and liquid are dissolved iin each other.
 */
template<class Grid, class Scalar>
void Multiphysics2p2c<Grid, Scalar>::assemble(bool first, const Scalar t=0)
{
    // initialization: set matrix A to zero
    A = 0;
    f = 0;

    double maxErr = fabs(problem.variables.volErr.infinity_norm());

    int size = indexset.size(0);
    double dV_[size][3];
    for (int i = 0; i < size; i++)
    {
        dV_[i][0] = 0;
    }

    // iterate over all cells in the grid
    ElementIterator eItEnd = grid_.template lend<0>(level_);
    ElementIterator eIt = grid_.template lbegin<0>(level_);
    for (; eIt != eItEnd; ++eIt)
    {
        // get geometry infos about the cell...
        GeometryType gt = eIt->geometry().type(); // cell geometry type
        const FieldVector<ct,dim>&
            localPos = ReferenceElements<ct,dim>::general(gt).position(0,0); // cell center in reference element
        FieldVector<ct,dim> globalPos = eIt->geometry().global(localPos);             // global coordinate of cell center
        double volume = eIt->geometry().integrationElement(localPos)
            *ReferenceElements<ct,dim>::general(gt).volume(); // cell volume

        double Ti = problem.temp(globalPos, *eIt, localPos);

        // cell index in host grid and subdomain flag
        int globalIdxi = indexset.index(*eIt);
        bool subI = subdomain[globalIdxi];

        // get absolute permeability
        FieldMatrix<ct,dim,dim> Ki(this->problem.soil.K(globalPos,*eIt,localPos));

        // get porosity
        double poroI = problem.soil.porosity(globalPos, *eIt, localPos);

        // get the cell's saturation
        double satI = problem.variables.saturation[globalIdxi];

        // derivatives of the fluid volume with respect to mass of compnents and pressure
        double dV_dm1 = 0;
        double dV_dm2 = 0;
        double dV_dp = 0;
        int numberOfFaces = ReferenceElements<ct,dim>::general(gt).size(1);

        // specific volume of the phases
        double Vg = 1. / problem.variables.density_nonwet[globalIdxi];
        double Vw = 1. / problem.variables.density_wet[globalIdxi];

        double lambdaI, fw_I, fn_I;

        if (first)
        {
            // get mobilities and fractional flow factors
            lambdaI = problem.variables.mobility_wet[globalIdxi] + problem.variables.mobility_nonwet[globalIdxi];
            fw_I = problem.variables.mobility_wet[globalIdxi] / lambdaI;
            fn_I = problem.variables.mobility_nonwet[globalIdxi] / lambdaI;

            FieldVector<Scalar,2> q = problem.source(globalPos,*eIt,localPos);
            f[globalIdxi] = volume * (q[0] * Vw + q[1] * Vg);
        }
        else if (!subI)
        {
            // assume that the source does not lead to a second phase!
            FieldVector<Scalar,2> q = problem.source(globalPos,*eIt,localPos);
            if (satI == 1)
                f[globalIdxi] = volume * (q[0] + q[1]) * Vw;
            else if (satI == 0)
                f[globalIdxi] = volume * (q[0] + q[1]) * Vg;
        }
        else
        {
            if (dV_[globalIdxi][0] == 0) volumeDerivatives(globalPos, *eIt, localPos, Ti, dV_[globalIdxi][0], dV_[globalIdxi][1], dV_[globalIdxi][2]);
            dV_dm1 = dV_[globalIdxi][0];
            dV_dm2 = dV_[globalIdxi][1];
            dV_dp = dV_[globalIdxi][2];
            // right hand side entry: sources
            FieldVector<Scalar,2> q = problem.source(globalPos,*eIt,localPos);
            f[globalIdxi] = volume * (dV_dm1 * q[0] + dV_dm2 * q[1]);
        }

        // iterate over all faces of the cell
        IntersectionIterator isItEnd = eIt->ilevelend();
        for (IntersectionIterator isIt = eIt->ilevelbegin(); isIt!=isItEnd; ++isIt)
        {
            // some geometry informations of the face
            GeometryType gtf = isIt->geometryInInside().type(); // get geometry type of face
            const FieldVector<ct,dim-1>& facelocalPos = ReferenceElements<ct,dim-1>::general(gtf).position(0,0); // center in face's reference element
            const FieldVector<ct,dim>& facelocalDim = ReferenceElements<ct,dim>::general(gtf).position(isIt->indexInInside(),1); // center of face inside volume reference element
            FieldVector<ct,dimworld> unitOuterNormal = isIt->unitOuterNormal(facelocalPos);// get normal vector
            FieldVector<ct,dimworld> integrationOuterNormal = isIt->integrationOuterNormal(facelocalPos);
            integrationOuterNormal *= ReferenceElements<ct,dim-1>::general(gtf).volume(); //normal vector scaled with volume
            double faceVol = isIt->geometry().volume(); // get face volume

            // compute directed permeability vector Ki.n
            FieldVector<ct,dim> Kni(0);
            Ki.umv(unitOuterNormal, Kni);

            // handle interior face
            if (isIt->neighbor())
            {
                // acces neighbor
                EntityPointer outside = isIt->outside();

                // cell index and subdomain flag of neighbor
                int globalIdxj = indexset.index(*outside);

                // some geometry infos of the neighbor
                GeometryType nbgt = outside->geometry().type();
                const FieldVector<ct,dim>& nblocalPos = ReferenceElements<ct,dim>::general(nbgt).position(0,0);
                FieldVector<ct,dimworld>
                    nbGlobalPos = outside->geometry().global(nblocalPos); // neighbor cell center in global coordinates

                // distance vector between barycenters
                FieldVector<ct,dimworld> distVec = nbGlobalPos - globalPos;

                // compute distance between cell centers
                double dist = distVec.two_norm();

                // get absolute permeability
                FieldMatrix<ct,dim,dim> Kj(problem.soil.K(nbGlobalPos, *outside, nblocalPos));

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

                // phase densities in cell in neighbor
                double rho_w_I = 1 / Vw;
                double rho_n_I = 1 / Vg;
                double rho_w_J = problem.variables.density_wet[globalIdxj];
                double rho_n_J = problem.variables.density_nonwet[globalIdxj];
                double rho_w = (rho_w_I + rho_w_J) / 2;
                double rho_n = (rho_n_I + rho_n_J) / 2;

                // update diagonal entry
                double entry;
                if (first) // first guess pressure calculation
                {
                    double lambdaJ = problem.variables.mobility_wet[globalIdxj] + problem.variables.mobility_nonwet[globalIdxj];
                    double fw_J = problem.variables.mobility_wet[globalIdxj] / lambdaJ;
                    double fn_J = problem.variables.mobility_nonwet[globalIdxj] / lambdaJ;

                    entry = fabs(0.5 * (lambdaI + lambdaJ) * faceVol * (K * distVec) / (dist * dist));
                    double factor = (fw_I + fw_J) * (rho_w) / 2 + (fn_I + fn_J) * (rho_n) / 2;
                    f[globalIdxi] -= factor * 0.5 * (lambdaI + lambdaJ) * faceVol * (K * gravity);
                }
                else if (!subI) //the current cell is not contained in subgrid
                {
                    double potW = (unitOuterNormal * distVec) * (problem.variables.pressure[globalIdxi] - problem.variables.pressure[globalIdxj]) / (dist * dist);
                    double potN = potW + rho_n * (unitOuterNormal * gravity);
                    potW += rho_w * (unitOuterNormal * gravity);

                    double lambdaW, lambdaN;

                    if (potW >= 0.)
                        lambdaW = problem.variables.mobility_wet[globalIdxi];
                    else
                        lambdaW = problem.variables.mobility_wet[globalIdxj];

                    if (potN >= 0.)
                        lambdaN = problem.variables.mobility_nonwet[globalIdxi];
                    else
                        lambdaN = problem.variables.mobility_nonwet[globalIdxj];

                    entry  = (lambdaW + lambdaN) * faceVol * (K * distVec) / (dist * dist);
                    double rightentry = (rho_w * lambdaW + rho_n * lambdaN) * faceVol * (K * gravity);
                    f[globalIdxi] -= rightentry;
                }
                else // current cell is contained in subgrid
                {
                    double potW = (unitOuterNormal * distVec) * (problem.variables.pressure[globalIdxi] - problem.variables.pressure[globalIdxj]) / (dist * dist);
                    double potN = potW + rho_n * (unitOuterNormal * gravity);
                    potW += rho_w * (unitOuterNormal * gravity);

                    double lambdaW, lambdaN;
                    double dV_w, dV_n;

                    if (subdomain[globalIdxj])
                    {
                        double Tj = problem.temp(nbGlobalPos, *outside, nblocalPos);
                        if (dV_[globalIdxj][0] == 0) volumeDerivatives(globalPos, *outside, localPos, Tj, dV_[globalIdxj][0], dV_[globalIdxj][1], dV_[globalIdxj][2]);
                        dV_dm1 = (dV_[globalIdxi][0] + dV_[globalIdxj][0]) / 2;
                        dV_dm2 = (dV_[globalIdxi][1] + dV_[globalIdxj][1]) / 2;

                        double graddV_dm1 = (dV_[globalIdxj][0] - dV_[globalIdxi][0]) / dist;
                        double graddV_dm2 = (dV_[globalIdxj][1] - dV_[globalIdxi][1]) / dist;

                        double gV_w, gV_n;

                        if (potW >= 0.)
                        {
                            dV_w = rho_w_I * (dV_dm1 * problem.variables.wet_X1[globalIdxi] + dV_dm2 * (1. - problem.variables.wet_X1[globalIdxi]));
                            lambdaW = problem.variables.mobility_wet[globalIdxi];
                            gV_w = rho_w_I * (graddV_dm1 * problem.variables.wet_X1[globalIdxi] + graddV_dm2 * (1. - problem.variables.wet_X1[globalIdxi]));
                        }
                        else
                        {
                            dV_w = rho_w_J * (dV_dm1 * problem.variables.wet_X1[globalIdxj] + dV_dm2 * (1. - problem.variables.wet_X1[globalIdxj]));
                            lambdaW = problem.variables.mobility_wet[globalIdxj];
                            gV_w = rho_w_J * (graddV_dm1 * problem.variables.wet_X1[globalIdxj] + graddV_dm2 * (1. - problem.variables.wet_X1[globalIdxj]));
                        }
                        if (potN >= 0.)
                        {
                            dV_n = rho_n_I * (dV_dm1 * problem.variables.nonwet_X1[globalIdxi] + dV_dm2 * (1. - problem.variables.nonwet_X1[globalIdxi]));
                            lambdaN = problem.variables.mobility_nonwet[globalIdxi];
                            gV_n = rho_n_I * (graddV_dm1 * problem.variables.nonwet_X1[globalIdxi] + graddV_dm2 * (1. - problem.variables.nonwet_X1[globalIdxi]));
                        }
                        else
                        {
                            dV_n = rho_n_J * (dV_dm1 * problem.variables.nonwet_X1[globalIdxj] + dV_dm2 * (1. - problem.variables.nonwet_X1[globalIdxj]));
                            lambdaN = problem.variables.mobility_nonwet[globalIdxj];
                            gV_n = rho_n_J * (graddV_dm1 * problem.variables.nonwet_X1[globalIdxj] + graddV_dm2 * (1. - problem.variables.nonwet_X1[globalIdxj]));
                        }

                        entry = faceVol * (lambdaW * dV_w + lambdaN * dV_n) - volume / numberOfFaces * (lambdaW * gV_w + lambdaN * gV_n);
                        entry *= fabs((K*distVec)/(dist*dist));
                        double rightentry = faceVol * (rho_w * lambdaW * dV_w + rho_n * lambdaN * dV_n) - volume / numberOfFaces * (rho_w * lambdaW * gV_w + rho_n * lambdaN * gV_n);
                        rightentry *= (K * gravity);
                        f[globalIdxi] -= rightentry;
                    }
                    else
                    {
                        dV_dm1 = dV_[globalIdxi][0];
                        dV_dm2 = dV_[globalIdxi][1];

                        if (potW >= 0.)
                        {
                            dV_w = rho_w_I * (dV_dm1 * problem.variables.wet_X1[globalIdxi] + dV_dm2 * (1. - problem.variables.wet_X1[globalIdxi]));
                            lambdaW = problem.variables.mobility_wet[globalIdxi];
                        }
                        else
                        {
                            dV_w = rho_w_J * (dV_dm1 * problem.variables.wet_X1[globalIdxj] + dV_dm2 * (1. - problem.variables.wet_X1[globalIdxj]));
                            lambdaW = problem.variables.mobility_wet[globalIdxj];
                        }
                        if (potN >= 0.)
                        {
                            dV_n = rho_n_I * (dV_dm1 * problem.variables.nonwet_X1[globalIdxi] + dV_dm2 * (1. - problem.variables.nonwet_X1[globalIdxi]));
                            lambdaN = problem.variables.mobility_nonwet[globalIdxi];
                        }
                        else
                        {
                            dV_n = rho_n_J * (dV_dm1 * problem.variables.nonwet_X1[globalIdxj] + dV_dm2 * (1. - problem.variables.nonwet_X1[globalIdxj]));
                            lambdaN = problem.variables.mobility_nonwet[globalIdxj];
                        }

                        entry = lambdaW * dV_w + lambdaN * dV_n;
                        entry *= fabs(faceVol*(K*distVec)/(dist*dist));
                        double rightentry = rho_w * lambdaW * dV_w + rho_n * lambdaN * dV_n;
                        rightentry *= faceVol * (K * gravity);
                        f[globalIdxi] -= rightentry;
                    }
                }

                // set diagonal entry
                A[globalIdxi][globalIdxi] += entry;

                // set off-diagonal entry
                A[globalIdxi][globalIdxj] = -entry;
            }

            // boundary face
            else
            {
                dV_dm1 = dV_[globalIdxi][0];
                dV_dm2 = dV_[globalIdxi][1];

                // center of face in global coordinates
                FieldVector<ct,dimworld> faceGlobalPos = isIt->geometry().global(facelocalPos);

                //get boundary condition for boundary face center
                BoundaryConditions::Flags bctype = problem.press_bc_type(faceGlobalPos, *eIt, facelocalDim);


                if (bctype == BoundaryConditions::dirichlet) //dirichlet boundary
                {
                    double pressBC = problem.dirichlet(faceGlobalPos, *eIt, facelocalDim);

                    // distance vector to boundary face
                    FieldVector<ct,dimworld> distVec(faceGlobalPos - globalPos);
                    double dist = distVec.two_norm();
                    if (first)
                    {
                        A[globalIdxi][globalIdxi] += lambdaI * faceVol * (Kni * distVec) / (dist * dist);
                        f[globalIdxi] += lambdaI * faceVol * pressBC * (Kni * distVec) / (dist * dist);
                        f[globalIdxi] -= (fw_I * problem.variables.density_wet[globalIdxi] + fn_I * problem.variables.density_nonwet[globalIdxi]) * lambdaI * faceVol * (Kni * gravity);
                    }
                    else
                    {
                        double satBound, C1Bound, C2Bound, Xw1Bound, Xn1Bound;
                        double Tb = problem.temp(faceGlobalPos, *eIt, facelocalDim);

                        //get boundary condition type for compositional transport
                        BoundaryConditions2p2c::Flags bctype = problem.bc_type(faceGlobalPos, *eIt, facelocalDim);
                        if (bctype == BoundaryConditions2p2c::saturation) // saturation given
                        {
                            satBound = problem.dirichletSat(faceGlobalPos, *eIt, facelocalDim);
                            satFlash(satBound, pressBC, Tb, problem.soil.porosity(globalPos, *eIt, localPos), C1Bound, C2Bound, Xw1Bound, Xn1Bound);
                        }
                        if (bctype == BoundaryConditions2p2c::concentration) // mass fraction given
                        {
                            double Z1Bound = problem.dirichletConcentration(faceGlobalPos, *eIt, facelocalDim);
                            flashCalculation(Z1Bound, pressBC, Tb, problem.soil.porosity(globalPos, *eIt, localPos), satBound, C1Bound, C2Bound, Xw1Bound, Xn1Bound);
                        }

                        // fluid properties on boundary
                        std::vector<double> kr = problem.materialLaw.kr(satBound, globalPos, *eIt, localPos, Tb);
                        double viscosityW = problem.liquidPhase.viscosity(Tb, pressBC , 1. - Xw1Bound);
                        double viscosityN = problem.gasPhase.viscosity(Tb, pressBC, Xn1Bound);

                        // phase densities in cell and on boundary
                        double rho_w_I = 1 / Vw;
                        double rho_n_I = 1 / Vg;
                        double rho_w_J = problem.liquidPhase.density(Tb, pressBC, 1. - Xw1Bound);
                        double rho_n_J = problem.gasPhase.density(Tb, pressBC, Xn1Bound);
                        double rho_w = (rho_w_I + rho_w_J) / 2;
                        double rho_n = (rho_n_I + rho_n_J) / 2;

                        if (!subI)
                        {
                            // gradient of total potential
                            double potW = (unitOuterNormal * distVec) * (problem.variables.pressure[globalIdxi] - pressBC) / (dist * dist);
                            double potN = potW + rho_n * (unitOuterNormal * gravity);
                            potW += rho_w * (unitOuterNormal * gravity);

                            double lambdaW, lambdaN;

                            if (potW >= 0.)
                                lambdaW = problem.variables.mobility_wet[globalIdxi];
                            else
                                lambdaW = satBound / viscosityW;

                            if (potN >= 0.)
                                lambdaN = problem.variables.mobility_nonwet[globalIdxi];
                            else
                                lambdaN = (1. - satBound) / viscosityN;

                            A[globalIdxi][globalIdxi] += (lambdaW + lambdaN) * faceVol * (Kni * distVec) / (dist * dist);
                            double rightentry = (rho_w * lambdaW + rho_n * lambdaN) * faceVol * (Kni * gravity);
                            f[globalIdxi] -= rightentry;
                            f[globalIdxi] += (lambdaW + lambdaN) * faceVol * pressBC * (Kni * distVec) / (dist * dist);
                        }
                        else
                        {
                            // velocities
                            double potW = (unitOuterNormal * distVec) * (problem.variables.pressure[globalIdxi] - pressBC) / (dist * dist);
                            double potN = potW + rho_n * (unitOuterNormal * gravity);
                            potW += rho_w * (unitOuterNormal * gravity);

                            double lambdaW, lambdaN;
                            double dV_w, dV_n;

                            if (potW >= 0.)
                            {
                                dV_w = rho_w_I * (dV_dm1 * problem.variables.wet_X1[globalIdxi] + dV_dm2 * (1. - problem.variables.wet_X1[globalIdxi]));
                                lambdaW = problem.variables.mobility_wet[globalIdxi];
                            }
                            else
                            {
                                dV_w = rho_w_J * (dV_dm1 * Xw1Bound + dV_dm2 * (1. - Xw1Bound));
                                lambdaW = satBound / viscosityW;
                            }
                            if (potN >= 0.)
                            {
                                dV_n = rho_n_I * (dV_dm1 * problem.variables.nonwet_X1[globalIdxi] + dV_dm2 * (1. - problem.variables.nonwet_X1[globalIdxi]));
                                lambdaN = problem.variables.mobility_nonwet[globalIdxi];
                            }
                            else
                            {
                                dV_n = rho_n_J * (dV_dm1 * Xn1Bound + dV_dm2 * (1. - Xn1Bound));
                                lambdaN = (1. - satBound) / viscosityN;
                            }

                            double entry = lambdaW * dV_w + lambdaN * dV_n;
                            entry *= fabs(faceVol * (Kni * distVec) / (dist * dist));
                            double rightentry = rho_w * lambdaW * dV_w + rho_n *lambdaN * dV_n;
                            rightentry *= faceVol * (Kni * gravity);

                            // set diagonal entry and right hand side entry
                            A[globalIdxi][globalIdxi] += entry;
                            f[globalIdxi] += entry * pressBC;
                            f[globalIdxi] -= rightentry;
                        }
                    }
                }
                else //neumann boundary
                {
                    FieldVector<Scalar,2> J = problem.neumann(faceGlobalPos, *eIt, facelocalDim);
                    if (first)
                        f[globalIdxi] += faceVol * (J[0] * Vw + J[1] * Vg);
                    else if (!subI)
                    {
                        if (satI == 1)
                            f[globalIdxi] += faceVol * (J[0] + J[1]) * Vw;
                        if (satI == 0)
                            f[globalIdxi] += faceVol * (J[0] + J[1]) * Vg;
                    }
                    else
                    {
//                        f[globalIdxi] += faceVol* (J[0] * dV_dm1 + J[1] * dV_dm2);
                    }
                }
            }
        } // end all intersections

        // compressibility term
        if (!first && subI)
        {
            if (timestep != 0.)
            {
                A[globalIdxi][globalIdxi] -= dV_dp / timestep;
                f[globalIdxi] -= problem.variables.pressure[globalIdxi] *dV_dp / timestep;
            }

            // error reduction routine: volumetric error is damped and inserted to right hand side
            // if damping is not done, the solution method gets unstable!
            double erri = fabs(problem.variables.volErr[globalIdxi]);
            double x_lo = 0.6;
            double x_mi = 0.9;
            double fac  = 0.3;
            double lofac = 0.;
            double hifac = 0.;
            hifac /= fac;

            if (erri*timestep > 5e-4)
                if (erri > x_lo * maxErr)
                {
                    if (erri <= x_mi * maxErr)
                        f[globalIdxi] += fac* (1-x_mi*(lofac/fac-1)/(x_lo-x_mi) + (lofac/fac-1)/(x_lo-x_mi)*erri/maxErr)* problem.variables.volErr[globalIdxi]*volume;
                    else
                        f[globalIdxi] += fac * (1 + x_mi - hifac*x_mi/(1-x_mi) + (hifac/(1-x_mi)-1)*erri/maxErr) * problem.variables.volErr[globalIdxi]*volume;
                }
        }
    } // end grid traversal

//                printmatrix(std::cout,A,"stiffnesmatrix","row");
//                printvector(std::cout,f,"right hand side","row");

    return;
} // end function assemble


template<class Grid, class Scalar>
void Multiphysics2p2c<Grid, Scalar>::solve()
{
    typedef FieldVector<double, 1> VB;
    typedef BlockVector<VB> Vector;

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

/*####################################################*
 *                                                    *
 *     FUNCTION DEFINITIONS 2: TRANSPORT EQUATION     *
 *                                                    *
 *####################################################*/

template<class Grid, class Scalar>
void Multiphysics2p2c<Grid,Scalar>::initialguess()
{
    // iterate through leaf grid an evaluate c0 at cell center
    ElementIterator eItEnd = grid_.template lend<0>(level_);
    for (ElementIterator eIt = grid_.template lbegin<0>(level_); eIt != eItEnd; ++eIt)
    {
        int globalIdx = indexset.index(*eIt);

        // get geometry information of cell
        GeometryType gt = eIt->geometry().type();
        const FieldVector<ct,dim>&
            localPos = ReferenceElements<ct,dim>::general(gt).position(0,0);
        FieldVector<ct,dimworld> globalPos = eIt->geometry().global(localPos);

        double T = problem.temp(globalPos, *eIt, localPos);

        // initial conditions
        double sat_0;
        double Z1_0;
        BoundaryConditions2p2c::Flags ictype = problem.initcond_type(globalPos, *eIt, localPos); // get type of initial condition

        if (ictype == BoundaryConditions2p2c::saturation)// saturation initial condition
            sat_0 = problem.initSat(globalPos, *eIt, localPos);
        else if(ictype == BoundaryConditions2p2c::concentration)            // saturation initial condition
        {
            Z1_0 = problem.initConcentration(globalPos, *eIt, localPos);
            double rho_l = problem.liquidPhase.density(T, problem.variables.pressure[globalIdx], 0.);
            sat_0 = Z1_0 / rho_l;
            sat_0 /= Z1_0 / rho_l + (1 - Z1_0) * problem.materialLaw.nonwettingPhase.density(T, problem.variables.pressure[globalIdx], 0.);
        }
        else
        {
            DUNE_THROW(Dune::NotImplemented, "Boundary condition " << ictype);
        }

        // initialize cell saturation
        this->problem.variables.saturation[globalIdx][0] = sat_0;

        // get relative permeabilities and set mobilities.
        // get viscosities
        double viscosityL = problem.liquidPhase.viscosity(T, problem.variables.pressure[globalIdx], 0.);
        double viscosityG = problem.gasPhase.viscosity(T, problem.variables.pressure[globalIdx], 0.);
        std::vector<double> kr = problem.materialLaw.kr(problem.variables.saturation[globalIdx], globalPos, *eIt, localPos, T);
        problem.variables.mobility_wet[globalIdx] = kr[0] / viscosityL;
        problem.variables.mobility_nonwet[globalIdx] = kr[1] / viscosityG;
        problem.variables.density_wet[globalIdx] = problem.liquidPhase.density(T, problem.variables.pressure[globalIdx], 0.);
        problem.variables.density_nonwet[globalIdx] = problem.gasPhase.density(T, problem.variables.pressure[globalIdx], 0.);
    }

    return;
}//end function initialguess


template<class Grid, class Scalar>
void Multiphysics2p2c<Grid,Scalar>::transportInitial()
{
    // iterate through grid an evaluate c0 at cell center
    ElementIterator eItEnd = grid_.template lend<0>(level_);
    for (ElementIterator eIt = grid_.template lbegin<0>(level_); eIt != eItEnd; ++eIt)
    {
        int globalIdx = indexset.index(*eIt);

        // get geometry information of cell
        GeometryType gt = eIt->geometry().type();
        const FieldVector<ct,dim>&
            localPos = Dune::ReferenceElements<ct,dim>::general(gt).position(0,0);
        FieldVector<ct,dimworld> globalPos = eIt->geometry().global(localPos);

        double T = problem.temp(globalPos, *eIt, localPos);

        // initial conditions
        double sat_0, C1_0, C2_0, viscosityL, viscosityG;
        Dune::BoundaryConditions2p2c::Flags ictype = problem.initcond_type(globalPos, *eIt, localPos);            // get type of initial condition

        if (ictype == Dune::BoundaryConditions2p2c::saturation)  // saturation initial condition
        {
            sat_0 = problem.initSat(globalPos, *eIt, localPos);
            satFlash(sat_0, problem.variables.pressure[globalIdx], T, problem.soil.porosity(globalPos, *eIt, localPos), C1_0, C2_0, problem.variables.wet_X1[globalIdx][0], problem.variables.nonwet_X1[globalIdx][0]);
        }
        else if (ictype == Dune::BoundaryConditions2p2c::concentration) // concentration initial condition
        {
            double Z1_0 = problem.initConcentration(globalPos, *eIt, localPos);
            flashCalculation(Z1_0, problem.variables.pressure[globalIdx], T, problem.soil.porosity(globalPos, *eIt, localPos), sat_0, C1_0, C2_0, problem.variables.wet_X1[globalIdx][0], problem.variables.nonwet_X1[globalIdx][0]);
        }
        // initialize cell concentration
        problem.variables.totalConcentration[globalIdx] = C1_0;
        problem.variables.totalConcentration[globalIdx + indexset.size(0)] = C2_0;
        problem.variables.saturation[globalIdx][0] = sat_0;
        problem.variables.density_wet[globalIdx] = problem.liquidPhase.density(T, problem.variables.pressure[globalIdx], 1. - problem.variables.wet_X1[globalIdx]);
        problem.variables.density_nonwet[globalIdx] = problem.gasPhase.density(T, problem.variables.pressure[globalIdx], problem.variables.nonwet_X1[globalIdx][0]);

        // get viscosities
        viscosityL = problem.liquidPhase.viscosity(T, problem.variables.pressure[globalIdx], 1. - problem.variables.wet_X1[globalIdx]);
        viscosityG = problem.gasPhase.viscosity(T, problem.variables.pressure[globalIdx], problem.variables.nonwet_X1[globalIdx]);

        // get relative permeabilities and set mobilities.
        std::vector<double> kr = problem.materialLaw.kr(problem.variables.saturation[globalIdx], globalPos, *eIt, localPos, T);
        problem.variables.mobility_wet[globalIdx] = kr[0] / viscosityL;
        problem.variables.mobility_nonwet[globalIdx] = kr[1] / viscosityG;
    }
    problem.variables.volErr = 0;
    return;
} //end function transportInitial

template<class Grid, class Scalar>
int Multiphysics2p2c<Grid,Scalar>::concentrationUpdate(const Scalar t, Scalar& dt, RepresentationType& updateVec)
{
    // initialize timestep dt very large
    dt = 1E100;

    // set update vector to zero
    updateVec = 0;

    int which;

    // compute update vector
    ElementIterator eItEnd = grid_.template lend<0>(level_);
    for (ElementIterator eIt = grid_.template lbegin<0>(level_); eIt != eItEnd; ++eIt)
    {
        // get cell geometry informations
        GeometryType gt = eIt->geometry().type(); //geometry type
        const FieldVector<ct,dim>& localPos = ReferenceElements<ct,dim>::general(gt).position(0,0); // cell center in reference element
        FieldVector<ct,dimworld> globalPos = eIt->geometry().global(localPos); // cell center in global coordinates
        double volume = eIt->geometry().integrationElement(localPos) * ReferenceElements<ct,dim>::general(gt).volume(); // cell volume, assume linear map here

        // cell index
        int globalIdxi = indexset.index(*eIt);

        double pressi = problem.variables.pressure[globalIdxi];

        // get saturation and concentration value at cell center
        double satI = problem.variables.saturation[globalIdxi];
        double Xw1_I = problem.variables.wet_X1[globalIdxi];
        double Xn1_I = problem.variables.nonwet_X1[globalIdxi];

        // phase densities in cell
        double rho_w_I = problem.variables.density_wet[globalIdxi];
        double rho_n_I = problem.variables.density_nonwet[globalIdxi];

        // some variables for time step calculation
        double sumfactorin = 0;
        double sumfactorout = 0;

        // get absolute permeability
        FieldMatrix<ct,dim,dim> Ki(problem.soil.K(globalPos,*eIt,localPos));

        // run through all intersections with neighbors and boundary
        IntersectionIterator isItEnd = eIt->ilevelend();
        for (IntersectionIterator isIt = eIt->ilevelbegin(); isIt!=isItEnd; ++isIt)
        {
            // get geometry informations of face
            GeometryType gtf = isIt->geometryInInside().type(); //geometry type
            const FieldVector<ct,dim-1>& faceLocalPos = ReferenceElements<ct,dim-1>::general(gtf).position(0,0); // center in face's reference element
            const FieldVector<ct,dim>& facelocalDim = ReferenceElements<ct,dim>::general(gtf).position(isIt->indexInInside(),1); // center of face inside volume reference element
            double faceVol = isIt->geometry().volume(); // get face volume
            FieldVector<ct,dimworld> unitOuterNormal = isIt->unitOuterNormal(faceLocalPos);

            // get normal vector scaled with volume of face
            FieldVector<ct,dimworld> integrationOuterNormal = isIt->integrationOuterNormal(faceLocalPos);
            integrationOuterNormal *= ReferenceElements<ct,dim-1>::general(gtf).volume();

            // variables for timestep calculation
            double factor[2], factorC1, factorC2;

            if (isIt->neighbor()) // handle interior face
            {
                // access neighbor
                EntityPointer outside = isIt->outside();
                int globalIdxj = indexset.index(*outside);

                // neighbor geometry informations
                GeometryType nbgt = outside->geometry().type();
                const FieldVector<ct,dim>& nbLocalPos = ReferenceElements<ct,dim>::general(nbgt).position(0,0);
                FieldVector<ct,dimworld> nbGlobalPos = outside->geometry().global(nbLocalPos); // neighbor cell center in global coordinates

                // distance vector between barycenters
                FieldVector<ct,dimworld> distVec = nbGlobalPos - globalPos;

                // distance between barycenters
                double dist = distVec.two_norm();

                double pressj = problem.variables.pressure[globalIdxj];

                // get saturation and concentration value at neighbor cell center
                double Xw1_J = problem.variables.wet_X1[globalIdxj];
                double Xn1_J = problem.variables.nonwet_X1[globalIdxj];

                // phase densities in neighbor
                double rho_w_J = problem.variables.density_wet[globalIdxj];
                double rho_n_J = problem.variables.density_nonwet[globalIdxj];

                // get absolute permeability
                FieldMatrix<ct,dim,dim> Kj(problem.soil.K(nbGlobalPos, *outside, nbLocalPos));

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

                double rho_w = (rho_w_I + rho_w_J) / 2;
                double rho_n = (rho_n_I + rho_n_J) / 2;

                // velocities
                double potW = (unitOuterNormal * distVec) * (pressi - pressj) / (dist * dist);
                double potN = potW + rho_n * (unitOuterNormal * gravity);
                potW += rho_w * (unitOuterNormal * gravity);
                //                  lambda = 2 * lambdaI * lambdaJ / (lambdaI + lambdaJ);

                double lambdaW, lambdaN;

                if (potW >= 0.)
                    lambdaW = problem.variables.mobility_wet[globalIdxi];
                else
                    lambdaW = problem.variables.mobility_wet[globalIdxj];

                if (potN >= 0.)
                    lambdaN = problem.variables.mobility_nonwet[globalIdxi];
                else
                    lambdaN = problem.variables.mobility_nonwet[globalIdxj];

                double velocityW = lambdaW * ((K * distVec) * (pressi - pressj) / (dist * dist) + (K * gravity) * rho_w);
                double velocityN = lambdaN * ((K * distVec) * (pressi - pressj) / (dist * dist) + (K * gravity) * rho_n);

                // standardized velocity
                double velocityJIw = std::max(-velocityW * faceVol / volume, 0.0);
                double velocityIJw = std::max( velocityW * faceVol / volume, 0.0);
                double velocityJIn = std::max(-velocityN * faceVol / volume, 0.0);
                double velocityIJn = std::max( velocityN * faceVol / volume, 0.0);

                // for timestep control
                factor[0] = velocityJIw + velocityJIn;
                double foutw = velocityIJw / (satI - problem.soil.Sr_w(globalPos, *eIt, localPos));
                double foutn = velocityIJn / (1 - satI - problem.soil.Sr_n(globalPos, *eIt, localPos));
                if (isnan(foutw) || isinf(foutw) || foutw < 0) foutw = 0;
                if (isnan(foutn) || isinf(foutn) || foutn < 0) foutn = 0;
                factor[1] = foutw + foutn;

                factorC1 =
                    velocityJIw * Xw1_J * rho_w_J
                    - velocityIJw * Xw1_I * rho_w_I
                    + velocityJIn * Xn1_J * rho_n_J
                    - velocityIJn * Xn1_I * rho_n_I;
                factorC2 =
                    velocityJIw * (1. - Xw1_J) * rho_w_J
                    - velocityIJw * (1. - Xw1_I) * rho_w_I
                    + velocityJIn * (1. - Xn1_J) * rho_n_J
                    - velocityIJn * (1. - Xn1_I) * rho_n_I;
            }

            else // handle boundary face
            {
                // face center in globel coordinates
                FieldVector<ct,dim> faceGlobalPos = isIt->geometry().global(faceLocalPos);

                // distance vector between cell and face center
                FieldVector<ct,dimworld> distVec = faceGlobalPos - globalPos;
                double dist = distVec.two_norm();

                // compute directed permeability vector Ki.n
                FieldVector<ct,dim> Kni(0);
                Ki.umv(unitOuterNormal, Kni);

                //get boundary conditions
                BoundaryConditions::Flags pressBCtype = problem.press_bc_type(faceGlobalPos, *eIt, facelocalDim);
                if (pressBCtype == BoundaryConditions::dirichlet)
                {
                    double pressBound = problem.dirichlet(faceGlobalPos, *eIt, facelocalDim);

                    double Tb = problem.temp(faceGlobalPos, *eIt, facelocalDim);

                    double satBound, C1Bound, C2Bound, Xw1Bound, Xn1Bound;
                    BoundaryConditions2p2c::Flags bctype = problem.bc_type(faceGlobalPos, *eIt, facelocalDim);
                    if (bctype == BoundaryConditions2p2c::saturation)
                    {
                        satBound = problem.dirichletSat(faceGlobalPos, *eIt, facelocalDim);
                        satFlash(satBound, pressBound, Tb, problem.soil.porosity(globalPos, *eIt, localPos), C1Bound, C2Bound, Xw1Bound, Xn1Bound);
                    }
                    if (bctype == BoundaryConditions2p2c::concentration)
                    {
                        double Z1Bound = problem.dirichletConcentration(faceGlobalPos, *eIt, facelocalDim);
                        flashCalculation(Z1Bound, pressBound, Tb, problem.soil.porosity(globalPos, *eIt, localPos), satBound, C1Bound, C2Bound, Xw1Bound, Xn1Bound);
                    }

                    // phase densities on boundary
                    double rho_w_Bound = problem.liquidPhase.density(Tb, pressBound, 1. - Xw1Bound);
                    double rho_n_Bound = problem.gasPhase.density(Tb, pressBound, Xn1Bound);

                    // neighbor cell
                    double viscosityW = problem.liquidPhase.viscosity(Tb, pressBound , 1. - Xw1Bound);
                    double viscosityN = problem.gasPhase.viscosity(Tb, pressBound, Xn1Bound);

                    double rho_w = (rho_w_I + rho_w_Bound) / 2;
                    double rho_n = (rho_n_I + rho_n_Bound) / 2;

                    // velocities
                    double potW = (unitOuterNormal * distVec) * (pressi - pressBound) / (dist * dist);
                    double potN = potW + rho_n * (unitOuterNormal * gravity);
                    potW += rho_w * (unitOuterNormal * gravity);

                    double lambdaW, lambdaN;

                    if (potW >= 0.)
                        lambdaW = problem.variables.mobility_wet[globalIdxi];
                    else
                        lambdaW = satBound / viscosityW;

                    if (potN >= 0.)
                        lambdaN = problem.variables.mobility_nonwet[globalIdxi];
                    else
                        lambdaN = (1. - satBound) / viscosityN;

                    double velocityW = lambdaW * ((Kni * distVec) * (pressi - pressBound) / (dist * dist) + (Kni * gravity) * rho_w);
                    double velocityN = lambdaN * ((Kni * distVec) * (pressi - pressBound) / (dist * dist) + (Kni * gravity) * rho_n);

                    // standardized velocity
                    double velocityJIw = std::max(-velocityW * faceVol / volume, 0.0);
                    double velocityIJw = std::max( velocityW * faceVol / volume, 0.0);
                    double velocityJIn = std::max(-velocityN * faceVol / volume, 0.0);
                    double velocityIJn = std::max( velocityN* faceVol / volume, 0.0);

                    // for timestep control
                    factor[0] = velocityJIw + velocityJIn;
                    double foutw = velocityIJw / (satI - problem.soil.Sr_w(globalPos, *eIt, localPos));
                    double foutn = velocityIJn / (1 - satI - problem.soil.Sr_n(globalPos, *eIt, localPos));
                    if (isnan(foutw) || isinf(foutw) || foutw < 0) foutw = 0;
                    if (isnan(foutn) || isinf(foutn) || foutn < 0) foutn = 0;
                    factor[1] = foutw + foutn;

                    factorC1 =
                        + velocityJIw * Xw1Bound * rho_w_Bound
                        - velocityIJw * Xw1_I * rho_w_I
                        + velocityJIn * Xn1Bound * rho_n_Bound
                        - velocityIJn * Xn1_I * rho_n_I ;
                    factorC2 =
                        velocityJIw * (1. - Xw1Bound) * rho_w_Bound
                        - velocityIJw * (1. - Xw1_I) * rho_w_I
                        + velocityJIn * (1. - Xn1Bound) * rho_n_Bound
                        - velocityIJn * (1. - Xn1_I) * rho_n_I ;
                }
                else if (pressBCtype == BoundaryConditions::neumann)
                {
                    FieldVector<Scalar,2> J = problem.neumann(faceGlobalPos, *eIt, facelocalDim);
                    double faceVol = integrationOuterNormal.two_norm();
                    factorC1 = J[0] * faceVol / volume;
                    factorC2 = J[1] * faceVol / volume;

                    // for timestep control
                    factor[0] = 0;
                    factor[1] = 0;
                }

                else DUNE_THROW(NotImplemented, "there is no process boundary condition implemented");
            }

            // add to update vector
            updateVec[globalIdxi] += factorC1;
            updateVec[indexset.size(0)+globalIdxi] += factorC2;

            // for time step calculation
            sumfactorin += factor[0];
            sumfactorout += factor[1];

        } // end all intersections

        // handle source term
        FieldVector<double,2> q = problem.source(globalPos, *eIt, localPos);
        updateVec[globalIdxi] += q[0];
        updateVec[globalIdxi +  indexset.size(0)] += q[1];

        // account for porosity
        double poro = problem.soil.porosity(globalPos, *eIt, localPos);

        sumfactorin = std::max(sumfactorin,sumfactorout) / poro;

        if ( 1./sumfactorin < dt)
        {
            dt = 1./sumfactorin;
            which= globalIdxi;
        }

    } // end grid traversal
    return which;
} // end function "concentrationUpdate"

template<class Grid, class Scalar>
void Multiphysics2p2c<Grid,Scalar>::flashCalculation(double Z1, double p, double temp, double poro, double& sat, double& C1, double& C2, double& Xw1, double& Xn1)
{
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

    double nu2 = 0;

    if (Z1 > Xn1 && Z1 < Xw1)
        nu2 = -((K1-1)*Z1 + (K2-1)*(1-Z1)) / (K1-1) / (K2-1);

    else if (Z1 < Xn1)
    {
        nu2 = 1;
        Xn1 = Z1;
    }
    else if (Z1 > Xw1)
    {
        nu2 = 0;
        Xw1 = Z1;
    }

    double rho_w = problem.liquidPhase.density(temp, p, 1. - Xw1);
    double rho_n = problem.gasPhase.density(temp, p, Xn1);

    sat = (1-nu2) / rho_w;
    sat /= ((1-nu2)/rho_w + nu2/rho_n);

    C1 = poro * (Xw1 * sat * rho_w + Xn1 * (1-sat) * rho_n);
    C2 = poro * ((1. - Xw1) * sat * rho_w + (1. - Xn1) * (1-sat) * rho_n);

} // end function flashCalculation

template<class Grid, class Scalar>
void Multiphysics2p2c<Grid,Scalar>::satFlash(double sat, double p, double temp, double poro, double& C1, double& C2, double& Xw1, double& Xn1)
{
    if (sat <= 0 || sat >= 1)
        DUNE_THROW(RangeError,
                   "Multiphysics2p2c :: saturation initial and boundary conditions may not equal zero or one!");
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

template<class Grid, class Scalar>
void Multiphysics2p2c<Grid,Scalar>::volumeDerivatives(FieldVector<ct,dim> globalPos, EntityPointer ep, FieldVector<ct,dim> localPos, double T, double& dV_dm1, double& dV_dm2, double& dV_dp)
{
    int globalIdxi = indexset.index(*ep);

    double poroI = problem.soil.porosity(globalPos, *ep, localPos);

    GeometryType gt = ep->geometry().type();
    double volume = ep->geometry().integrationElement(localPos)
        *ReferenceElements<ct,dim>::general(gt).volume();

    double Vw = 1. / problem.variables.density_wet[globalIdxi];
    double Vg = 1. / problem.variables.density_nonwet[globalIdxi];
    double sati = problem.variables.saturation[globalIdxi];
    // mass of components inside the cell
    double m1 = problem.variables.totalConcentration[globalIdxi] * volume;
    double m2 = problem.variables.totalConcentration[globalIdxi+indexset.size(0)] * volume;
    // mass fraction of wetting phase
    double nuw1 =  sati/ Vw / (sati/Vw + (1-sati)/Vg);
    // actual fluid volume
    double volalt = (m1+m2) * (nuw1 * Vw + (1-nuw1) * Vg);

    // increments for numerical derivatives
    double inc1 = (fabs(upd[globalIdxi][0]) > 1e-8 / Vw) ?  upd[globalIdxi][0] : 1e-8/Vw;
    double inc2 =(fabs(upd[globalIdxi+indexset.size(0)][0]) > 1e-8 / Vg) ?  upd[globalIdxi+indexset.size(0)][0] : 1e-8 / Vg;
    inc1 *= volume;
    inc2 *= volume;

    // numerical derivative of fluid volume with respect to mass of component 1
    m1 +=  inc1;
    double Z1 = m1 / (m1 + m2);
    double dummy1, dummy2, dummy3, dummy4, satt;
    flashCalculation(Z1, problem.variables.pressure[globalIdxi], T, poroI, satt, dummy1, dummy2, dummy3, dummy4);
    double nuw = satt / Vw / (satt/Vw + (1-satt)/Vg);
    dV_dm1 = ((m1+m2) * (nuw * Vw + (1-nuw) * Vg) - volalt) /inc1;
    m1 -= inc1;

    // numerical derivative of fluid volume with respect to mass of component 2
    m2 += inc2;
    Z1 = m1 / (m1 + m2);
    flashCalculation(Z1, problem.variables.pressure[globalIdxi], T, poroI, satt, dummy1, dummy2, dummy3, dummy4);
    nuw = satt / Vw / (satt/Vw + (1-satt)/Vg);
    dV_dm2 = ((m1+m2) * (nuw * Vw + (1-nuw) * Vg) - volalt)/ inc2;
    m2 -= inc2;

    // numerical derivative of fluid volume with respect to pressure
    double incp = 1e-2;
    double p_ = problem.variables.pressure[globalIdxi] + incp;
    double Vg_ = 1. / problem.gasPhase.density(T, p_, problem.variables.nonwet_X1[globalIdxi]);
    double Vw_ = 1. / problem.liquidPhase.density(T, p_, 1. - problem.variables.wet_X1[globalIdxi]);
    dV_dp = ((m1+m2) * (nuw1 * Vw_ + (1-nuw1) * Vg_) - volalt) /incp;
}

template<class Grid, class Scalar>
void Multiphysics2p2c<Grid,Scalar>::postProcessUpdate(double t, double dt)
{
    problem.variables.volErr = 0.;
    int size = indexset.size(0);
    BlockVector<FieldVector<bool,1> > nextsubdomain;
    nextsubdomain.resize(size);
    nextsubdomain = false;
    double viscosityL, viscosityG;

    // iterate through leaf grid an evaluate c0 at cell center
    ElementIterator eItEnd = grid_.template lend<0>(level_);
    for (ElementIterator eIt = grid_.template lbegin<0>(level_); eIt != eItEnd; ++eIt)
    {
        int globalIdxi = indexset.index(*eIt);
        // get cell geometry informations
        GeometryType gt = eIt->geometry().type(); //geometry type
        const FieldVector<ct,dim>& localPos = ReferenceElements<ct,dim>::general(gt).position(0,0); // cell center in reference element
        FieldVector<ct,dimworld> globalPos = eIt->geometry().global(localPos); // cell center in global coordinates

        double T = problem.temp(globalPos, *eIt, localPos);

        double poro = problem.soil.porosity(globalPos, *eIt, localPos);
        double C1 = problem.variables.totalConcentration[globalIdxi][0];
        double C2 = problem.variables.totalConcentration[indexset.size(0)+globalIdxi][0];
        double Z1 = C1 / (C1 + C2);

        if (subdomain[globalIdxi])
        {
            // get cell geometry informations
            GeometryType gt = eIt->geometry().type(); //geometry type
            const FieldVector<ct,dim>& localPos = ReferenceElements<ct,dim>::general(gt).position(0,0); // cell center in reference element
            FieldVector<ct,dimworld> globalPos = eIt->geometry().global(localPos); // cell center in global coordinates

            flashCalculation(Z1, problem.variables.pressure[globalIdxi], T, poro, problem.variables.saturation[globalIdxi][0], C1, C2, problem.variables.wet_X1[globalIdxi][0], problem.variables.nonwet_X1[globalIdxi][0]);

            double rho_l = problem.liquidPhase.density(T, problem.variables.pressure[globalIdxi][0],1. - problem.variables.wet_X1[globalIdxi][0]);
            double rho_g = problem.gasPhase.density(T, problem.variables.pressure[globalIdxi][0],problem.variables.nonwet_X1[globalIdxi][0]);
            problem.variables.density_wet[globalIdxi][0] = rho_l;
            problem.variables.density_nonwet[globalIdxi][0] = rho_g;

            // get viscosities
            viscosityL = problem.liquidPhase.viscosity(T, problem.variables.pressure[globalIdxi], 1. - problem.variables.wet_X1[globalIdxi][0]);
            viscosityG = problem.gasPhase.viscosity(T, problem.variables.pressure[globalIdxi], problem.variables.nonwet_X1[globalIdxi][0]);

            double nuw = problem.variables.saturation[globalIdxi] * rho_l / (problem.variables.saturation[globalIdxi] * rho_l + (1-problem.variables.saturation[globalIdxi]) * rho_g);
            double massw = (problem.variables.totalConcentration[globalIdxi][0] + problem.variables.totalConcentration[indexset.size(0)+globalIdxi][0]) * nuw;
            double massn = (problem.variables.totalConcentration[globalIdxi][0] + problem.variables.totalConcentration[indexset.size(0)+globalIdxi][0]) * (1-nuw);
            double vol = massw / rho_l + massn / rho_g;
            problem.variables.volErr[globalIdxi] = (vol - poro) / dt;

            if (problem.variables.saturation[globalIdxi] != 1. ||fabs(problem.variables.volErr[globalIdxi] * dt) > 5e-4)
            {
                nextsubdomain[globalIdxi] = true;
                IntersectionIterator isItEnd = eIt->ilevelend();
                for (IntersectionIterator isIt = eIt->ilevelbegin();    isIt!=isItEnd; ++isIt)
                {
                    if (isIt->neighbor())
                    {
                        int globalIdxj = indexset.index(*(isIt->outside()));
                        if (!nextsubdomain[globalIdxj] && !subdomain[globalIdxj])
                        {
                            // get cell geometry informations
                            GeometryType gt = eIt->geometry().type(); //geometry type
                            const FieldVector<ct,dim>& localPos = ReferenceElements<ct,dim>::general(gt).position(0,0); // cell center in reference element
                            FieldVector<ct,dimworld> globalPos = eIt->geometry().global(localPos); // cell center in global coordinates
                            double rho_lj = problem.liquidPhase.density(T, problem.variables.pressure[globalIdxi][0], 0.);
                            double poroJ = problem.soil.porosity(globalPos, *(isIt->outside()), localPos);
                            problem.variables.wet_X1[globalIdxj] = 1. - problem.variables.totalConcentration[globalIdxj + size] / poroJ / rho_lj;
                            problem.variables.totalConcentration[globalIdxj] = problem.variables.wet_X1[globalIdxj] * poroJ * rho_lj;
                            problem.variables.saturation[globalIdxj] = 1.;
                            problem.variables.nonwet_X1[globalIdxj] = 0.;
                        }
                        nextsubdomain[globalIdxj] = true;
                    }
                }
            }
        }
        else
        {
            if (problem.variables.saturation[globalIdxi] == 1.)
            {
                problem.variables.wet_X1[globalIdxi] = Z1;
                problem.variables.nonwet_X1[globalIdxi] = 0.;
            }
            else
            {
                problem.variables.wet_X1[globalIdxi] = 0;
                problem.variables.nonwet_X1[globalIdxi] = 1. - Z1;
            }

            // get viscosities
            viscosityL = problem.liquidPhase.viscosity(T, problem.variables.pressure[globalIdxi], 0.);
            viscosityG = problem.gasPhase.viscosity(T, problem.variables.pressure[globalIdxi], 0.);

            problem.variables.density_wet[globalIdxi] = problem.liquidPhase.density(T, problem.variables.pressure[globalIdxi], 0.);
            problem.variables.density_nonwet[globalIdxi] = problem.gasPhase.density(T, problem.variables.pressure[globalIdxi], 0.);
        }

        // get relative permeabilities and set mobilities.
        std::vector<double> kr = problem.materialLaw.kr(problem.variables.saturation[globalIdxi], globalPos, *eIt, localPos, T);
        problem.variables.mobility_wet[globalIdxi] = kr[0] / viscosityL;
        problem.variables.mobility_nonwet[globalIdxi] = kr[1] / viscosityG;

    }// end grid traversal

    subdomain = nextsubdomain;

    timestep = dt;
}

}//end namespace Dune

#endif /*MULTIPHYSICS2P2C_HH_*/
