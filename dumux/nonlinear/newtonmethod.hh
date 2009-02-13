// $Id$

#ifndef DUNE_NEWTONMETHOD_HH
#define DUNE_NEWTONMETHOD_HH

namespace Dune {

/** \todo Please doc me! */

template<class G, class Model> class NewtonMethod {
    typedef typename Model::FunctionType FunctionType;
    typedef typename Model::OperatorAssembler OperatorAssembler;
public:
    void execute(bool verbose = true) {
        if (verbose)
            std::cout.setf(std::ios_base::scientific, std::ios_base::floatfield);

        double uNorm = std::max((*u).two_norm(), 1.0);
        grid.comm().sum(&uNorm, 1);
        double oneByMagnitude = 1.0/uNorm;
        double relDiff = 1e100;
        double globalResiduum = 1e100;
        bool divided = false;
        double dt = model.getDt();
//        char buf[128];

        while (dt > minDt && (relDiff > diffTolerance || globalResiduum > resTolerance)) {
            relDiff = 1e100;
            int iter = 0;
            int relDiffIncreased = 0;
            double lambda = 1.0;
            double lambdaOld = lambda;
            model.globalDefect(defectGlobal);
            globalResiduum = 0.5*(*defectGlobal).two_norm();
            grid.comm().sum(&globalResiduum, 1);
//            globalResiduum = model.residual(defectGlobal);
            double residuumWeight = 1.0/std::max(globalResiduum, 1.0e-8);
            globalResiduum *= residuumWeight;
            double globalResiduumOld = globalResiduum;
            if (grid.comm().rank() == 0)
                std::cout << "initial residual = " << globalResiduum << std::endl;
            //printvector(std::cout, *defectGlobal, "global Defect", "row", 200, 1, 3);
            while ((relDiff > diffTolerance || globalResiduum > resTolerance)
                    && relDiffIncreased < maxIncreased
                    && iter < maxIter)
            {
                iter ++;
                double relDiffOld = relDiff;
                //int iiter = 0;
                num = 0;
                *uOldNewtonStep = *u;
                globalResiduumOld = globalResiduum;
                lambda = lambdaOld;
                //printvector(std::cout, *uOldNewtonStep, "uOldNewtonStep", "row", 200, 1, 3);
                //printvector(std::cout, *(model.uOldTimeStep), "uOldTimeStep", "row", 200, 1, 3);
                model.assemble();
//                *(uOldNewtonStep) = 0.0;
//                *(model.uOldTimeStep) = 1.0;
//                (*A).mv(*(model.uOldTimeStep), *uOldNewtonStep);
//                sprintf(buf, "rank %d, A*1: ", grid.comm().rank());
//                printvector(std::cout, *uOldNewtonStep, buf, "row", 200, 1, 3);

                if (grid.comm().rank() == 1) {
                    printmatrix(std::cout, *A, "global stiffness matrix", "row", 11, 4);
                    printvector(std::cout, *uOldNewtonStep, "uOldNewtonStep", "row", 3, 1, 3);
                    printvector(std::cout, *f, "right hand side", "row", 3, 1, 3);
                }

                model.solve();
                relDiff = oneByMagnitude*((*u).two_norm());
//                printvector(std::cout, *u, "update", "row", 3, 1, 3);
                *u *= -lambda; // hm, lambda is always 1.0, right???
                *u += *uOldNewtonStep;
//                sprintf(buf, "rank %d, solution: ", grid.comm().rank());
//                printvector(std::cout, *u, buf, "row", 3, 1, 3);
                model.globalDefect(defectGlobal);
                globalResiduum = residuumWeight*0.5*(*defectGlobal).two_norm();
                grid.comm().sum(&globalResiduum, 1);
//                globalResiduum = residuumWeight*model.residual(defectGlobal);
//                sprintf(buf, "rank %d, global Defect: ", grid.comm().rank());
//                printvector(std::cout, *defectGlobal, buf, "row", 200, 1, 3);

                if (verbose && grid.comm().rank() == 0)
                    std::cout << "Newton step "<< iter << ", residual = "
                    << globalResiduum << ", difference = "<< relDiff << std::endl;

                if (relDiff > relDiffOld)
                    relDiffIncreased++;
            }
            if (relDiff > diffTolerance || globalResiduum > resTolerance) {
                dt *= 0.5;
                if (grid.comm().rank() == 0) {
                    std::cout << "NewtonMethod::execute(), tolerances = "
                    << diffTolerance << " , " << resTolerance;
                    if (relDiffIncreased > 1)
                        std::cout << ": relative difference increased " << relDiffIncreased
                            << " time(s)."<< std::endl;
                    else
                        std::cout << ": did not converge in "<< iter << " iterations."<< std::endl;
                    std::cout << "Retry same time step with reduced size of " << dt << std::endl;
                }
                model.setDt(dt);
                *u = *(model.uOldTimeStep);
                divided = true;
            }
            else {
                if (grid.comm().rank() == 0)
                    std::cout << "Converged. Residual = "<< globalResiduum
                    << ", difference = "<< relDiff << std::endl;
                if ((!divided) && iter < goodIter) {
                    dt *= 2;
                    if (grid.comm().rank() == 0)
                        std::cout << "Below "<< goodIter
                        << " Newton iterations. Initial size for the next time step doubled to "
                        << dt << std::endl;
                }
                model.setDt(dt);

                return;
            }
        }

        if (dt <= minDt)
            DUNE_THROW(MathError, "NewtonMethod:: time step size below minimum " << minDt << ".");
        return;
    }

    NewtonMethod(const G& g, Model& mod, double dtol = 1e-4,
            double rtol = 1e-7, int maxIt = 10, double mindt = 1e-5,
            int goodIt = 4, int maxInc = 2)
    : grid(g), model(mod), u(mod.u), f(mod.f), A(mod.A),
    uOldNewtonStep(g, g.overlapSize(0)==0),
    diffTolerance(dtol), resTolerance(rtol), maxIter(maxIt),
    defectGlobal(g, g.overlapSize(0)==0), num(0), minDt(mindt), goodIter(goodIt), maxIncreased(maxInc)
    {}

    NewtonMethod(const G& g, Model& mod, int level, double dtol = 1e-8,
            double rtol = 1e-5, int maxIt = 12, double mindt = 1e-5,
            int goodIt = 3, int maxInc = 2)
    : grid(g), model(mod), u(mod.u), f(mod.f), A(mod.A),
    uOldNewtonStep(g, level, g.overlapSize(0)==0),
    diffTolerance(dtol), resTolerance(rtol), maxIter(maxIt),
    defectGlobal(g, level, g.overlapSize(0)==0), num(0), minDt(mindt), goodIter(goodIt), maxIncreased(maxInc)
    {}

private:
    const G& grid;
    Model& model;
    FunctionType& u;
    FunctionType& f;
    OperatorAssembler& A;
    FunctionType uOldNewtonStep;
    double diffTolerance;
    double resTolerance;
    int maxIter;
    FunctionType defectGlobal;
    int num;
    double minDt;
    int goodIter;
    int maxIncreased;
};
}
#endif
