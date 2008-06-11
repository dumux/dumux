#ifndef DUNE_NEWTONMETHOD_HH
#define DUNE_NEWTONMETHOD_HH

namespace Dune {
template<class G, class Model> class NewtonMethod {
	typedef typename Model::FunctionType FunctionType;
	typedef typename Model::OperatorAssembler OperatorAssembler;
	typedef typename Model::LocalJacobian LocalJacobian;
public:
	void execute(bool verbose = true) {
		double uNorm = std::max((*u).two_norm(), 1.0);
		grid.comm().sum(&uNorm, 1);		
		double oneByMagnitude = 1.0/uNorm;
		double error = 1e100;
		double globalResiduum = 1e100;
		bool divided = false;
		double dt = localJacobian.getDt();


		while (dt > minDt && (error > difftolerance || globalResiduum > restolerance)) {
			double error = 1e100;
			int iter = 0;
			int iiter =0;
			double lambda = 1;
			double lambdaOld=lambda;			
			model.globalDefect(defectGlobal);
			double globalResiduum = 0.5*(*defectGlobal).two_norm();
			double globalResiduumOld=globalResiduum;
		    grid.comm().sum(&globalResiduum, 1);	    
		    if (grid.comm().rank() == 0)
				std::cout << "initial residual = " << globalResiduum << std::endl;		    	
			while ((error > difftolerance || globalResiduum > restolerance) && iter < maxIter) {
				iter ++;
				iiter=0;
				num=0;
				*uOldNewtonStep = *u;
				globalResiduumOld=globalResiduum;
				lambda=lambdaOld;
				//printvector(std::cout, *uOldNewtonStep, "uOldNewtonStep", "row", 200, 1, 3);
				//printvector(std::cout, *(model.uOldTimeStep), "uOldTimeStep", "row", 200, 1, 3);
				*f = 0;
				localJacobian.clearVisited();
				A.assemble(localJacobian, u, f);
				//printmatrix(std::cout, *A, "global stiffness matrix", "row", 11, 3);
				//printvector(std::cout, *f, "right hand side", "row", 200, 1, 3);

				model.solve();
				error = oneByMagnitude*((*u).two_norm());
				//printvector(std::cout, *u, "update", "row", 200, 1, 3);
				*u *= -lambda;
				*u += *uOldNewtonStep;
				//printvector(std::cout, *u, "u", "row", 200, 1, 3);
				model.globalDefect(defectGlobal);
				globalResiduum=0.5*(*defectGlobal).two_norm();
				grid.comm().sum(&globalResiduum, 1);
				//printvector(std::cout, *defectGlobal, "global Defect", "row", 200, 1, 3);

				while (globalResiduum >= globalResiduumOld && iiter < (maxIter
						/2)) {
					iiter++;
					*uOldNewtonStep = *u;
					globalResiduumOld=globalResiduum;
					if (num==0)
						lambdaOld=lambda;
					if (num<4) {
						lambda*= 0.5;
						num++;
					} else
						lambda*=0.99;
					*f = 0;
					localJacobian.clearVisited();
					A.assemble(localJacobian, u, f);
					model.solve();
					error = oneByMagnitude*((*u).two_norm());
					*u *= -lambda;
					*u += *uOldNewtonStep;
					model.globalDefect(defectGlobal);
					globalResiduum=0.5*(*defectGlobal).two_norm();
					grid.comm().sum(&globalResiduum, 1);
				}
				if (verbose && grid.comm().rank() == 0)
					std::cout << "Newton step "<< iter << ", residual = "
							<< globalResiduum << ", difference = "<< error
							<< std::endl;
			}
			if (error > difftolerance || globalResiduum > restolerance) {
				if (grid.comm().rank() == 0)
					std::cout << "NewtonMethod::execute(), tolerances = "
							<< difftolerance << " , " << restolerance << ": did not converge in "<< iter
							<< " iterations"<< std::endl;
				dt *= 0.5;
				if (grid.comm().rank() == 0)
					std::cout << "Retry same time step with reduced size of "
							<< dt << std::endl;
				localJacobian.setDt(dt);
				*u = *(model.uOldTimeStep);
				divided = true;
			} else {
				if (grid.comm().rank() == 0)
					std::cout << "Converged. Residual = "<< globalResiduum
							<< ", difference = "<< error << std::endl;
				if ((!divided) && iter < goodIter) {
					dt *= 2;
					if (grid.comm().rank() == 0)
						std::cout << "Below "<< goodIter
								<< " Newton iterations. Initial size for the next time step doubled to "
								<< dt << std::endl;
				}
				localJacobian.setDt(dt);

				return;
			}
		}

		if (dt <= minDt)
			DUNE_THROW(MathError,
					"NewtonMethod:: time step size below minimum " << minDt
							<< ".");
		return;
	}

	NewtonMethod(const G& g, Model& mod, double dtol = 1e-8,
			double rtol = 1e-2, int maxIt = 20, double mindt = 1e-5,
			int goodIt = 3) :
		grid(g), model(mod), u(mod.u), f(mod.f), A(mod.A),
				localJacobian(mod.localJacobian), uOldNewtonStep(g),
				difftolerance(dtol), restolerance(rtol), maxIter(maxIt),
				minDt(mindt), goodIter(goodIt), defectGlobal(g), num(0) {
	}

	NewtonMethod(const G& g, Model& mod, int level, double dtol = 1e-8,
			double rtol = 1e-5, int maxIt = 12, double mindt = 1e-5,
			int goodIt = 3) :
		grid(g), model(mod), u(mod.u), f(mod.f), A(mod.A),
				localJacobian(mod.localJacobian), uOldNewtonStep(g, level),
				difftolerance(dtol), restolerance(rtol), maxIter(maxIt),
				minDt(mindt), goodIter(goodIt), defectGlobal(g, level), num(0)

	{
	}

private:
	const G& grid;
	Model& model;
	FunctionType& u;
	FunctionType& f;
	OperatorAssembler& A;
	LocalJacobian& localJacobian;
	FunctionType uOldNewtonStep;
	double difftolerance;
	double restolerance;
	int maxIter;
	FunctionType defectGlobal;
	int num;
	double minDt;
	int goodIter;
};
}
#endif
