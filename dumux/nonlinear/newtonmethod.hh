#ifndef DUNE_NEWTONMETHOD_HH
#define DUNE_NEWTONMETHOD_HH

namespace Dune {
	template<class G, class Model>
	class NewtonMethod 
	{
		typedef typename Model::FunctionType FunctionType;
		typedef typename Model::OperatorAssembler OperatorAssembler;
		typedef typename Model::LocalJacobian LocalJacobian;
	public:
		void execute(bool verbose = true)
		{
			double oneByMagnitude = 1.0/std::max((*u).two_norm(), 1e-5);
			double error = 1e100;
			int iter = 0;
			while (error > tolerance && iter < maxIter) {
				iter ++;
				*uOldNewtonStep = *u;
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
				*u *= -1;
				*u += *uOldNewtonStep;
				//printvector(std::cout, *u, "u", "row", 200, 1, 3);
				if (verbose)
					std::cout << "Newton step " << iter << ", defect = " << error << std::endl;
			}
			
			if (error > tolerance) {
				char message[120];
				sprintf(message, "NewtonMethod::execute(), tolerance = %g: did not converge in %d iterations", 
						tolerance, iter);
				DUNE_THROW(MathError, message);	
			}
			
			return;
		}
		
		NewtonMethod(const G& g, Model& mod, double tol = 1e-5, int maxIt = 30)
		: grid(g), model(mod), u(mod.u), f(mod.f), A(mod.A), localJacobian(mod.localJacobian), 
		  uOldNewtonStep(g), tolerance(tol), maxIter(maxIt)
		{ }
		
		NewtonMethod(const G& g, Model& mod, int level, double tol = 1e-5, int maxIt = 30)
		: grid(g), model(mod), u(mod.u), f(mod.f), A(mod.A), localJacobian(mod.localJacobian), 
		  uOldNewtonStep(g, level), tolerance(tol), maxIter(maxIt)
		{ }
		
	private:
		const G& grid;
		Model& model;
		FunctionType& u;
		FunctionType& f;
		OperatorAssembler& A;
		LocalJacobian& localJacobian;
		FunctionType uOldNewtonStep;
		double tolerance;
		int maxIter;
	};
}
#endif
