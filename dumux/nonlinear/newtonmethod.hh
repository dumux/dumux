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
//			int components = Model::m;
//			int size = (*u).size();
//			double oneByMagnitude[components];
//			double initialNorm[components];
//			for (int comp = 0; comp < components; comp++)
//				initialNorm[comp] = 0;
//			for (int i = 0; i < size; i++)
//				for (int comp = 0; comp < components; comp++) {
//					initialNorm[comp] += (*u)[i][comp]*(*u)[i][comp];
//			}
//			for (int comp = 0; comp < components; comp++) {
//				initialNorm[comp] = sqrt(initialNorm[comp]);
//				oneByMagnitude[comp] = 1.0/std::max(initialNorm[comp], 1e-6);
//			}
			
			double initialTotalNorm = std::max((*u).two_norm(), 1e-6);
			double error = 1e100;			
			double dt = localJacobian.getDt();
			bool divided = false;

		     while (dt > minDt && error > tolerance) {
				double lambdaOld;
				double error = 1e100;
				int iter = 0;
				model.globalDefect(defectGlobal);
				double globalResiduumOld = (*defectGlobal).two_norm();
				while (error > tolerance && iter < maxIter) {
					iter ++;
			    	double lambda = 1;
					//printvector(std::cout, *uOldNewtonStep, "uOldNewtonStep", "row", 200, 1, 3);
					//printvector(std::cout, *(model.uOldTimeStep), "uOldTimeStep", "row", 200, 1, 3);
					*f = 0;
					localJacobian.clearVisited();
					*update = *u;
					A.assemble(localJacobian, update, f);
					//printmatrix(std::cout, *A, "global stiffness matrix", "row", 11, 3);
					std::cout << "matrix norm: " << (*A).infinity_norm() << std::endl;
					//printvector(std::cout, *f, "right hand side", "row", 200, 1, 3);
					model.solve();
//					double updateNorm[components];
//					for (int comp = 0; comp < components; comp++)
//						updateNorm[comp] = 0;
//					for (int i = 0; i < size; i++)
//						for (int comp = 0; comp < components; comp++) {
//							updateNorm[comp] += (*u)[i][comp]*(*u)[i][comp];
//					}
//					error = 0;
//					for (int comp = 0; comp < components; comp++) {
//						updateNorm[comp] = sqrt(updateNorm[comp]);
//						error += oneByMagnitude[comp]*updateNorm[comp];
//					}
					error = (*update).two_norm()/initialTotalNorm;
					std::cout << "update norm = " << error << std::endl;
					
					//printvector(std::cout, *u, "update", "row", 200, 1, 3);
					*update *= -lambda;
					*u += *update;

					//printvector(std::cout, *u, "u", "row", 200, 1, 3);
					model.globalDefect(defectGlobal);
					double globalResiduum = (*defectGlobal).two_norm();
					std::cout << "old res = " << globalResiduumOld << ", new res = " << globalResiduum << std::endl;
					while (globalResiduum > globalResiduumOld){
					    lambda *= 0.5;
					    std::cout << "Linesearch: lambda divided to " << lambda << std::endl;
					    if (lambda < 1e-2)
					    	break;
					    *update *= 0.5;
					    *u -= *update;
					    globalResiduumOld = globalResiduum; 
						model.globalDefect(defectGlobal);
						globalResiduum = (*defectGlobal).two_norm();					    
						std::cout << "old res = " << globalResiduumOld << ", new res = " << globalResiduum << std::endl;
					}
					
					globalResiduumOld = globalResiduum;
					
					if (verbose)
						std::cout << "Newton step " << iter << ", defect = " << error << std::endl;
				}
				
				if (error > tolerance) {
					std::cout << "NewtonMethod::execute(), tolerance = " << tolerance 
						<< ": did not converge in " << iter << " iterations" << std::endl; 
					dt  = 0.5*dt;
					std::cout << "retry with reduced time step size of " << dt << std::endl;
					localJacobian.setDt(dt);
					*u = *(model.uOldTimeStep);
					divided = true;
				}
				else { 
					if (!divided && iter < goodIter) {
						dt = 2.0*dt;
						std::cout << "Below " << goodIter 
							<< " Newton iterations. Time step size doubled to " << dt << std::endl;
					}
					localJacobian.setDt(dt);
						
					return;
				}
			}
			
			if (dt <= minDt) 
				DUNE_THROW(MathError, "NewtonMethod:: time step size below minimum " << minDt << ".");
			
			
			return;
		}
		
		NewtonMethod(const G& g, Model& mod, double tol = 1e-5, int maxIt = 12, double mind = 1e-5, int goodIt = 5)
		: grid(g), model(mod), u(mod.u), f(mod.f), A(mod.A), localJacobian(mod.localJacobian), 
		  tolerance(tol), maxIter(maxIt), minDt(mind), 
		  goodIter(goodIt),defectGlobal(g), update(g)
		{ }
		
		NewtonMethod(const G& g, Model& mod, int level, double tol = 1e-5, int maxIt = 12, double mind = 1e-5, int goodIt = 5)
		: grid(g), model(mod), u(mod.u), f(mod.f), A(mod.A), localJacobian(mod.localJacobian), 
		  tolerance(tol), maxIter(maxIt), minDt(mind), 
		  goodIter(goodIt),defectGlobal(g,level), update(g, level)
		{ }
		
	private:
		const G& grid;
		Model& model;
		FunctionType& u;
		FunctionType& f;
		OperatorAssembler& A;
		LocalJacobian& localJacobian;
		double tolerance;
		int maxIter;
		double minDt;
		int goodIter;
	    FunctionType  defectGlobal;
	    FunctionType  update;
	};
}
#endif
