#ifndef DUNE_IMPES_HH
#define DUNE_IMPES_HH

#include <dune/grid/io/file/vtk/vtkwriter.hh>
#include "dumux/fractionalflow/fractionalflow.hh"

/**
 * @file
 * @brief  IMPES scheme
 * @author Bernd Flemisch, last changed by Markus Wolff
 */

namespace Dune {
/**
 * \ingroup fracflow
 * @brief IMplicit Pressure Explicit Saturation (IMPES) scheme for the solution of
 * coupled diffusion/transport problems
 */

template<class G, class Diffusion, class Transport, class VC> class IMPES :
	public FractionalFlow<G, Diffusion, Transport, VC> {
	typedef typename Diffusion::RepresentationType PressType;
	typedef typename Diffusion::NumberType RT;

public:
	typedef typename Transport::RepresentationType RepresentationType;

	virtual void totalVelocity(const RT t=0) {
		this->calcTotalVelocity(t);
		return;
	}

	VC variables() const{
		return this->transproblem.variables;
	}

	virtual void initial() {
		double t = 0;
		this->initialTransport();
		if (this->calcpressure) {
			this->pressure(t);
			totalVelocity(t);
		}
//		return;
	}

	virtual int update(const RT t, RT& dt, RepresentationType& updateVec,
			RT cFLFactor = 1) {
		int pressSize = variables().pressure.size();
		PressType pressOldIter(variables().pressure);
		PressType pressHelp(pressSize);
		int satSize = variables().saturation.size();
		RepresentationType saturation(variables().saturation);
		RepresentationType satOldIter(variables().saturation);
		RepresentationType satHelp(satSize);
		RepresentationType satDiff(satSize);
		RepresentationType updateOldIter(satSize);
		RepresentationType updateHelp(satSize);
		RepresentationType updateDiff(satSize);

		bool converg = false;
		int iter = 0;
		int iterTot = 0;
		updateOldIter = 0;
		while (!converg) {
			iter++;
			iterTot++;
			if (!this->diffproblem.materialLaw.isLinear()
					|| this->diffproblem.capillary|| this->calcpressure) { // update pressure 
				pressure(t);
				totalVelocity(t);
			}
			Transport::update(t, dt, updateVec);
			if (iterFlag) { // only needed if iteration has to be done
				if (this->calcpressure) {
					variables().pressure *= omega;
					pressHelp = pressOldIter;
					pressHelp *= (1-omega);
					variables().pressure += pressHelp;
					pressOldIter = variables().pressure;
				}
				updateHelp = updateVec;
				saturation = variables().saturation;
				saturation += (updateHelp *= dt*cFLFactor);
				saturation *= omega;
				satHelp = satOldIter;
				satHelp *= (1-omega);
				saturation += satHelp;
				updateDiff = updateVec;
				updateDiff -= updateOldIter;
				satOldIter = saturation;
				updateOldIter = updateVec;
			}
			// break criteria for iteration loop
			if (iterFlag==2&& dt*updateDiff.two_norm()/(saturation).two_norm() <= maxDefect )
				converg = true;
			else if (iterFlag==2&& iter > nIter ) {
				std::cout << "Nonlinear loop in IMPES.update exceeded nIter = "
						<< nIter << " iterations."<< std::endl;
				return 1;
			} else if (iterFlag==1&& iter > nIter )
				converg = true;
			else if (iterFlag==0)
				converg = true;
		}
		// outputs
		if (iterFlag==2)
			std::cout << "Iteration steps: "<< iterTot << std::endl;
		std::cout.setf(std::ios::scientific, std::ios::floatfield);

		return 0;
	}

	virtual void vtkout(const char* name, int k) const {
		variables().vtkout(name, k);
		return;
	}

	//! Construct an IMPES object.
	IMPES(Diffusion& diff, Transport& trans, int flag = 0, int nIt = 2,
			double maxDef = 1e-5, double om = 1) :
		FractionalFlow<G, Diffusion, Transport, VC>(diff, trans),
				iterFlag(flag), nIter(nIt), maxDefect(maxDef), omega(om) {
	}

protected:
	const int iterFlag;
	const int nIter;
	const double maxDefect;
	const double omega;
};
}
#endif
