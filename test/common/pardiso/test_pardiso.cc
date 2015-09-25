// $Id$
/*****************************************************************************
 *   Copyright (C) 2008 by Bernd Flemisch                                    *
 *   Institute of Hydraulic Engineering                                      *
 *   University of Stuttgart, Germany                                        *
 *   email: andreas.lauser _at_ iws.uni-stuttgart.de                         *
 *                                                                           *
 *   This program is free software: you can redistribute it and/or modify    *
 *   it under the terms of the GNU General Public License as published by    *
 *   the Free Software Foundation, either version 2 of the License, or       *
 *   (at your option) any later version.                                     *
 *                                                                           *
 *   This program is distributed in the hope that it will be useful,         *
 *   but WITHOUT ANY WARRANTY; without even the implied warranty of          *
 *   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the           *
 *   GNU General Public License for more details.                            *
 *                                                                           *
 *   You should have received a copy of the GNU General Public License       *
 *   along with this program.  If not, see <http://www.gnu.org/licenses/>.   *
 *****************************************************************************/
/*!
 * \file
 *
 * \brief test for the Pardiso preconditioner
 */
#include "config.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "config.h"
#include <dune/istl/bcrsmatrix.hh>
#include <dune/istl/io.hh>
#include <dune/istl/operators.hh>
#include <dune/istl/solvers.hh>
#include "dumux/common/pardiso.hh"

int main(int argc, char** argv)
{
    try
    {
#ifdef HAVE_PARDISO
        /* Matrix data. */
        int n = 8;
        int ia[ 9] = { 0, 4, 7, 9, 11, 12, 15, 17, 20 };
        int ja[20] = { 0, 2, 5, 6,
                          1, 2, 4,
                          2, 7,
                          3, 6,
                          4,
                          2, 5, 7,
                          1, 6,
                          2, 6, 7 };
        double a[20] = { 7.0, 1.0, 2.0, 7.0,
                          -4.0, 8.0, 2.0,
                          1.0, 5.0,
                          7.0, 9.0,
                          -4.0,
                          7.0, 3.0, 8.0,
                          1.0, 11.0,
                          -3.0, 2.0, 5.0 };

        int nnz = ia[n];
        int      mtype = 11;        /* Real unsymmetric matrix */

        /* RHS and solution vectors. */
        double b[8], x[8];
        int      nrhs = 1;          /* Number of right hand sides. */

        /* Internal solver memory pointer pt,                  */
        /* 32-bit: int pt[64]; 64-bit: long int pt[64]         */
        /* or void *pt[64] should be OK on both architectures  */
        void    *pt[64];

        /* Pardiso control parameters. */
        int iparm[64];
        double dparm[64];
        int maxfct, mnum, phase, error, msglvl;

        /* Number of processors. */
        //int num_procs;

        /* Auxiliary variables. */
        //char    *var;
        int i;

        double   ddum;              /* Double dummy */
        int      idum;              /* Integer dummy. */
	int solver = 0;

        /* -------------------------------------------------------------------- */
        /* ..  Setup Pardiso control parameters.                                */
        /* -------------------------------------------------------------------- */

        F77_FUN(pardisoinit) (pt, &mtype, &solver, iparm, dparm, &error);
        if (error != 0) {
            printf("\nERROR during initialization: %d", error);
            exit(3);
        }

        iparm[2]  = 1;

        maxfct = 1;        /* Maximum number of numerical factorizations.  */
        mnum   = 1;         /* Which factorization to use. */

        msglvl = 0;         /* Print statistical information  */
        error  = 0;         /* Initialize error flag */

        /* -------------------------------------------------------------------- */
        /* ..  Convert matrix from 0-based C-notation to Fortran 1-based        */
        /*     notation.                                                        */
        /* -------------------------------------------------------------------- */
        for (i = 0; i < n+1; i++) {
            ia[i] += 1;
        }
        for (i = 0; i < nnz; i++) {
            ja[i] += 1;
        }

        phase = 13;

        iparm[7] = 1;       /* Max numbers of iterative refinement steps. */

        /* Set right hand side to one. */
        for (i = 0; i < n; i++) {
            b[i] = 1;
        }

        F77_FUN(pardiso) (pt, &maxfct, &mnum, &mtype, &phase,
                           &n, a, ia, ja, &idum, &nrhs,
                           iparm, &msglvl, b, x, &error, dparm);

        if (error != 0) {
            printf("\nERROR during solution: %d", error);
            exit(3);
        }

        printf("\nSolve completed ... ");
        printf("\nThe solution of the system is: ");
        for (i = 0; i < n; i++) {
            printf("\n x [%d] = % f", i, x[i] );
        }
        printf ("\n\n");

        /* -------------------------------------------------------------------- */
        /* ..  Termination and release of memory.                               */
        /* -------------------------------------------------------------------- */
        phase = -1;                 /* Release internal memory. */

        F77_FUN(pardiso) (pt, &maxfct, &mnum, &mtype, &phase,
                           &n, &ddum, ia, ja, &idum, &nrhs,
                           iparm, &msglvl, &ddum, &ddum, &error, dparm);

        typedef Dune::FieldMatrix<double,1,1> M;
        Dune::BCRSMatrix<M> B(8,8,Dune::BCRSMatrix<M>::random);

        // initially set row size for each row
        B.setrowsize(0,4);
        B.setrowsize(1,3);
        B.setrowsize(2,2);
        B.setrowsize(3,2);
        B.setrowsize(4,1);
        B.setrowsize(5,3);
        B.setrowsize(6,2);
        B.setrowsize(7,3);

        // finalize row setup phase
        B.endrowsizes();

        // add column entries to rows
        B.addindex(0,0); B.addindex(0,2); B.addindex(0,5); B.addindex(0,6);
        B.addindex(1,1); B.addindex(1,2); B.addindex(1,4);
        B.addindex(2,2); B.addindex(2,7);
        B.addindex(3,3); B.addindex(3,6);
        B.addindex(4,4);
        B.addindex(5,2); B.addindex(5,5); B.addindex(5,7);
        B.addindex(6,1); B.addindex(6,6);
        B.addindex(7,2); B.addindex(7,6); B.addindex(7,7);

        // finalize column setup phase
        B.endindices();

        // set entries using the random access operator
        B[0][0] = 7; B[0][2] = 1; B[0][5] = 2; B[0][6] = 7;
        B[1][1] = -4; B[1][2] = 8; B[1][4] = 2;
        B[2][2] = 1; B[2][7] = 5;
        B[3][3] = 7; B[3][6] = 9;
        B[4][4] = -4;
        B[5][2] = 7; B[5][5] = 3; B[5][7] = 8;
        B[6][1] = 1; B[6][6] = 11;
        B[7][2] = -3; B[7][6] = 2; B[7][7] = 5;

        //printmatrix(std::cout, B, "matrix B", "row", 9, 1);

        typedef Dune::FieldVector<double, 1> VB;
        typedef Dune::BlockVector<VB> Vector;
        typedef Dune::BCRSMatrix<M> Matrix;
        Dune::MatrixAdapter<Matrix,Vector,Vector> op(B);        // make linear operator from A
        Dune::SeqPardiso<Matrix,Vector,Vector> pardiso(B);        // preconditioner object
        Dune::LoopSolver<Vector> loop(op, pardiso, 1E-14, 2, 1);         // an inverse operator

        Vector f(n);
        f = 1;
        Vector y(n);
        y = 0;
        Dune::InverseOperatorResult r;
        loop.apply(y, f, r);

        std::cout << "\nSolve completed ... ";
        std::cout << "\nThe solution of the system is: ";
        for (i = 0; i < n; i++) {
            std::cout << "\n x [" << i << "] = " << y[i];
        }
        std::cout << "\n";

#else
        DUNE_THROW(Dune::NotImplemented, "no Pardiso library available, reconfigure with correct --with-pardiso options");
#endif


        return 0;
    }
    catch (Dune::Exception &e){
        std::cerr << "Dune reported error: " << e << std::endl;
    }
    catch (...){
        std::cerr << "Unknown exception thrown!" << std::endl;
    }
}
