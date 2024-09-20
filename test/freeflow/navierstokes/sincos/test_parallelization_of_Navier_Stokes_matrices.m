## Test the parallelization of a Navier-Stokes test case with a
## staggered grid discretization for a grid with 2x4 cells.
##
## Sequential    Parallel subgrids
##              process P1   process P2
##   - -                       - -
##  | | |                     | | |
##   - -           - -         - -
##  | | |         | | |       | | |
##   - -           - -         - -
##  | | |         | | |       | | |
##   - -           - -         - -
##  | | |         | | |
##   - -           - -
##
## The local assembly of a sequential matrix M on a parallel process
## should be performed as
##     ---------------------
##    |                |    |
##    |                |    |
##    |                |    |
##    |       M_ii     |M_ij|  owned dofs
##    |                |    |
##    |                |    |
##    |                |    |
##     ---------------------
##    |        0       | I  |  not-owned dofs
##    |                |    |
##     ---------------------
## according to the article of Bastian and Blatt "On the Generic
## Parallelisation of Iterative Solvers for the Finite Element Method"
##
## How to use this test:
## - Start a sequential run with 1 process that assembles and stores the
##  matrices on the 2x4 grid as above:
##   ~$ mpirun -np 1 ./test_ff_navierstokes_sincos_uzawa ./params_for_parallel_matrix_test.input
## - Then start the same program in parallel with 2 processes:
##   ~$ mpirun -np 2 ./test_ff_navierstokes_sincos_uzawa ./params_for_parallel_matrix_test.input
## - Finally call this test function:
##   ~$ octave test_parallelization_of_Navier_Stokes_matrices.m

function [A,A0,A1,B,B0,B1,C,C0,C1,D,D0,D1]=test_parallelization_of_Navier_Stokes_matrices()
  disp('############################################################### ')
  disp('## Test the parallelization of a Navier-Stokes test problem  ## ')
  disp('## discretized by a staggered grid discretization for a grid ## ')
  disp('## with 2x4 cells.                                           ## ')
  disp('############################################################### ')
  disp('The test of the parallelization is performed by comparing the   ')
  disp('sequential block matrix with the block matrices of the two      ')
  disp('sub-processes P0 and P1, i.e.,                                  ')
  disp('                              process P1         process P2     ')
  disp('      / A B \        vs.       / A0 B0 \         / A1 B1 \      ')
  disp('      \ C D /                  \ C0 D0 /         \ C1 D1 /      ')
  disp(' Sequential matrix                Parallel sub-matrices         ')
  disp('                                                                ')
  disp('The matrices are assemebled on each of the following sequential ')
  disp('grid and parallel sub-grids:                                    ')
  disp('                              process P1         process P2     ')
  disp('    --9-- -10--                                 --7-- --8--     ')
  disp('   |     |     |                               |     |     |    ')
  disp('  20    21    22                              15    16    17    ')
  disp('   |     |     |                               |     |     |    ')
  disp('    --7-- --8--               --7-- --8--       --5-- --6--     ')
  disp('   |     |     |             |     |     |     |     |     |    ')
  disp('  17    18    19            15    16    17    12    13    14    ')
  disp('   |     |     |             |     |     |     |     |     |    ')
  disp('    --5-- --6--      vs.      --5-- --6--       --3-- --4--     ')
  disp('   |     |     |             |     |     |     |     |     |    ')
  disp('  14    15    16            12    13    14     9    10    11    ')
  disp('   |     |     |             |     |     |     |     |     |    ')
  disp('    --3-- --4--               --3-- --4--       --1-- --2--     ')
  disp('   |     |     |             |     |     |                      ')
  disp('  11    12    13             9    10    11                      ')
  disp('   |     |     |             |     |     |                      ')
  disp('    --1-- --2--               --1-- --2--                       ')
  disp('  Sequential grid                   Parallel sub-grids          ')
  disp('                                                                ')
  disp('The numbering is with respect to the DoFs in the staggered grid.')
  disp('                                                                ')
  disp('Keep in mind, that in Octave/Matlab the numbering of arrays     ')
  disp('starts with 1,2,... instead of 0,1,... (as e.g. in C++). The    ')
  disp('output of this test programm is accordingly.                    ')
  disp('                                                                ')
  disp('The local assembly of a sequential matrix M on a parallel       ')
  disp('process is according to the article of Bastian and Blatt "On the')
  disp('Generic Parallelisation of Iterative Solvers for the Finite     ')
  disp('Element Method" and should be performed as                      ')
  disp('    ---------------------                                       ')
  disp('   |                |    |                                      ')
  disp('   |                |    |                                      ')
  disp('   |                |    |                                      ')
  disp('   |       M_ii     |M_ij|  owned dofs                          ')
  disp('   |                |    |                                      ')
  disp('   |                |    |                                      ')
  disp('   |                |    |                                      ')
  disp('    ---------------------                                       ')
  disp('   |        0       | I  |  not-owned dofs                      ')
  disp('   |                |    |                                      ')
  disp('    ---------------------                                       ')
  disp('      owned dofs     not-owned                                  ')
  disp('                     dofs                                       ')
  disp('')

  ## Verbosity of the test output
  verbosity = 1;

  ## Formatting
  format('short')
  output_precision(2);

  ## Load all matrices of the sequential and parallel runs
  disp('Load matrices...')
  load matrixA_s1_r0_iter0.dat	# the sequential matrix
  load matrixA_s2_r0_iter0.dat	# the parallel matrix of process P0
  load matrixA_s2_r1_iter0.dat	# the parallel matrix of process P1
  A =_get_matrix(matrixA_s1_r0_iter0);
  A0=_get_matrix(matrixA_s2_r0_iter0);
  A1=_get_matrix(matrixA_s2_r1_iter0);
  load matrixB_s1_r0_iter0.dat
  load matrixB_s2_r0_iter0.dat
  load matrixB_s2_r1_iter0.dat
  B =_get_matrix(matrixB_s1_r0_iter0);
  B0=_get_matrix(matrixB_s2_r0_iter0);
  B1=_get_matrix(matrixB_s2_r1_iter0);
  load matrixC_s1_r0_iter0.dat
  load matrixC_s2_r0_iter0.dat
  load matrixC_s2_r1_iter0.dat
  C =_get_matrix(matrixC_s1_r0_iter0);
  C0=_get_matrix(matrixC_s2_r0_iter0);
  C1=_get_matrix(matrixC_s2_r1_iter0);
  load matrixD_s1_r0_iter0.dat
  load matrixD_s2_r0_iter0.dat
  load matrixD_s2_r1_iter0.dat
  D =_get_matrix(matrixD_s1_r0_iter0);
  D0=_get_matrix(matrixD_s2_r0_iter0);
  D1=_get_matrix(matrixD_s2_r1_iter0);

  ## Define mapping from the owned dofs of the parallel/local processes
  ## to the owned dofs of the global matrix and vice versa such that
  ## restrictions and prolongations can be performed
  disp('Define mappings...')

  ## Mapping: global indices of the parallel process in the sequential
  ## matrix (owned dofs)
  idx_u0=[1:6,11:16];
  idx_u1=[7:10,17:22];
  idx_p0=[1:4];
  idx_p1=[5:8];
  ## Mapping: global indices of the parallel process in the sequential
  ## matrix (not owned dofs)
  idx_u0_not_owned=[7:8,17:19];
  idx_u1_not_owned=[3:6,14:16];
  idx_p0_not_owned=[5:6];
  idx_p1_not_owned=[3:4];
  ## Mapping: global indices of the parallel process in the sequential
  ## matrix (owned dofs + not owned dofs)
  idx_u0_all=[1:8,11:19];
  idx_u1_all=[3:10,14:22];
  idx_p0_all=[1:6];
  idx_p1_all=[3:8];
  ## Mapping: local indices of the parallel process in the parallel
  ## matrix (owned dofs)
  local_idx_u0 = [1:6,9:14];
  local_idx_u1 = [5:8,12:17];
  local_idx_p0 = [1:4];
  local_idx_p1 = [3:6];
  ## Mapping: global indices of the parallel process in the sequential
  ## matrix (not owned dofs)
  local_idx_u0_not_owned=[7:8,15:17];
  local_idx_u1_not_owned=[1:4,9:11];
  local_idx_p0_not_owned=[5:6];
  local_idx_p1_not_owned=[1:2];

  ## Collect the mappings in structs
  ## global indices
  idx_u = {idx_u0, idx_u1};
  idx_p = {idx_p0, idx_p1};
  idx_u_not_owned = {idx_u0_not_owned, idx_u1_not_owned};
  idx_p_not_owned = {idx_p0_not_owned, idx_p1_not_owned};
  idx_u_all = {idx_u0_all, idx_u1_all};
  idx_p_all = {idx_p0_all, idx_p1_all};
  ## local indices
  local_idx_u = {local_idx_u0, local_idx_u1};
  local_idx_p = {local_idx_p0, local_idx_p1};
  local_idx_u_not_owned = {local_idx_u0_not_owned, local_idx_u1_not_owned};
  local_idx_p_not_owned = {local_idx_p0_not_owned, local_idx_p1_not_owned};

  disp('Start tests...')
  disp('')
  ## Compare matrix entries of the local/parallel matrices with the
  ## global/sequential matrix
  disp('############')
  disp('#  TEST A  #')
  disp('############')
  _test_matrix_entries(A, idx_u, idx_u, idx_u_not_owned, idx_u_not_owned,
		       {A0, A1},
		       local_idx_u, local_idx_u,
		       local_idx_u_not_owned, local_idx_u_not_owned,
		       idx_u_all, idx_u_all,
		       verbosity)

  disp('############')
  disp('#  TEST B  #')
  disp('############')
  _test_matrix_entries(B, idx_u, idx_p, idx_u_not_owned, idx_p_not_owned,
		       {B0, B1},
		       local_idx_u, local_idx_p,
		       local_idx_u_not_owned, local_idx_p_not_owned,
		       idx_u_all, idx_p_all,
		       verbosity)

  disp('############')
  disp('#  TEST C  #')
  disp('############')
  _test_matrix_entries(C, idx_p, idx_u, idx_p_not_owned, idx_u_not_owned,
		       {C0, C1},
		       local_idx_p, local_idx_u,
		       local_idx_p_not_owned, local_idx_u_not_owned,
		       idx_p_all, idx_u_all,
		       verbosity)

  disp('############')
  disp('#  TEST D  #')
  disp('############')
  _test_matrix_entries(D, idx_p, idx_p, idx_p_not_owned, idx_p_not_owned,
		       {D0, D1},
		       local_idx_p, local_idx_p,
		       local_idx_p_not_owned, local_idx_p_not_owned,
		       idx_p_all, idx_p_all,
		       verbosity)
endfunction

######################
## Helper functions ##
######################

function M=_get_matrix(mat)
  M=full(spconvert(mat));
endfunction

function idx=_get_indices_neq_zero(M)
  [m,n] = size(M);
  idx={};
  for i=1:m
    for j=1:n
      if (M(i,j) != 0)
	idx(end+1)={[i,j]};
      end
    end
  end
endfunction

function _display_indices_neq_zero(str_what, local_M, global_M,
				   diff_M,
				   local_idx_m,local_idx_n,
				   global_idx_m,global_idx_n)
  idx=_get_indices_neq_zero(diff_M);
  str_local  = '     ';
  str_global = '     ';
  [k,l]=size(idx);
  for (i=1:l)
    local_i = local_idx_m(idx{i}(1));
    local_j = local_idx_n(idx{i}(2));
    str_local  = cstrcat(str_local,' [',
			 num2str(local_i), ' ',
			 num2str(local_j), '] ',
			 num2str(local_M(local_i,local_j)), ', ');
    global_i = global_idx_m(idx{i}(1));
    global_j = global_idx_n(idx{i}(2));
    str_global = cstrcat(str_global,' [',
			 num2str(global_i), ' ',
			 num2str(global_j), '] ',
			 num2str(global_M(global_i,global_j)), ', ');
  endfor
  disp(cstrcat(str_what,': Different matrix entries!'))
  disp('    - Indices and value of local (parallel) matrix:   ')
  disp(str_local)
  disp('    - Indices and value of global (sequential) matrix:')
  disp(str_global)
endfunction

## - M: the sequential (global) matrix
## - idx_m, idx_n: the indices of the local dofs in the sequential
##   matrix
## - local_M: the local parallel matrix
## - local_idx_m, local_idx_n: the indices of the local dofs in the
##   parallel (local) matrix
function _compare_matrix_entries(str_what,
				 M, idx_m, idx_n,
				 local_M, local_idx_m, local_idx_n)
  ## Difference between sequential and parallel matrix restricted to
  ## the indices of the local process
  diff = M(idx_m,idx_n) - local_M(local_idx_m,local_idx_n);
  ## Largest error
  error = max(max(abs(diff)));
  ## Output of wrong indices, if there are any
  if (error)
    _display_indices_neq_zero(str_what, local_M, M,
			      diff,
			      local_idx_m, local_idx_n,
			      idx_m, idx_n)
  else
    disp(cstrcat(str_what,': matrix entries are correct'))
  endif
endfunction

function GM=_global_matrix(local_matrix,
			   global_idx_rows, global_idx_cols,
			   global_rows, global_cols)
  GM=zeros(global_rows,global_cols);
  [rows,cols] = size(local_matrix);
  for i=1:rows
    for j=1:cols
      GM(global_idx_rows(i),global_idx_cols(j)) = local_matrix(i,j);
    endfor
  endfor
endfunction

function spM=_summed_prolongated_matrix(M0, M1,
					global_idx_0_rows, global_idx_0_cols,
					global_idx_1_rows, global_idx_1_cols,
					global_rows, global_cols)
  prolongated_M0 = _global_matrix(M0, global_idx_0_rows, global_idx_0_cols,
				  global_rows, global_cols);
  prolongated_M1 = _global_matrix(M1, global_idx_1_rows, global_idx_1_cols,
				  global_rows, global_cols);
  spM = prolongated_M0 + prolongated_M1;
endfunction

function _compare_summed_prolongated_matrix(str_what,
					    M, M0, M1,
					    global_idx_0_rows,global_idx_0_cols,
					    global_idx_1_rows,global_idx_1_cols)
  [global_rows,global_cols] = size(M);
  spM = _summed_prolongated_matrix(M0, M1,
				   global_idx_0_rows, global_idx_0_cols,
				   global_idx_1_rows, global_idx_1_cols,
				   global_rows, global_cols);
  global_idx_rows = [1:global_rows];
  global_idx_cols = [1:global_cols];
  _compare_matrix_entries(str_what,
			  M, global_idx_rows, global_idx_cols,
			  spM, global_idx_rows, global_idx_cols);
endfunction

## - M is the sequential (global) matrix
## - local_M is an array with all parallel matrices
## - All indices are structs in this function, e.g.,
##   - idx_m = {idx_x0, idx_x1}
##   - local_idx_m = {local_idx_x0, local_idx_x1}
function _test_matrix_entries(M,
			      idx_m, idx_n,
			      idx_m_not_owned, idx_n_not_owned,
			      local_M,
			      local_idx_m, local_idx_n,
			      local_idx_m_not_owned, local_idx_n_not_owned,
			      local_idx_m_all, local_idx_n_all,
			      verbosity=1)
  for i = 1:length(local_M)
    disp('----------------------------------------------------------------')
    disp(cstrcat('Compare local matrix of process P',
		 num2str(i), ' with restricted global matrix'))
    disp('----------------------------------------------------------------')
    _compare_matrix_entries('Owned rows and owned cols',
			    M, idx_m{i}, idx_n{i},
			    local_M{i}, local_idx_m{i}, local_idx_n{i});
    _compare_matrix_entries('Owned rows and not-owned cols',
			    M, idx_m{i}, idx_n_not_owned{i},
			    local_M{i}, local_idx_m{i}, local_idx_n_not_owned{i});

    if (verbosity > 1)
      _compare_matrix_entries('Not-owned rows and owned cols (compare with global matrix)',
    			      M, ## zeros(size(M)),
    			      idx_m_not_owned{i}, idx_n{i},
    			      local_M{i}, local_idx_m_not_owned{i}, local_idx_n{i});
    endif
    if (verbosity > 2)
      _compare_matrix_entries('Not-owned rows and owned cols (compare with zero matrix)',
    			      zeros(size(M)),
    			      idx_m_not_owned{i}, idx_n{i},
    			      local_M{i}, local_idx_m_not_owned{i}, local_idx_n{i});
    endif

    if (verbosity > 1)
      _compare_matrix_entries('Not-owned dofs only (compare with global matrix)',
    			      M, ## eye(size(M)),
    			      idx_m_not_owned{i}, idx_n_not_owned{i},
    			      local_M{i},
    			      local_idx_m_not_owned{i}, local_idx_n_not_owned{i});
    endif
    if (verbosity > 2)
      _compare_matrix_entries('Not-owned dofs only (compare with identity matrix)',
    			      eye(size(M)),
    			      idx_m_not_owned{i}, idx_n_not_owned{i},
    			      local_M{i},
    			      local_idx_m_not_owned{i}, local_idx_n_not_owned{i});
    endif
  endfor

  if (verbosity > 2)
    disp('-----------------------------------------------------------------------------')
    disp(cstrcat('Compare global matrix with summed prolongated local ',
		 'matrices of all processes'))
    disp('-----------------------------------------------------------------------------')
    _compare_summed_prolongated_matrix('All dofs',
				       M, local_M{1}, local_M{2},
				       local_idx_m_all{1}, local_idx_n_all{1},
				       local_idx_m_all{2}, local_idx_n_all{2})
    endif
  disp('')  
endfunction

