/* config.h.  Generated from config.h.in by configure.  */
/* config.h.in.  Generated from configure.ac by autoheader.  */

/* Dimension of ALBERTA grid */
#define ALBERTA_DIM 2

/* Define to the version of _dune_module. */
#define DUNE_COMMON_VERSION "1.2svn"

/* how to create a deprecated warning */
#define DUNE_DEPRECATED 

/* Define to the version of _dune_module. */
#define DUNE_DISC_VERSION "1.1svn"

/* Define to 1 if the experimental expression templates should be used */
/* #undef DUNE_EXPRESSIONTEMPLATES */

/* Define to the version of _dune_module. */
#define DUNE_GRID_VERSION "1.2svn"

/* Define to the version of _dune_module. */
#define DUNE_ISTL_VERSION "1.2svn"

/* Standard debug streams with a level below will collapse to doing nothing */
#define DUNE_MINIMAL_DEBUG_LEVEL 4

/* Define to dummy `main' function (if any) required to link to the Fortran
   libraries. */
/* #undef F77_DUMMY_MAIN */

/* Define if F77 and FC dummy `main' functions are identical. */
/* #undef FC_DUMMY_MAIN_EQ_F77 */

/* This is only true if alberta-library was found by configure _and_ if the
   application uses the ALBERTA_CPPFLAGS */
#define HAVE_ALBERTA ENABLE_ALBERTA

/* This is only true if alugrid-library was found by configure _and_ if the
   application uses the ALUGRID_CPPFLAGS */
#define HAVE_ALUGRID ENABLE_ALUGRID

/* Define to 1 if you have the <alugrid_parallel.h> header file. */
/* #undef HAVE_ALUGRID_PARALLEL_H */

/* Define to 1 if you have the <alugrid_serial.h> header file. */
#define HAVE_ALUGRID_SERIAL_H 1

/* Define to 1 if amiramesh-library is found */
/* #undef HAVE_AMIRAMESH */

/* Use the Apple OpenGL framework. */
/* #undef HAVE_APPLE_OPENGL_FRAMEWORK */

/* Define to 1 if you have the <array> header file. */
/* #undef HAVE_ARRAY */

/* Define if you have a BLAS library. */
#define HAVE_BLAS 1

/* Define to 1 if you have the <dlfcn.h> header file. */
#define HAVE_DLFCN_H 1

/* Define to 1 if _dune_module was found */
#define HAVE_DUNE_COMMON 1

/* Define to 1 if _dune_module was found */
#define HAVE_DUNE_DISC 1

/* Define to 1 if _dune_module was found */
#define HAVE_DUNE_GRID 1

/* Define to 1 if _dune_module was found */
#define HAVE_DUNE_ISTL 1

/* Define to 1 if grape-library is found */
/* #undef HAVE_GRAPE */

/* Define to 1 if hdf5 was found */
/* #undef HAVE_HDF5 */

/* Define to 1 if you have the <inttypes.h> header file. */
#define HAVE_INTTYPES_H 1

/* Define if you have LAPACK library. */
#define HAVE_LAPACK 1

/* Define to 1 if you have the `m' library (-lm). */
#define HAVE_LIBM 1

/* Define to 1 if your system has a GNU libc compatible `malloc' function, and
   to 0 otherwise. */
#define HAVE_MALLOC 1

/* Define to 1 if you have the <malloc.h> header file. */
#define HAVE_MALLOC_H 1

/* Define to 1 if you have the <memory.h> header file. */
#define HAVE_MEMORY_H 1

/* Define if you have METIS library */
/* #undef HAVE_METIS */

/* Define if you have the MPI library. This is only true if MPI was found by
   configure _and_ if the application uses the MPI_CPPFLAGS */
/* #undef HAVE_MPI */

/* Define if you have a PARDISO library. */
#define HAVE_PARDISO 1

/* Define to 1 if PARMETIS is found */
/* #undef HAVE_PARMETIS */

/* Define to 1 if psurface-library is found */
/* #undef HAVE_PSURFACE */

/* Define if you have POSIX threads libraries and header files. */
#define HAVE_PTHREAD 1

/* Define to 1 if you have the <rpc/rpc.h> header file. */
#define HAVE_RPC_RPC_H 1

/* Define to 1 if you have the `sqrt' function. */
#define HAVE_SQRT 1

/* Define to 1 if you have the <stdint.h> header file. */
#define HAVE_STDINT_H 1

/* Define to 1 if you have the <stdlib.h> header file. */
#define HAVE_STDLIB_H 1

/* Define to 1 if you have the `strchr' function. */
#define HAVE_STRCHR 1

/* Define to 1 if you have the <strings.h> header file. */
#define HAVE_STRINGS_H 1

/* Define to 1 if you have the <string.h> header file. */
#define HAVE_STRING_H 1

/* Define to 1 if SUPERLU is found */
/* #undef HAVE_SUPERLU */

/* Define to 1 if SUPERLU_DIST is found */
/* #undef HAVE_SUPERLU_DIST */

/* Define to 1 if you have the <sys/stat.h> header file. */
#define HAVE_SYS_STAT_H 1

/* Define to 1 if you have the <sys/types.h> header file. */
#define HAVE_SYS_TYPES_H 1

/* Define to 1 if you have the <tr1/array> header file. */
/* #undef HAVE_TR1_ARRAY */

/* Define to 1 if you have the <tr1/tuple> header file. */
/* #undef HAVE_TR1_TUPLE */

/* Define to 1 if you have the <tr1/type_traits> header file. */
/* #undef HAVE_TR1_TYPE_TRAITS */

/* Define to 1 if you have the <tuple> header file. */
/* #undef HAVE_TUPLE */

/* Define to 1 if you have the <type_traits> header file. */
/* #undef HAVE_TYPE_TRAITS */

/* This is only true if UG was found by configure _and_ if the application
   uses the UG_CPPFLAGS */
#define HAVE_UG ENABLE_UG

/* Define to 1 if you have the <unistd.h> header file. */
#define HAVE_UNISTD_H 1

/* Define to 1 if you have the <windows.h> header file. */
/* #undef HAVE_WINDOWS_H */

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT ""

/* Define to the full name of this package. */
#define PACKAGE_NAME ""

/* Define to the full name and version of this package. */
#define PACKAGE_STRING ""

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME ""

/* Define to the version of this package. */
#define PACKAGE_VERSION ""

/* Define to necessary symbol if this constant uses a non-standard name on
   your system. */
/* #undef PTHREAD_CREATE_JOINABLE */

/* Define to 1 if you have the ANSI C header files. */
#define STDC_HEADERS 1

/* define to 1 if there is a header slu_ddefs.h in SuperLU */
/* #undef SUPERLU_POST_2005_VERSION */

/* Define to 1 if your <sys/time.h> declares `struct tm'. */
/* #undef TM_IN_SYS_TIME */

/* use UG LGM domain */
/* #undef UG_LGMDOMAIN */

/* Define to 1 if the X Window System is missing or not being used. */
/* #undef X_DISPLAY_MISSING */

/* Define to empty if `const' does not conform to ANSI C. */
/* #undef const */

/* Define to `__inline__' or `__inline' if that's what the C compiler
   calls it, or to nothing if 'inline' is not supported under any name.  */
#ifndef __cplusplus
/* #undef inline */
#endif

/* Define to rpl_malloc if the replacement function should be used. */
/* #undef malloc */

/* Define to `unsigned int' if <sys/types.h> does not define. */
/* #undef size_t */
