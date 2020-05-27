/* begin dumux
   put the definitions for config.h specific to
   your project here. Everything above will be
   overwritten
*/

/* begin private */
/* Name of package */
#define PACKAGE "@DUNE_MOD_NAME@"

/* Define to the address where bug reports for this package should be sent. */
#define PACKAGE_BUGREPORT "@DUNE_MAINTAINER@"

/* Define to the full name of this package. */
#define PACKAGE_NAME "@DUNE_MOD_NAME@"

/* Define to the full name and version of this package. */
#define PACKAGE_STRING "@DUNE_MOD_NAME@ @DUNE_MOD_VERSION@"

/* Define to the one symbol short name of this package. */
#define PACKAGE_TARNAME "@DUNE_MOD_NAME@"

/* Define to the home page for this package. */
#define PACKAGE_URL "@DUNE_MOD_URL@"

/* Define to the version of this package. */
#define PACKAGE_VERSION "@DUNE_MOD_VERSION@"

/* end private */

/* Define to the version of dumux */
#define DUMUX_VERSION "${DUMUX_VERSION}"

/* Define to the major version of dumux */
#define DUMUX_VERSION_MAJOR ${DUMUX_VERSION_MAJOR}

/* Define to the minor version of dumux */
#define DUMUX_VERSION_MINOR ${DUMUX_VERSION_MINOR}

/* Define to the revision of dumux */
#define DUMUX_VERSION_REVISION ${DUMUX_VERSION_REVISION}

/* Define the path to dumux */
#define DUMUX_SOURCE_DIR "${CMAKE_SOURCE_DIR}"

/* Define to 1 if gnuplot was found */
#cmakedefine HAVE_GNUPLOT 1

/* Define path to gnuplot executable */
#cmakedefine GNUPLOT_EXECUTABLE "@GNUPLOT_EXECUTABLE@"

/* Define to 1 if gstat was found */
#cmakedefine HAVE_GSTAT 1

/* Define path to gstat executable */
#cmakedefine GSTAT_EXECUTABLE "@GSTAT_EXECUTABLE@"

/* Define to 1 if gmsh was found */
#cmakedefine HAVE_GMSH 1

/* Define path to gmsh executable */
#cmakedefine GMSH_EXECUTABLE "@GMSH_EXECUTABLE@"

/* Defines whether pvpython has been found */
#cmakedefine HAVE_PVPYTHON 1

/* Define the path to pvpython */
#define PVPYTHON_EXECUTABLE "${PVPYTHON_EXECUTABLE}"

/* Define to 1 if quadmath was found */
#cmakedefine HAVE_QUAD 1

/* end dumux
   Everything below here will be overwritten
*/
