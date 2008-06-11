#/bin/bash

# This file is a simplistic wrapper script which compiles the actual
# simulation binary if necessary, runs it and filters its output to
# something readable.
SIMNAME=($1)
SIMARGS=""
ARGS=($@)

N=1
while test "$N" -lt ${#ARGS[*]}; do
    SIMARGS="$SIMARGS ${ARGS[$N]}"
    N=$((N+1))
done

DUMUXDIR=~/temp/DUMUX
SOURCEDIR=$HOME/Stupid
BUILDDIR=$HOME/sb

CONFDIR="$SOURCEDIR/scripts"
CONFFILE="$CONFDIR/$SIMNAME.conf"

test -f $CONFFILE && source $CONFFILE

if test -z "$SIMRESULT"; then SIMRESULT="$SIMNAME.pvd"; fi
if test -z "$SIMBINARY"; then SIMBINARY="$BUILDDIR/Problems/$SIMNAME/$SIMNAME"; fi
if test -z "$FILTER"; then FILTER="timestep|took"; fi
if test -z "$TARGET"; then TARGET="release"; fi

cd $BUILDDIR
CMAKE="cmake -DCMAKE_BUILD_TYPE=$TARGET \
       -DDUNE_DIR=$DUMUXDIR \
       -DUG_DIR=$DUMUXDIR/external/ug \
       $SOURCEDIR"
echo $CMAKE
$CMAKE || exit 1

MAKE="make $SIMNAME"
$MAKE || exit 2

# remove simulation results from previous runs
ls *.vtp *.vtu *.pvd 2> /dev/null | xargs -d '\n' rm -f 

# run the simulation
set -o pipefail
$SIMBINARY $SIMARGS | tee "$SIMNAME.log" | grep -i -E "$FILTER" 
FAILED=$?

# start paraview
if test $FAILED -eq "0" && test -f "$SIMRESULT"; then
    paraview --data="$SIMRESULT"
fi
