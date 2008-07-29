#!/bin/sh
# : duneproject 4887 2007-04-02 10:06:45Z christi $

# barf on errors
set -e

usage () {
    echo "Usage: ./autogen.sh DUNE_COMMON_DIR [options]"
    echo "  --ac=, --acversion=VERSION   use a specific VERSION of autoconf"
    echo "  --am=, --amversion=VERSION   use a specific VERSION of automake"
    echo "  -h,    --help                you already found this :-)"
}

for OPT in "$@"; do
    set +e
    # stolen from configure...
    # when no option is set, this returns an error code
    arg=`expr "x$OPT" : 'x[^=]*=\(.*\)'`
    set -e

    case "$OPT" in
	--ac=*|--acversion=*)
			if test "x$arg" = "x"; then
				usage; 
				exit 1;
			fi
			ACVERSION=$arg
			;;
	--am=*|--amversion=*)
			if test "x$arg" = "x"; then
				usage; 
				exit 1;
			fi
			AMVERSION=$arg
			;;
	-h|--help) usage ; exit 0 ;;
	*)
            if test -d "$OPT/m4"; then
              ACLOCAL_FLAGS="$ACLOCAL_FLAGS -I $OPT/m4"
            fi
            if test -d "$OPT/am"; then
              am_dir="$OPT/am"
            fi
            if test -d "$OPT/share/aclocal"; then
              ACLOCAL_FLAGS="$ACLOCAL_FLAGS -I $OPT/share/aclocal"
            fi
            if test -d "$OPT/share/dune-common/am"; then
              am_dir="$OPT/share/dune-common/am"
            fi
            PATH="$PATH:$OPT/bin"
            ;;
    esac
done

## report parameters
if test "x$ACVERSION" != "x"; then
	echo "Forcing autoconf version �$ACVERSION�"
	if ! which autoconf$ACVERSION > /dev/null; then
		echo
		echo "Error: Could not find autoconf$ACVERSION"
		echo "       Did you specify a wrong version?"
		exit 1
	fi
fi
if test "x$AMVERSION" != "x"; then
	echo "Forcing automake version �$AMVERSION�"
	if ! which automake$AMVERSION > /dev/null; then
		echo
		echo "Error: Could not find automake$AMVERSION"
		echo "       Did you specify a wrong version?"
		exit 1
	fi
fi


## run autotools

echo "--> libtoolize..."
# this script won't rewrite the files if they already exist. This is a
# PITA when you want to upgrade libtool, thus I'm setting --force
libtoolize --force

# prepare everything
echo "--> aclocal..."
aclocal$AMVERSION $ACLOCAL_FLAGS

# applications should provide a config.h for now
echo "--> autoheader..."
autoheader$ACVERSION

# create a link to the dune-common am directory
if [ -n "$am_dir" ] && [ -d $am_dir ]; then
  echo "--> linking dune-common/am..."
  rm -f am
  ln -s $am_dir am
else
  echo
  echo "Error: Could not find dune-common/am!"
  usage
  exit 1
fi

# call automake/conf
echo "--> automake..."
automake$AMVERSION --add-missing

echo "--> autoconf..."
autoconf$ACVERSION

## tell the user what to do next
echo "Now run ./configure "

