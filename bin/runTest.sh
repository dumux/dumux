#! /bin/bash
#
# Runs a test from the test directory and compare the resulting VTU files.
#
# Usage:
#
# runTest.sh COMPAREFLAG/SCRIPT REFERENCE_RESULT_FILE TEST_RESULT_FILE TEST_BINARY TEST_ARGS
#

function usage() {
    echo "Usage:"
    echo
    echo "runTest.sh COMPAREFLAG/SCRIPT REFERENCE_RESULT_FILE TEST_RESULT_FILE TEST_BINARY [TEST_ARGS]"
    echo "COMPAREFLAG: fuzzyvtu         - uses the fuzzycomparevtu.py script"
    echo "             exact            - uses diff for exact compare of two files"
    echo "             PATH_TO_SCRIPT   - tries to execute the custom script"
};

CMAKE_SOURCE_DIR=$(dirname "$0")
COMPARE_FLAG="$1"
REFERENCE_RESULT="$2"
TEST_RESULT="$3"
TEST_BINARY="$4"
TEST_ARGS="${@:5:100}"
rm -fv $TEST_RESULT

# make sure we have at least 3 parameters
if test "$#" -lt 4; then
    echo "Wrong number of parameters"
    echo
    usage
    exit 1
fi

# make sure the reference result exists
if ! test -r "$REFERENCE_RESULT"; then
    echo "File $REFERENCE_RESULT does not exist or is not readable"
    echo
    usage
    exit 1
fi

#run the test
echo "######################"
echo "# Running test"
echo "######################"
$TEST_BINARY $TEST_ARGS
RETURN_CODE_TEST_BINARY=$?
if test "$RETURN_CODE_TEST_BINARY" != "0"; then
    echo "Return code: $RETURN_CODE_TEST_BINARY"
    exit $RETURN_CODE_TEST_BINARY
fi

# compare the results
echo "######################"
echo "# Comparing results"
echo "######################"
if ! test -r "$TEST_RESULT"; then
    echo "File $TEST_RESULT does not exist or is not readable"
    exit 1
fi

# running the compare script
NOT_EQUAL=false
if [ "$COMPARE_FLAG" = "fuzzyvtu" -o \
     "$COMPARE_FLAG" = "$CMAKE_SOURCE_DIR/fuzzycomparevtu.py" ]; then
    if ! python $CMAKE_SOURCE_DIR/fuzzycomparevtu.py "$REFERENCE_RESULT" "$TEST_RESULT"; then
        NOT_EQUAL=true
    fi
elif [ "$COMPARE_FLAG" = "exact" ]; then
    if ! diff "$REFERENCE_RESULT" "$TEST_RESULT"; then
        NOT_EQUAL=true
    fi
elif [ -e "$COMPARE_FLAG" ]; then
    if ! $COMPARE_FLAG "$REFERENCE_RESULT" "$TEST_RESULT"; then
        NOT_EQUAL=true
    fi
else
    echo
    echo "ERROR: $0 was not able to run the compare script:"
    echo "       $COMPARE_FLAG"
    echo
    exit 2
fi

# printing error message in case of failure
if [ "$NOT_EQUAL" = "true" ]; then
    echo "The files \"$TEST_RESULT\" and \"$REFERENCE_RESULT\" are different."
    echo "Make sure the contents of \"$TEST_RESULT\" are still valid and "
    echo "make it the reference result if necessary."
    exit 1
fi

# SUCCESS!!!!!!
echo "Result and reference result are identical"
exit 0
