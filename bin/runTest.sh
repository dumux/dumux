#! /bin/bash
#
# Runs a test from the test directory and compare the resulting VTU files.
#
# Usage:
#
# runTest.sh FUZZY_COMPARE REFERENCE_RESULT_FILE TEST_RESULT_FILE TEST_BINARY TEST_ARGS
#

function usage() {
    echo "Usage:"
    echo
    echo "runTest.sh FUZZY_COMPARE REFERENCE_RESULT_FILE TEST_RESULT_FILE TEST_BINARY [TEST_ARGS]"
};

FUZZY_COMPARE="$1"
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

if ! python "$FUZZY_COMPARE" "$REFERENCE_RESULT" "$TEST_RESULT"; then
    echo "The files \"$TEST_RESULT\" and \"$REFERENCE_RESULT\" are different."
    echo "Make sure the contents of \"$TEST_RESULT\" are still valid and "
    echo "make it the reference result if necessary."
    exit 1
fi

# SUCCESS!!!!!!
echo "Result and reference result are identical" 
exit 0
