#! /bin/bash

for TMP in $@; do
    echo "File: $TMP"
    sed -i \
"
s/\t/    /g;  # replace all tabs by 4 spaces
s/ *$//g; # remove all spaces at the end of a line
" \
$TMP

done
