#!/bin/bash
#
# Fix the indention of all source files given on the command
# line. This script requires GNU emacs to be installed.

function usage() {
    echo "$0 FILE_LIST"
    echo
    echo "Fix the indention of all source files given on the command"
    echo "line. This script requires GNU emacs to be installed."
    exit 1
}

if test "$#" -lt 1; then
    usage
fi

for TMP in $@; do
    echo "Processing file: $TMP..."
    echo "Sanitizing white space"
    sed -i \
"
s/^\(.*\)/ \1/; # add a space at the beginning of each line. this will
                # later be correctly indented by emacs
s/\t/    /g;  # replace all tabs by 4 spaces
s/ *$//g; # remove all spaces at the end of a line
" \
$TMP

    echo "Re-indenting everything"
    emacs --batch $TMP --eval "
(progn

(setq tab-width 4) ; make tabs 4 characters wide
(setq-default indent-tabs-mode nil) ; don't use tabs for indentation
(c-set-offset 'substatement-open 0) ; make braces for ifs and fors 
                                    ; start at the same column as the statement
(c-set-offset 'inline-open 0) ; don't indent inline opens
(c-set-offset 'innamespace 0) ; don't indent namespaces
(c-set-offset 'case-label 0) ; don't indent case labels further
                             ; than the rest of the function
                              ; body
(setq c-basic-offset 4)

;; Always ending a file with a newline may be useful in some scripts
(setq require-final-newline t)

(indent-region (point-min) (point-max) nil)
(save-buffer)
)"
done
