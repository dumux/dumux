#!/bin/bash

help()
{
   echo "Prints all commit authors (contributors) between two commits."
   echo "Usage: $0 -from <commit sha or tag> -to <commit sha or tag>"
   echo "Example: $0 -from 2.12.0 -to d9efb9d80b1f794d976f6b5ae1ffb9f71a7dcdf9"
   echo -e "\t-from From where: The commit sha or other tags working with git log"
   echo -e "\t-to To where: The commit sha or other tags working with git log"
   exit 1 # Exit script after printing help
}

if [[ $# -eq 0 ]]
then help;
fi

while [[ $# -gt 0 ]]
do
key="$1"

case "$key" in
    -from|--from)
    from_sha="$2"
    shift
    shift
    ;;
    -to|--to)
    to_sha="$2"
    shift
    shift
    ;;
    ?)
    help # print help
    ;;
esac
done

echo "Contributors from commit ${from_sha} to ${to_sha}:"

git log $from_sha..$to_sha "$@" | grep ^Author: | sed 's/ <.*//; s/^Author: //' | sort | uniq -c | sort -nr
