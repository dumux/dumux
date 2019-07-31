#!/bin/sh

if [[ ${1: -3} == ".hh" ]]; then
  sed '1,/define/d' $1 >tmpfile
else
  sed '1,/\*\*\*\*\*/d' $1 >tmpfile
fi

insidecodeblock=false

while IFS= read -r line
do
  strippedline=$(echo $line | sed "s/^[ \t]*//")
  if [[ ${strippedline:0:2} == "//" ]]; then
    if [[ $insidecodeblock == true ]]; then
      echo "\`\`\`"
      insidecodeblock=false
    fi
    echo "${strippedline:2}"
  else
    if [[ $insidecodeblock == false ]]; then
      echo "\`\`\`cpp"
      insidecodeblock=true
    fi
    echo "$line"
  fi
done < tmpfile

if [[ $insidecodeblock == true ]]; then
  echo "\`\`\`"
  insidecodeblock=false
fi
