#!/bin/sh

# 1. obtain a list new_parameters.csv of all current parameters
# retrieve all occurrences of GET_PARAM and GET_RUNTIME_PARAM
find dumux/ -name '*.[ch][ch]' -exec grep 'GET_PARAM' {} \; | sort -u >new_parameters.csv
find dumux/ -name '*.[ch][ch]' -exec grep 'GET_RUNTIME_PARAM' {} \; | sort -u >>new_parameters.csv
# remove #define's
sed -i '/#define/d' new_parameters.csv
# remove everything before GET_PARAM and GET_RUNTIME_PARAM
sed -i 's/^.*GET_PARAM/GET_PARAM/' new_parameters.csv
sed -i 's/^.*GET_RUNTIME_PARAM/GET_RUNTIME_PARAM/' new_parameters.csv
# remove everything before and including first comma
awk '{print substr($0,index($0,",")+1)}' new_parameters.csv >tmp.txt
mv tmp.txt new_parameters.csv
# remove leading whitespace
sed -i 's/^[ \t]*//' new_parameters.csv
# remove everything after and including paranthesis
sed -i 's/).*//' new_parameters.csv
# sort uniquely
sort -u new_parameters.csv -o new_parameters.csv
# keep only the lines containing two commas
sed -i '/,.*,/!d' new_parameters.csv
# append types to the end of the lines
sed -i '/^bool,/ s/$/, bool/' new_parameters.csv
sed -i '/^double,/ s/$/, double/' new_parameters.csv
sed -i '/^int,/ s/$/, int/' new_parameters.csv
sed -i '/^Scalar,/ s/$/, Scalar/' new_parameters.csv
sed -i '/^std::string,/ s/$/, std::string/' new_parameters.csv
# remove types from the beginning of the lines
sed -i 's/bool,//' new_parameters.csv
sed -i 's/double,//' new_parameters.csv
sed -i 's/int,//' new_parameters.csv
sed -i 's/Scalar,//' new_parameters.csv
sed -i 's/std::string,//' new_parameters.csv
# remove all blanks
sed -i 's/ //g' new_parameters.csv
# add blanks after every comma
sed -i 's/,/, /g' new_parameters.csv
# sort uniquely
sort -u new_parameters.csv -o new_parameters.csv
# remove lines containing no parameter names
sed -i '/, , /d' new_parameters.csv
#adapt to doxygen format
sed -i 's/^/ * | /' new_parameters.csv
sed -i 's/, / | /g' new_parameters.csv

# 2. obtain a list old_parameters.csv of all old parameters
cp -v doc/doxygen/extradoc/parameterlist.txt old_parameters.csv
# truncate the test from file
sed -i -e '1,12d' old_parameters.csv
# remove lines containing no parameter names
sed -i '/| |  | /d' old_parameters.csv
# remove bold format for parameter groups
sed -i 's/ \\b / /g' old_parameters.csv
# loop over all lines, add group name if necessary
word=""
while read line; do
  word=$(echo "$line" | awk '{print $3;}')
  #echo "$word"
  firstletter=$(echo $word | head -c 1)
  if [ "$firstletter" != "|" ]; then
   group=$word
  fi
  if [ "$group" = "" ]; then
   echo "No group found in line $line."
  fi
  linewithoutgroup=$(echo "$line" | cut -d " " -f2-)
  echo "$group $linewithoutgroup" >>tmp.csv
done <old_parameters.csv
mv tmp.csv old_parameters.csv
# truncate the default and description information
cat old_parameters.csv | while read line
do
  TEMP=`echo "${line% |*}"`
  echo "${TEMP% |*}" >> tmp.csv
done
mv tmp.csv old_parameters.csv
#adapt to doxygen format
sed -i 's/| | /| /g' old_parameters.csv
sed -i 's/^/ * | /' old_parameters.csv

# 3. compare lists of old and new parameters
# treat additions
diff -u old_parameters.csv new_parameters.csv | grep -E "^\+" | grep -v "\++" >added.csv
if [[ -s added.csv ]]; then
  echo ""
  echo "Compared to parameters.ods, the following parameters have been _added_:"
  echo ""
  echo "Group | Parameter | Type:"
  echo "-------------------------"
  sed -i 's/+//g' added.csv
  sed -i 's/ \* | //g' added.csv
  cat added.csv
  echo ""
  echo "Search for those parameters and their default values. If default values"
  echo "are set, you can find them by searching for GroupNameParameterName (one"
  echo "word only). Decide whether the parameter should be included in"
  echo "doc/doxygen/extradoc/parameters.ods, possibly by discussing with the responsible guy."
fi
# treat deletions
diff -u old_parameters.csv new_parameters.csv | grep -E "^\-" | grep -v "\--" >deleted.csv
if [[ -s deleted.csv ]]; then
  echo ""
  echo "Compared to parameters.ods, the following parameters have been _deleted_:"
  echo ""
  echo "Group | Parameter | Type:"
  echo "-------------------------"
  sed -i 's/-//g' deleted.csv
  sed -i 's/ \* | //g' deleted.csv
  sed -i '/^\*/d' deleted.csv
  cat deleted.csv
  echo ""
  echo "Check whether those parameters really don't exist anymore by grepping for"
  echo "their name. If so, delete the corresponding rows in doc/doxygen/extradoc/parameters.ods."
fi
# final remark
if [[ ! -s added.csv && ! -s deleted.csv ]]; then
  echo ""
  echo "Compared to parameters.ods, no parameters have been added or deleted."
  echo ""
else
  echo ""
  echo "_Important:_ If you change parameters.ods, be sure to export it to"
  echo "parameters.html."
fi
# clean up
rm -f old_parameters.csv new_parameters.csv added.csv deleted.csv
