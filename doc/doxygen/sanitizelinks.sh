# sanitizes the links to the given lists, because doxygen somehow links
# to a page with a wrong index
sanitizelinks () {
  NEW_FILE=`grep -l "\"title\">$1" html/*html | egrep -o [0-9]+`
  OLD_FILE=`awk -v a=$NEW_FILE 'BEGIN {printf("%05d", a-1)}'`
  sed -i "s/$OLD_FILE/$NEW_FILE/g" html/*html
}

sanitizelinks "Todo List"
sanitizelinks "Bug List"
sanitizelinks "Deprecated List"
sanitizelinks "Bibliography"
