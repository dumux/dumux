# sanitizes the links to the given lists, because doxygen somehow links
# to a page with a wrong index
sanitizelinks () {
  NEW_FILE=`grep -l "\"title\">$1" html/*html | egrep -o [0-9]+`
  OLD_FILE=`awk -v a=$NEW_FILE 'BEGIN {printf("%05d", a-1)}'`
#   echo $1: $OLD_FILE $NEW_FILE
  sed -i "s#$OLD_FILE#$NEW_FILE#g" html/a*html
  sed -i "s# $1 # <a href=\"a$NEW_FILE.html\">$1</a> #g" html/index.html
  sed -i "s# $1,# <a href=\"a$NEW_FILE.html\">$1</a>,#g" html/index.html
}

sanitizelinks "Deprecated List" # has to be called before Todo List
sanitizelinks "Todo List" # has to be called before Warning List
sanitizelinks "Warning List"
sanitizelinks "Bug List"
sanitizelinks "Bibliography"
