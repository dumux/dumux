#! /bin/bash

if [ $# -eq 0 ]
  then
    echo "No arguments supplied."
    echo "Usage: bash ./replace_property_macros.sh file1 [file2 ...]"
fi

for TMP in $@; do

echo "File: $TMP"

sed -i \
"
s/typename[ ]*GET_PROP_TYPE[ ]*(\([^,]*\),[ ]*\([^)]*\))::/typename GetPropType<\1, Properties::\2>::/g;
s/typename[ ]*GET_PROP_TYPE[ ]*(\([^,]*\),[ ]*\([^)]*\))/GetPropType<\1, Properties::\2>/g;
s/GET_PROP_TYPE[ ]*(\([^,]*\),[ ]*\([^)]*\))/GetPropType<\1, Properties::\2>/g;
s/NEW_TYPE_TAG[ ]*(\([^,]*\),[ ]*INHERITS_FROM[ ]*(\([^)]*\)).*/struct TTag::\1 { using InheritsFrom = std::tuple<\2>; };/g;
s/using InheritsFrom = std::tuple<\([^,]*\),[ ]*\([^>]*\)>/using InheritsFrom = std::tuple<\2, \1>/g;
" \
$TMP

gawk -i inplace '/struct TTag::/ && !x {print "// Create new type tags\nnamespace TTag {"; x=1} 1' $TMP
gawk -i inplace '/struct TTag::/{seen++} seen && !/struct TTag::/{print "} // end namespace TTag"; seen=0} 1' $TMP

sed -i \
"
s/struct TTag::/struct /g;
s/GET_PROP_VALUE[ ]*(\([^,]*\),[ ]*\([^)]*\))/getPropValue<\1, Properties::\2>()/g;
s/SET_BOOL_PROP[ ]*(\([^,]*\),[ ]*\([^,]*\),[ ]*\([^)]*\))/template<class TypeTag>\nstruct \2<TypeTag, TTag::\1> { static constexpr bool value = \3; }/g;
s/SET_INT_PROP[ ]*(\([^,]*\),[ ]*\([^,]*\),[ ]*\([^)]*\))/template<class TypeTag>\nstruct \2<TypeTag, TTag::\1> { static constexpr int value = \3; }/g;
s/SET_TYPE_PROP[ ]*(\([^,]*\),[ ]*\([^,]*\),[ ]*\([^)]*\))/template<class TypeTag>\nstruct \2<TypeTag, TTag::\1> { using type = \3; }/g;
s/SET_PROP[ ]*(\([^,]*\),[ ]*\([^)]*\))/template<class TypeTag>\nstruct \2<TypeTag, TTag::\1>/g;
s/NEW_PROP_TAG[ ]*(\([^)]*\))/template<class TypeTag, class MyTypeTag>\nstruct \1 { using type = UndefinedProperty; }/g;
s/TTAG[ ]*(\([^)]*\))/Properties::TTag::\1/g;
" \
$TMP

done

echo "Property macros have been removed."
echo "Manual tweaking might be necessary, especially if"
echo ""
echo "- NEW_TYPE_TAG uses have not been one line after another. In this case,"
echo "  . remove superfluous lines \"} // end namespace TTag\" and/or"
echo "  . add additional lines \"namespace TTag {\"."
echo ""
echo "- Usages of SET_TYPE_PROP or other macros extend over more than one line."
echo "  In this case, replace the usages manually."
echo ""
echo "- Macros have been used outside of the namespace Dumux and without a"
echo "  corresponding alias. In this case, prepend things with a \"Dumux::\"."
