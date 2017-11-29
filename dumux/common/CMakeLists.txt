#install headers
install(FILES
basicproperties.hh
boundaryconditions.hh
boundarytypes.hh
dimensionlessnumbers.hh
eigenvalues.hh
exceptions.hh
fixedlengthspline_.hh
intersectionmapper.hh
math.hh
parameters.hh
propertysystem.hh
quad.hh
reorderingdofmapper.hh
splinecommon_.hh
spline.hh
start.hh
tabulated2dfunction.hh
timemanager.hh
valgrind.hh
variablelengthspline_.hh
DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}/dumux/common)
