set(TEX_INPUTS
  designpatterns.tex
  DumuxFlow.tex
  fluidframework.tex
  getting-started.tex
  guidelines.tex
  install.tex
  intro.tex
  models.tex
  newFolder.tex
  dumux-handbook.tex
  NewtonInANutshell.tex
  parameterTree.tex
  propertysystem.tex
  quick-install.tex
  quickstart-guide.tex
  spatialdiscretization.tex
  structure.tex
  TipsNTricks.tex
  tutorial-coupled.tex
  tutorial-decoupled.tex
  tutorial-newmodel.tex
  tutorial.tex
  ModelDescriptions/1p2cimplicitmodel.tex
  ModelDescriptions/1pdecoupledmodel.tex
  ModelDescriptions/1pimplicitmodel.tex
  ModelDescriptions/2p2cdecoupledpressuremodel.tex
  ModelDescriptions/2p2cdecoupledtransportmodel.tex
  ModelDescriptions/2p2cimplicitmodel.tex
  ModelDescriptions/2pdecoupledpressuremodel.tex
  ModelDescriptions/2pdecoupledsaturationmodel.tex
  ModelDescriptions/2pdfmimplicitmodel.tex
  ModelDescriptions/2pimplicitmodel.tex
  ModelDescriptions/3p3cimplicitmodel.tex
  ModelDescriptions/3pimplicitmodel.tex
  ModelDescriptions/co2implicitmodel.tex
  ModelDescriptions/el1p2cmodel.tex
  ModelDescriptions/el2pmodel.tex
  ModelDescriptions/elasticmodel.tex
  ModelDescriptions/mpncimplicitmodel.tex
  ModelDescriptions/nonisothermalimplicitmodel.tex
  ModelDescriptions/richardsimplicitmodel.tex
  ModelDescriptions/stokesimplicitmodel.tex
  ModelDescriptions/stokesncimplicitmodel.tex
  ModelDescriptions/stokesncniimplicitmodel.tex
  ../../tutorial/tutorial_coupled.cc
  ../../tutorial/tutorial_coupled.input
  ../../tutorial/tutorialproblem_coupled.hh
  ../../tutorial/tutorialspatialparams_coupled.hh
  ../../tutorial/tutorial_decoupled.cc
  ../../tutorial/tutorial_decoupled.input
  ../../tutorial/tutorialproblem_decoupled.hh
  ../../tutorial/tutorialspatialparams_decoupled.hh)

set(TEX_IMAGES
  PNG/box_disc.png
  PNG/cc_disc.png
  PNG/dumux_strucutre_flowchart_horizontal_explained.png
  PNG/dunedesign.png
  ../logo/dumux_logo_hires_whitebg.png)

dune_add_latex_document(dumux-handbook.tex
  FATHER_TARGET doc
  DEFAULT_PDF
  BIBFILES dumux-handbook.bib
  INPUTS ${TEX_INPUTS}
  IMAGES ${TEX_IMAGES})

create_doc_install(${CMAKE_CURRENT_BINARY_DIR}/dumux-handbook.pdf
  ${CMAKE_INSTALL_DOCDIR})
