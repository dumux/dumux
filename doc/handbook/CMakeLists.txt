set(TEX_INPUTS
  0_dumux-handbook.tex
  0_listingstyle.tex
  1_introduction.tex
  2_detailedinstall.tex
  2_quickinstall.tex
  2_quickstartguide.tex
  3_tutorial.tex
  3_tutorialimplicit.tex
  3_tutorialsequential.tex
  3_furtherpractice.tex
  4_assemblinglinearsystem.tex
  4_guidelines.tex
  4_developingdumux.tex
  4_externaltools.tex
  4_newfoldersetup.tex
  4_parameterfiles.tex
  4_restartsimulations.tex
  4_structure.tex
  5_stepsofasimulation.tex
  5_fluidframework.tex
  5_grids.tex
  5_models.tex
  5_propertysystem.tex
  5_spatialdiscretizations.tex
  ../../tutorial/tutorial_implicit.cc
  ../../tutorial/tutorial_implicit.input
  ../../tutorial/tutorialproblem_implicit.hh
  ../../tutorial/tutorialspatialparams_implicit.hh
  ../../tutorial/tutorial_sequential.cc
  ../../tutorial/tutorial_sequential.input
  ../../tutorial/tutorialproblem_sequential.hh
  ../../tutorial/tutorialspatialparams_sequential.hh)

set(TEX_IMAGES
  PNG/box_disc.png
  PNG/cc_disc.png
  PNG/dumux_strucutre_flowchart_horizontal_explained.png
  PNG/dunedesign.png
  ../logo/dumux_logo_hires_whitebg.png)

if(DUNE_COMMON_VERSION VERSION_LESS 2.5)
  # this only works before Dune 2.5
  dune_add_latex_document(0_dumux-handbook.tex
    FATHER_TARGET doc
    DEFAULT_PDF
    BIBFILES dumux-handbook.bib
    INPUTS ${TEX_INPUTS}
    IMAGES ${TEX_IMAGES})

  create_doc_install(${CMAKE_CURRENT_BINARY_DIR}/0_dumux-handbook.pdf
    ${CMAKE_INSTALL_DOCDIR})
else()
  # the Dune 2.5 way of calling these macros
  dune_add_latex_document(0_dumux-handbook.tex
    BIBFILES dumux-handbook.bib
    INPUTS ${TEX_INPUTS}
    IMAGES ${TEX_IMAGES})

  create_doc_install(${CMAKE_CURRENT_BINARY_DIR}/0_dumux-handbook.pdf
    ${CMAKE_INSTALL_DOCDIR} 0_dumux-handbook)
endif()
