set(TEX_INPUTS
  0_dumux-handbook.tex
  0_discussion.tex
  0_listingstyle.tex
  1_introduction.tex
  2_detailedinstall.tex
  2_quickinstall.tex
  2_quickstartguide.tex
  3_tutorial.tex
  3_tutorialcoupled.tex
  3_tutorialdecoupled.tex
  3_furtherpractice.tex
  4_guidelines.tex
  4_developingdumux.tex
  4_externaltools.tex
  4_newfoldersetup.tex
  4_newtoninanutshell.tex
  4_parameterfiles.tex
  4_restartsimulations.tex
  4_structure.tex
  5_stepsofasimulation.tex
  5_fluidframework.tex
  5_grids.tex
  5_models.tex
  5_propertysystem.tex
  5_spatialdiscretizations.tex
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

dune_add_latex_document(0_dumux-handbook.tex
  FATHER_TARGET doc
  DEFAULT_PDF
  BIBFILES dumux-handbook.bib
  INPUTS ${TEX_INPUTS}
  IMAGES ${TEX_IMAGES})

create_doc_install(${CMAKE_CURRENT_BINARY_DIR}/dumux-handbook.pdf
  ${CMAKE_INSTALL_DOCDIR})
