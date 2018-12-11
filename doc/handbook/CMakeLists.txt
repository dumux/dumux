set(TEX_INPUTS
  0_dumux-handbook.tex
  0_listingstyle.tex
  1_introduction.tex
  2_detailedinstall.tex
  2_quickinstall.tex
  3_course.tex
  3_furtherpractice.tex
  4_assemblinglinearsystem.tex
  4_developingdumux.tex
  4_externaltools.tex
  4_newfoldersetup.tex
  4_parameterfiles.tex
  4_restartsimulations.tex
  4_structure.tex
  5_inputoutput.tex
  5_models.tex
  5_parallel.tex
  5_propertysystem.tex
  5_spatialdiscretizations.tex
  5_stepsofasimulation.tex
  installDumux.sh)

set(TEX_IMAGES
  png/box_disc.png
  png/cc_disc.png
  png/dumux_strucutre_flowchart_horizontal_explained.png
  png/dunedesign.png
  ../logo/dumux_logo_hires_whitebg.png
  png/dalton1.png
  png/dalton2.png
  pdf/staggered_grid.pdf)

dune_add_latex_document(0_dumux-handbook.tex
  BIBFILES dumux-handbook.bib
  INPUTS ${TEX_INPUTS}
  IMAGES ${TEX_IMAGES})

create_doc_install(${CMAKE_CURRENT_BINARY_DIR}/0_dumux-handbook.pdf
  ${CMAKE_INSTALL_DOCDIR} 0_dumux-handbook)
