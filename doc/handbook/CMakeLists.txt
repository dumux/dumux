if(${DUNE_COMMON_VERSION} VERSION_GREATER_EQUAL "2.7")
# use latexmk
dune_add_latex_document(
  SOURCE 0_dumux-handbook.tex
  FATHER_TARGET doc
  INSTALL ${CMAKE_INSTALL_DOCDIR}/doc/handbook)

else()
# use UseLATEX.cmake
set(TEX_INPUTS
  0_dumux-handbook.tex
  0_listingstyle.tex
  1_introduction.tex
  2_quickinstall.tex
  3_detailedinstall.tex
  4_course.tex
  4_furtherpractice.tex
  5_assemblinglinearsystem.tex
  5_developingdumux.tex
  5_externaltools.tex
  5_newfoldersetup.tex
  5_parameterfiles.tex
  5_restartsimulations.tex
  5_structure.tex
  6_basics.tex
  6_inputoutput.tex
  6_parallel.tex
  6_propertysystem.tex
  6_spatialdiscretizations.tex
  6_stepsofasimulation.tex
  6_temporaldiscretizations.tex
  installDumux.sh)

set(TEX_IMAGES
  png/box_disc.png
  png/cctpfa.png
  png/dumux_strucutre_flowchart_horizontal_explained.png
  png/dunedesign.png
  ../logo/dumux_logo_hires_whitebg.png
  png/dalton1.png
  png/dalton2.png
  pdf/staggered_grid.pdf
  pdf/mpfa_iv.pdf)

dune_add_latex_document(0_dumux-handbook.tex
  BIBFILES dumux-handbook.bib
  INPUTS ${TEX_INPUTS}
  IMAGES ${TEX_IMAGES})

create_doc_install(${CMAKE_CURRENT_BINARY_DIR}/0_dumux-handbook.pdf
  ${CMAKE_INSTALL_DOCDIR} 0_dumux-handbook)
endif()
