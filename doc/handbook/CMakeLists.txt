# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later

# use latexmk
dune_add_latex_document(
  SOURCE 0_dumux-handbook.tex
  FATHER_TARGET doc
  INSTALL ${CMAKE_INSTALL_DOCDIR}/doc/handbook)
