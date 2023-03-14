#!/usr/bin/env perl
# SPDX-FileCopyrightInfo: Copyright © DuMux Project contributors, see AUTHORS.md in root folder
# SPDX-License-Identifier: GPL-3.0-or-later
# SPDX-FileCopyrightInfo: Copyright 2017 - 2023 Simulating eXtreme Spacetimes Collaboration
# SPDX-License-Identifier: MIT
# Based on the original by Simulating eXtreme Spacetimes Collaboration,
# licensed under MIT. All Changes are licensed under GPL-3.0-or-later.
#
#
# Doxygen filter to format markdown math for Doxygen
#
# In Doxygen documentation, including markdown files or Jupyter notebooks that
# get converted to markdown, you can use standard $...$ syntax for inline math,
# and $$...$$ or ```math...``` syntax for display-style math.
#
# This is needed until Doxygen adds support for markdown-style math, see
# https://github.com/doxygen/doxygen/issues/8009.
#
# Based on the original by Simulating eXtreme Spacetimes Collaboration
# distributed under the MIT License as specified below:
#
# ==============================================================================
# SpECTRE Release License
# ==============================================================================
#
# Copyright 2017 - 2023 Simulating eXtreme Spacetimes Collaboration
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.
#
# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.
#
# This is added to the original to allow for $$...$$ for display math
# Wrap inline math in \f$...\f$
# - Match code blocks (wrapped in ``` or \code...\endcode) and inline code
#   (wrapped in ` or ``), and print them directly to the output with no changes.
# - Match display math blocks (wrapped in \f[...\f], \f{...\f}, or
#   \begin...\end), and print them directly to the output with no changes.
# - Replace $$...$$ by \f{equation}{...\f}.
# - The '?' makes the pattern match lazily, i.e., match as few characters as
#   possible.
# - Modifiers:
#     s: Allow '.' to match newlines
#     g: Replace all occurrences
#     e: Evaluate replacement as perl code, so we can switch between
#        replacements depending on the matched alternative, using the '//'
#        (definedness) operator.
#     x: Ignore spaces for better readability
s{ (```.*?```)
   | (``.*?``)
   | (`.*?`)
   | (\\code.*?\\endcode)
   | (\\f\[.*?\\f\])
   | (\\f\{.*?\\f\})
   | (\\begin\{(.*?)\}(.*?)\\end\{\6\})
   | \${2}(.*?)\${2}
}{ $1 // $2 // $3 // $4 // $5 // $6 // $7 // "\\begin\{align\}$10\\end{align}" }sgex;

# This is added to the original to allow for $`...`$ for inline math (GitLab flavor)
# Wrap inline math in \f$...\f$
# - Match code blocks (wrapped in ``` or \code...\endcode) and inline code
#   (wrapped in ` or ``), and print them directly to the output with no changes.
# - Match display math blocks (wrapped in \f[...\f], \f{...\f}, or
#   \begin...\end), and print them directly to the output with no changes.
# - Replace $...$ by \f$...\f$, unless already preceded by '\f$'.
s{ (```.*?```)
   | (``.*?``)
   | (`.*?`)
   | (\\code.*?\\endcode)
   | (\\f\[.*?\\f\])
   | (\\f\{.*?\\f\})
   | (\\begin\{(.*?)\}(.*?)\\end\{\6\})
   | (?<!\\f)\$`(.*?)`\$
}{ $1 // $2 // $3 // $4 // $5 // $6 // $7 // "\\f\$$10\\f\$" }sgex;

# Wrap inline math in \f$...\f$
# - Match code blocks (wrapped in ``` or \code...\endcode) and inline code
#   (wrapped in ` or ``), and print them directly to the output with no changes.
# - Match display math blocks (wrapped in \f[...\f], \f{...\f}, or
#   \begin...\end), and print them directly to the output with no changes.
# - Replace $...$ by \f$...\f$, unless already preceded by '\f$'.
s{ (```.*?```)
   | (``.*?``)
   | (`.*?`)
   | (\\code.*?\\endcode)
   | (\\f\[.*?\\f\])
   | (\\f\{.*?\\f\})
   | (\\begin\{(.*?)\}(.*?)\\end\{\6\})
   | (?<!\\f)\$(.*?)\$
}{ $1 // $2 // $3 // $4 // $5 // $6 // $7 // "\\f\$$10\\f\$" }sgex;

# This is added to the original to allow for ```math...``` for display math (GitLab flavor)
# Wrap display math in \f{equation}{ ... \f}
# - Replace ```math...``` with \f{equation}{...\f}.
# - Match code blocks (wrapped in ``` or \code...\endcode) and Doxygen-style
#   display equations formatted either \f[...\f] or \f{...}{...\f}, and print
#   them to the output with no changes.
s{ ```math(.*?)```
   | (```.*?```)
   | (\\code.*?\\endcode)
   | (\\f\[.*?\\f\])
   | (\\f\{.*?\\f\})
}{ $2 // "\\f\{align\}\{$1\\f\}" // $3 // $4 // $5 }sgex;

# Wrap display math in \f{equation}{ ... \f}
# - Match code blocks (wrapped in ``` or \code...\endcode) and Doxygen-style
#   display equations formatted either \f[...\f] or \f{...}{...\f}, and print
#   them to the output with no changes.
# - Replace \begin{...}...\end{...} with \f{...}{...\f}.
s{ (```.*?```)
   | (\\code.*?\\endcode)
   | (\\f\[.*?\\f\])
   | (\\f\{.*?\\f\})
   | \\begin\{(.*?)\}(.*?)\\end\{\5\}
}{ $1 // $2 // $3 // $4 // "\\f\{$5\}\{$6\\f\}" }sgex;
