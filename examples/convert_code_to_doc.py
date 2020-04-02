#!/usr/bin/env python3

"""
A simple documentation generator
generating Markdown-documented examples
from annotated source code
"""
from pyparsing import *

# tell pyparsing to never ignore white space characters
ParserElement.setDefaultWhitespaceChars('')

def parseTaggedContent(keyname, action, open="[[", close="]]", endTag="/"):
    """
    Match content between [[keyname]] and [[/keyname]] and apply the action on it
    """
    start = LineStart() + ZeroOrMore(" ") + "//" + ZeroOrMore(" ") + (open + str(keyname) + close) + LineEnd()
    end = LineStart() + ZeroOrMore(" ") + "//" + ZeroOrMore(" ") + (open + endTag + str(keyname) + close) + LineEnd()
    return start.suppress() + SkipTo(end).setParseAction(action) + end.suppress()

def createMarkdownCode(markDownToken):
    """
    Put code into Markdown syntax for code blocks with syntax highlighting
    """
    def action(token):
        # only print code snippets with content
        if not token[0].rstrip():
            return ""
        else:
            return "\n```" + markDownToken + "\n" + token[0].rstrip() + "\n```\n\n"
    return action

def cppRules():
    """
    Define a list of rules to apply for cpp source code
    """
    suppressHeader = Suppress(Combine("// -*-" + SkipTo("*******/" + LineEnd(), include=True)))
    suppressHeaderGuard = Suppress("#ifndef" + Optional(restOfLine) + LineEnd() + "#define" + Optional(restOfLine))
    suppressEndHeaderGuard = Suppress("#endif" + Optional(restOfLine))

    # make a code block (possibly containing comments) between [[codeblock]] and [[/codeblock]]
    createCppBlock = createMarkdownCode("cpp")
    parseCodeblock = parseTaggedContent("codeblock", action=createCppBlock)

    # treat doc and code line
    parseDoc = LineStart() + Suppress(ZeroOrMore(" ") + "//" + ZeroOrMore(" ")) + Optional(restOfLine)
    parseCode = LineStart() + ~(ZeroOrMore(" ") + "//") + (SkipTo(parseDoc) | SkipTo(StringEnd()))
    parseCode.setParseAction(createCppBlock)
    docTransforms = parseCodeblock | parseDoc | parseCode

    return [suppressHeader, suppressHeaderGuard, suppressEndHeaderGuard, docTransforms]

def transformCode(code, rules, codeFileName):

    # exclude stuff between [[exclude]] and [[/exclude]]
    exclude = parseTaggedContent("exclude", action=replaceWith(""))
    code = exclude.transformString(code)

    # Enable toggling content between [[content]] and [[/content]]
    def wrapContentIntoDetails(token):
        beginDetails = "//\n// <details open>\n"
        summmary = "// <summary><b>Click to hide/show the file documentation</b> (or inspect the [source code]({}))</summary>\n//\n".format(codeFileName)
        endDetails = "\n//\n// </details>\n"
        return beginDetails + summmary + token[0] + endDetails
    wrapContent = parseTaggedContent("content", action=wrapContentIntoDetails)
    code = wrapContent.transformString(code)

    # Transform "[[details]] content" and "[[/details]]" to HTML
    transformDetailsBegin = LineStart() + Suppress(ZeroOrMore(" ") + "//" + ZeroOrMore(" ") + "[[details]]" + ZeroOrMore(" ")) + Optional(restOfLine)
    def detailBeginHTML(token):
        return "// <details><summary> Click to show " + token[0] + "</summary>\n"
    transformDetailsBegin.setParseAction(detailBeginHTML)
    code = transformDetailsBegin.transformString(code)
    transformDetailsEnd = LineStart() + Suppress(ZeroOrMore(" ") + "//" + ZeroOrMore(" ") + "[[/details]]") + Optional(restOfLine)
    transformDetailsEnd.setParseAction(replaceWith("// </details>\n"))
    code = transformDetailsEnd.transformString(code)

    for transform in rules:
        code = transform.transformString(code)
    return code
