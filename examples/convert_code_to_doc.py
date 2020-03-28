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
    start = LineStart() + ZeroOrMore(" ") + Literal("//") + ZeroOrMore(" ") + Literal(open + str(keyname) + close)
    end = LineStart() + ZeroOrMore(" ") + Literal("//") + ZeroOrMore(" ") + Literal(open + endTag + str(keyname) + close)
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
            return "```" + markDownToken + "\n" + token[0].rstrip() + "\n```\n"
    return action

def cppRules():
    """
    Define a list of rules to apply for cpp source code
    """
    header = Suppress(Combine(Literal("// -*-") + SkipTo(Literal("*******/") + LineEnd()) + Literal("*******/")))
    headerGuard = Suppress(Literal("#ifndef") + Optional(restOfLine) + LineEnd() + Literal("#define") + Optional(restOfLine))
    endHeaderGuard = Suppress(Literal("#endif") + Optional(restOfLine))

    # exclude stuff between [[exclude]] and [[/exclude]]
    exclude = parseTaggedContent("exclude", action=replaceWith(""))

    # make a code block (possibly containing comments) between [[codeblock]] and [[/codeblock]]
    action = createMarkdownCode("cpp")
    codeblock = parseTaggedContent("codeblock", action=action)

    # treat doc and code line
    doc = LineStart() + Suppress(ZeroOrMore(" ") + Literal("//") + ZeroOrMore(" ")) + Optional(restOfLine)
    code = LineStart() + ~(ZeroOrMore(" ") + Literal("//")) + (SkipTo(doc) | SkipTo(StringEnd()))
    code.setParseAction(action)
    docTransforms = codeblock | doc | code

    return [header, headerGuard, endHeaderGuard, exclude, docTransforms]

def transformCode(code, rules):
    for transform in rules:
        code = transform.transformString(code)
    return code
