#!/usr/bin/env python3
import sys
import argparse
import re

ignorable = " \t\n\r"
TEMPLATE = """
- name: {name}
  in:
    {_in}
  out:
    {out}
"""

class EOF(Exception):
    pass

def ignore_comments(s):
    if peek(s) == "%":
        while True:
            c = s.read(1)
            if c == "\n":
                return ignore_comments(s)
            elif c == '':
                raise EOF()

def peek(s):
    c = s.read(1)
    s.seek(s.tell() -1)
    return c

def ignore(s, chars):
    while True:
        c = s.read(1)
        if c not in chars:
            s.seek(s.tell() -1)
            break
        elif c == '':
            raise EOF()

def ignore_all(s):
    ignore(s, ignorable)
    ignore_comments(s)
    ignore(s, ignorable)

def string(s, ignoring):
    result = ""
    while True:
        c = s.read(1)
        if c and c in ignoring:
            s.seek(s.tell() -1)
            break
        elif c == '':
            raise EOF()
        else:
            result += c
    return result

def parse_name(s):
    ignore_all(s)
    return string(s, ignorable + "[]()")

def expect(s, c):
    cc = s.read(1)
    if cc != c:
        s.seek(s.tell() -1)
        raise Exception("Expecting %s but got %s" % (c, cc))

def parse_argument(s):
    ignore_all(s)
    if peek(s) == "(":
        expect(s, "(")
        name = parse_name(s)
        ignore_all(s)
        value = parse_name(s)
        ignore_all(s)
        if re.match(r"^[a-zA-Z]+", value):
            value = "${}".format(value)
        ignore_all(s)
        expect(s, ")")
    else:
        if peek(s) == "]":
            return None
        name = parse_name(s)
        value = "${}".format(name)
    ignore_all(s)
    return (name, value)

def parse_arguments(s):
    args = []
    ignore_all(s)
    expect(s, "[")
    ignore_all(s)
    while True:
        arg = parse_argument(s)
        if arg is None:
            break
        args += [(arg)]
    ignore_all(s)
    expect(s, "]")
    return args

def parse_step(s):
    ignore_all(s)
    _name = parse_name(s)
    ignore_all(s)
    _in = parse_arguments(s)
    ignore_all(s)
    _out = parse_arguments(s)
    ignore_all(s)
    expect(s, ".")
    return (_name, _in, _out)

def emit_step(step):
    _in = ("{}"
            if not step[1]
            else "\n    ".join("{}: {}".format(*k) for k in step[1]))
    out = ("{}"
            if not step[2]
            else "\n    ".join("{}: {}".format(*k) for k in step[2]))
    name = step[0]
    return TEMPLATE.format(**locals())

if __name__ == "__main__":
    while True:
        try:
            print(emit_step(parse_step(sys.stdin)))
        except EOF:
            print("# end of automatic translation")
            break
