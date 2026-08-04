#!/usr/bin/env python3
"""Remove unreferenced typedef'd structs, enums and flagsets.

deadcode.py and deadheader.py only look at functions.  Removing a command also
strands its *types* -- the thread-context struct it passed around, the flagset
its modifiers set.  Those keep other dead symbols alive by referencing them, so
they have to go too.

A type is dead when its name appears nowhere outside its own definition.

The count-based rule this used to apply -- "dead if the name appears at most
twice in the project, once for the tag and once for the closing `} Name;`" --
is wrong for an *anonymous* typedef:

    typedef union {          <- no tag, so the name appears ONCE, not twice
      CmpExprK k;
      ...
    } CmpExprArgs;

which leaves the threshold covering one genuine use. CmpExprArgs is used
exactly once, in adjutant_cmdline.h itself, and the old rule deleted it.
Excising the definition span and searching what remains has no such blind spot.

Usage:
  deadtype.py report <file>
  deadtype.py strip  <file>
"""
import glob
import re
import sys

sys.path.insert(0, __file__.rsplit('/', 1)[0])
import tools_strip as T
from deadcode import blank_noncode


def project_code(exclude=None):
    parts = []
    for f in sorted(glob.glob('adjutant*.cc') + glob.glob('adjutant*.h') +
                    glob.glob('include/*.cc') + glob.glob('include/*.h')):
        if f == exclude:
            continue
        parts.append(blank_noncode(open(f).read()))
    return '\n'.join(parts)


def enum_members(code, s, e):
    """Identifiers declared between the braces of an enum/flagset definition.

    A type name can be unreferenced while its constants are load-bearing:

        typedef enum {          <- 'ChrOffset' appears nowhere else...
          kChrOffsetX,
          ...
        } ChrOffset;            <- ...but kChrOffsetX is used all over

    Same for FLAGSET_DEF, whose members are named kf<Name>* rather than after
    the type. Deleting the definition takes the constants with it.
    """
    try:
        br = code.index('{', s)
        close = T.match_brace(code, br)
    except ValueError:
        return []
    body = code[br + 1:close]
    out, depth, cur = [], 0, []
    for ch in body:
        if ch in '([{':
            depth += 1
        elif ch in ')]}':
            depth -= 1
        if ch == ',' and depth == 0:
            out.append(''.join(cur))
            cur = []
        else:
            cur.append(ch)
    out.append(''.join(cur))
    names = []
    for part in out:
        m = re.match(r'\s*([A-Za-z_]\w*)', part.split('=')[0])
        if m:
            names.append(m.group(1))
    return names


def find_types(src, code):
    """name -> (start, end) covering the whole typedef/flagset definition."""
    out = {}

    # typedef struct Foo { ... } Foo;   /  typedef enum { ... } Foo;
    for m in re.finditer(r'\btypedef\s+(?:struct|union|enum)\b', code):
        try:
            br = code.index('{', m.start())
            close = T.match_brace(code, br)
            semi = code.index(';', close)
        except ValueError:
            continue
        tail = code[close + 1:semi].strip()
        if re.fullmatch(r'\w+', tail or ''):
            out.setdefault(tail, (m.start(), semi + 1))

    # FLAGSET_DEF_START() ... FLAGSET_DEF_END(Foo);
    for m in re.finditer(r'FLAGSET_DEF_END\((\w+)\)\s*;', code):
        name = m.group(1)
        start = code.rfind('FLAGSET_DEF_START', 0, m.start())
        if start >= 0:
            out.setdefault(name, (start, m.end()))
    for m in re.finditer(r'ENUM_U31_DEF_END\((\w+)\)\s*;', code):
        name = m.group(1)
        start = code.rfind('ENUM_U31_DEF_START', 0, m.start())
        if start >= 0:
            out.setdefault(name, (start, m.end()))
    return out


def main():
    # A mistyped mode must not reach the strip path below: it rewrites the
    # file in place with no backup.  See deadcode.py for the incident.
    if len(sys.argv) < 3 or sys.argv[1] not in ('report', 'strip'):
        sys.exit(__doc__)
    mode, path = sys.argv[1], sys.argv[2]
    src = open(path).read()
    code = blank_noncode(src)
    types = find_types(src, code)
    othercode = project_code(exclude=path)

    dead = []
    for name, (s, e) in types.items():
        outside = code[:s] + code[e:]
        # The type name, plus -- for enums and flagsets -- every constant the
        # definition declares. Any one of them still in use keeps the lot.
        probes = [name]
        if re.search(r'\b(?:typedef\s+enum|FLAGSET_DEF_START|ENUM_U31_DEF_START)',
                     code[s:e]):
            probes += enum_members(code, s, e)
        alive = False
        for p in probes:
            pat = r'\b%s\b' % re.escape(p)
            if re.search(pat, outside) or re.search(pat, othercode):
                alive = True
                break
        if not alive:
            dead.append(name)

    if mode == 'report':
        total = sum(src.count('\n', *types[d]) for d in dead)
        print('%s: %d types, %d unreferenced (~%d lines)'
              % (path, len(types), len(dead), total))
        for d in sorted(dead):
            print('    %s' % d)
        return 0

    base = T.balanced(src)
    removed = 0
    for d in sorted(dead, key=lambda x: -types[x][0]):
        s, e = types[d]
        cand = src[:s].rstrip('\n') + '\n' + src[e:].lstrip('\n')
        if T.balanced(cand) == base:
            src = cand
            removed += 1
    open(path, 'w').write(src)
    print('%s: removed %d/%d unreferenced types' % (path, removed, len(dead)))
    return 0


if __name__ == '__main__':
    sys.exit(main())
