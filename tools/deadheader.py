#!/usr/bin/env python3
"""Remove unreferenced inline function definitions from a header.

deadcode.py handles top-level functions in a .cc.  Headers are different: they
carry HEADER_INLINE / static inline bodies that no .o file records a symbol for,
so the linker view used for .cc pruning says nothing about them.

A header inline is dead when its name occurs nowhere in the project's code
outside its own definition.  Comments and string literals are blanked first, so
a mention in prose does not keep it alive.

Usage:
  deadheader.py report <file.h>
  deadheader.py strip  <file.h>
"""
import glob
import re
import sys

sys.path.insert(0, __file__.rsplit('/', 1)[0])
import tools_strip as T
from deadcode import blank_noncode, index_functions


def project_code(exclude):
    """Stripped code of every source file except the one being analysed."""
    parts = []
    for f in sorted(glob.glob('adjutant*.cc') + glob.glob('adjutant*.h') +
                    glob.glob('include/*.cc') + glob.glob('include/*.h')):
        if f == exclude:
            continue
        parts.append(blank_noncode(open(f).read()))
    return '\n'.join(parts)


def main():
    # A mistyped mode must not reach the strip path below: it rewrites the
    # header in place with no backup.  See deadcode.py for the incident.
    if len(sys.argv) < 3 or sys.argv[1] not in ('report', 'strip'):
        sys.exit(__doc__)
    mode, path = sys.argv[1], sys.argv[2]
    src = open(path).read()
    code = blank_noncode(src)
    defs = index_functions(src, code)
    outside_files = project_code(path)

    dead = []
    for name, (s, e) in defs.items():
        # references elsewhere in this same header (outside the definition)
        self_rest = code[:s] + code[e:]
        if re.search(r'\b%s\b' % re.escape(name), self_rest):
            continue
        if re.search(r'\b%s\b' % re.escape(name), outside_files):
            continue
        dead.append(name)

    if mode == 'report':
        total = sum(src.count('\n', *defs[d]) for d in dead)
        print('%s: %d inline defs, %d unreferenced (~%d lines)'
              % (path, len(defs), len(dead), total))
        for d in sorted(dead)[:40]:
            print('    %s' % d)
        return 0

    base = T.balanced(src)
    removed = 0
    for d in sorted(dead, key=lambda x: -defs[x][0]):
        s, e = defs[d]
        cand = src[:s].rstrip('\n') + '\n' + src[e:].lstrip('\n')
        if T.balanced(cand) == base:
            src = cand
            removed += 1
    open(path, 'w').write(src)
    print('%s: removed %d/%d unreferenced inline defs' % (path, removed, len(dead)))
    return 0


if __name__ == '__main__':
    sys.exit(main())
