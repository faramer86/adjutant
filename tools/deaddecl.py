#!/usr/bin/env python3
"""Prototypes left behind after their definitions were removed.

Deleting a function body and its declaration is two edits, and the second one
gets missed: a `.h` keeps `void ReflectFmatrix(...);` for a function that no
longer exists anywhere. Nothing complains -- a declaration for an undefined
function is legal C++ right up until someone calls it -- so these accumulate
silently and make the header describe an API the binary does not have.

adjutant_matrix.h had thirteen: SvdRect, ExtractEigvecs, TransposeMultiplySelfIncr
and friends from earlier rounds, plus InvertMatrixChecked from this one.

A declaration is dead when no definition of that name exists anywhere in the
tree AND nothing references the name outside the declaration itself. Comments
and string literals are blanked first.

Usage:
  deaddecl.py report [file.h ...]     default: every adjutant*.h and include/*.h
  deaddecl.py strip  [file.h ...]
"""
import glob
import re
import sys

sys.path.insert(0, __file__.rsplit('/', 1)[0])
from deadcode import blank_noncode, index_functions

KEYWORDS = {'sizeof', 'if', 'for', 'while', 'switch', 'return', 'defined',
            'alignof', 'static_assert', 'assert', 'typedef', 'operator'}

DECL = re.compile(r'^[A-Za-z_][^\n;{]*?\b([A-Za-z_]\w*)\s*\([^;{]*\)\s*;[ \t]*\n',
                  re.M)


def all_sources():
    return sorted(glob.glob('adjutant*.cc') + glob.glob('adjutant*.h') +
                  glob.glob('include/*.cc') + glob.glob('include/*.h'))


def main():
    # `report` stays the default for a bare invocation, but anything that is
    # not a documented mode must not reach the else-branch below, which
    # rewrites headers in place.  See deadcode.py for the incident.
    mode = sys.argv[1] if len(sys.argv) > 1 else 'report'
    if mode not in ('report', 'strip'):
        sys.exit(__doc__)
    targets = sys.argv[2:] or [f for f in all_sources() if f.endswith('.h')]

    code = {f: blank_noncode(open(f, errors='replace').read())
            for f in all_sources()}
    defined = set()
    for f, c in code.items():
        defined |= set(index_functions(open(f, errors='replace').read(), c))

    total = 0
    for path in targets:
        src = open(path).read()
        c = code[path]
        dead = []
        for m in DECL.finditer(c):
            name = m.group(1)
            if name in KEYWORDS or name in defined:
                continue
            # any mention of the name other than this declaration?
            used = False
            for f, cc in code.items():
                for mm in re.finditer(r'\b%s\b' % re.escape(name), cc):
                    if f == path and m.start() <= mm.start() < m.end():
                        continue
                    used = True
                    break
                if used:
                    break
            if not used:
                dead.append((m.start(), m.end(), name))

        if not dead:
            continue
        total += len(dead)
        if mode == 'report':
            print('%s: %d dead declarations' % (path, len(dead)))
            for s, _e, n in dead:
                print('    line %-5d %s' % (c.count('\n', 0, s) + 1, n))
        else:
            for s, e, _n in sorted(dead, reverse=True):
                src = src[:s] + src[e:]
            open(path, 'w').write(src)
            print('%s: removed %d dead declarations' % (path, len(dead)))
    if mode == 'report':
        print('%d dead declarations in total' % total)
    return 0


if __name__ == '__main__':
    sys.exit(main())
