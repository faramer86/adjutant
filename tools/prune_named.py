#!/usr/bin/env python3
"""Remove a *named* list of top-level functions from a .cc, and their
prototypes from the matching .h.

deadcode.py decides for itself what is dead, from reachability within one
translation unit.  That misses two things:

  * a cluster of dead functions that call *each other* looks alive to it --
    adjutant_matrix.cc's matrix-inversion family is 13 functions that form
    exactly such a ring;
  * a function whose only call site sits in a branch the optimiser proves
    unreachable (`if (randmem)` where randmem is a literal 0) looks alive,
    because the call really is there in the source.

Both are found by linkmap.py instead, which asks the linker.  This script is
the other half: it takes that verdict as an explicit list and applies it.

Usage:
  prune_named.py <file.cc> <name,name,...> [--header <file.h>] [--dry-run]
"""
import re
import sys

sys.path.insert(0, __file__.rsplit('/', 1)[0])
import deadcode as D
import tools_strip as T


COND = re.compile(r'^\s*#\s*(if|ifdef|ifndef|endif)\b', re.M)


def cond_balanced(span):
    """True if #if/#endif pair up inside `span`.

    Brace matching does not understand the preprocessor. plink2_bits.cc defines
    PopcountVecsNoAvx2Intersect and its AVX2 sibling in adjacent #ifdef arms; a
    span that starts in one arm can run past #else into the next, so deleting
    one function silently takes the other with it. The build then fails on a
    symbol nobody asked to remove. Refuse those spans.
    """
    depth = 0
    for m in COND.finditer(span):
        depth += -1 if m.group(1) == 'endif' else 1
        if depth < 0:
            return False
    return depth == 0


def strip_header_decls(hsrc, names):
    """Drop prototypes for `names`. Only whole single- or multi-line
    declarations ending in ';' -- never an inline body, which would leave
    a dangling block behind."""
    out, removed = hsrc, []
    for n in names:
        # A declaration runs from a line-start return type to the first ';'
        # with no intervening '{'.  The [^{;]* guard is what keeps this from
        # eating an inline definition.
        pat = re.compile(
            r'^[A-Za-z_][^\n;{]*\b%s\s*\([^;{]*\)\s*;[ \t]*\n' % re.escape(n),
            re.M | re.S)
        new = pat.sub('', out)
        if new != out:
            out, _ = new, removed.append(n)
    return out, removed


def main():
    if len(sys.argv) < 3:
        sys.exit(__doc__)
    path = sys.argv[1]
    names = [n for n in sys.argv[2].split(',') if n]
    hdr = None
    if '--header' in sys.argv:
        hdr = sys.argv[sys.argv.index('--header') + 1]
    elif path.endswith('.cc'):
        import os
        cand = path[:-3] + '.h'
        if os.path.exists(cand):
            hdr = cand
    dry = '--dry-run' in sys.argv

    src = open(path).read()
    code = D.blank_noncode(src)
    assert len(code) == len(src)
    defs = D.index_functions(src, code)

    missing = [n for n in names if n not in defs]
    if missing:
        print('  !! not found as top-level definitions in %s: %s'
              % (path, ' '.join(missing)))

    base = T.balanced(src)
    removed, lines = [], 0
    # Delete back-to-front so earlier offsets stay valid.
    for n in sorted((n for n in names if n in defs), key=lambda x: -defs[x][0]):
        s, e = defs[n]
        if not cond_balanced(src[s:e]):
            print('  !! refusing %s -- span straddles #if/#endif' % n)
            continue
        cand = src[:s].rstrip('\n') + '\n' + src[e:].lstrip('\n')
        if T.balanced(cand) != base:
            print('  !! refusing %s -- deletion would unbalance braces' % n)
            continue
        lines += src.count('\n', s, e)
        src = cand
        removed.append(n)

    hremoved = []
    if hdr:
        hsrc = open(hdr).read()
        hnew, hremoved = strip_header_decls(hsrc, removed)

    if dry:
        print('%s: would remove %d/%d functions (%d lines)'
              % (path, len(removed), len(names), lines))
        if hdr:
            print('%s: would remove %d declarations' % (hdr, len(hremoved)))
        return 0

    open(path, 'w').write(src)
    print('%s: removed %d/%d functions (%d lines)'
          % (path, len(removed), len(names), lines))
    if hdr:
        open(hdr, 'w').write(hnew)
        print('%s: removed %d declarations' % (hdr, len(hremoved)))
    left = [n for n in removed if n not in hremoved]
    if hdr and left:
        print('  (no prototype found for: %s)' % ' '.join(left))
    return 0


if __name__ == '__main__':
    sys.exit(main())
