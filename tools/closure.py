#!/usr/bin/env python3
"""What in a module is reachable from the symbols the rest of the tree uses?

deadcode.py answers a narrower question -- reachability inside one .cc, with
the matching header's inline bodies bolted on as extra roots. That misses the
case where a header inline and a .cc function keep *each other* alive while
nothing outside the pair calls either:

    adjutant_matrix.h:  HEADER_INLINE RowMajorMatrixMultiply(...)  -> calls
    adjutant_matrix.cc: ColMajorMatrixMultiply(...)

Neither is reachable, but each is referenced, so per-file analysis keeps both.
adjutant_matrix was 33 definitions of which exactly 3 were reachable; the
per-TU pass reported "0 dead".

This treats the .cc and .h as one unit, takes as roots only those names that
some *other* file mentions, and walks the reference graph from there. Comments
and string literals are blanked first, so prose cannot make a symbol look used.

Caveat: a root here is any mention outside the pair, including a mention from
code that is itself dead. Prune leaf-first, or re-run after each removal --
adjutant_common's GetMajIdx only became visible as dead once ComputeMajAlleles
was gone.

Macro bodies count as references. BigstackAllocX is reached only through
`#define BIGSTACK_ALLOC_X(...) plink2::BigstackAllocX<tt>(...)`, so a walk over
function bodies alone declares it dead and the build then fails in every file
that uses the macro.

Usage:
  closure.py <module>            e.g. closure.py adjutant_compress_stream
  closure.py <module> --names    just the dead names, comma-separated, for
                                 feeding straight into prune_named.py
  closure.py <module> --root X,Y extra roots. Needed for adjutant itself: no
                                 other file mentions main(), so without it the
                                 whole program looks dead.
"""
import glob
import re
import sys

sys.path.insert(0, __file__.rsplit('/', 1)[0])
from deadcode import blank_noncode, index_functions


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__)
    mod = sys.argv[1]
    want_names = '--names' in sys.argv

    pair = [f for f in (mod + '.cc', mod + '.h') if glob.glob(f)]
    if not pair:
        sys.exit('no such module: %s' % mod)

    # Every span per name, not just the first. plink2_bits.cc defines
    # PopcountWordsIntersect twice, once per #ifdef arm; indexing only the
    # first meant the walk scanned the AVX2 body and never the non-AVX2 one,
    # so the callee that only the second arm uses looked unreferenced.
    defs, srcs = {}, {}
    for f in pair:
        src = open(f).read()
        code = blank_noncode(src)
        srcs[f] = (src, code)
        for name, (s, e) in index_functions(src, code).items():
            defs.setdefault(name, []).append((f, s, e))

    others = []
    for f in sorted(glob.glob('adjutant*.cc') + glob.glob('adjutant*.h') +
                    glob.glob('include/*.cc') + glob.glob('include/*.h')):
        if f in pair:
            continue
        others.append(blank_noncode(open(f, errors='replace').read()))
    others = '\n'.join(others)

    roots = [n for n in defs if re.search(r'\b%s\b' % re.escape(n), others)]

    # Anything named in a #define inside the pair is a root: the macro is the
    # real call site, and it lives outside every function body the walk visits.
    # `#  define` -- hash, spaces, keyword -- is legal and pgenlib uses it, so
    # match the directive rather than a literal '#define' prefix.
    directive = re.compile(r'^\s*#\s*define\b')
    for f in pair:
        for line in srcs[f][1].split('\n'):
            if directive.match(line):
                roots += [t for t in re.findall(r'\b([A-Za-z_]\w*)\b', line)
                          if t in defs]

    if '--root' in sys.argv:
        roots += [r for r in sys.argv[sys.argv.index('--root') + 1].split(',')
                  if r in defs]

    keep, stack = set(), list(roots)
    while stack:
        n = stack.pop()
        if n in keep:
            continue
        keep.add(n)
        for f, s, e in defs[n]:
            for tok in set(re.findall(r'\b([A-Za-z_]\w*)\b', srcs[f][1][s:e])):
                if tok in defs and tok not in keep:
                    stack.append(tok)

    dead = sorted(set(defs) - keep)

    # Safety net, same as deadcode.py's: index_functions does not recognise
    # every signature style, so a caller can go unindexed and its body never
    # scanned. Before calling anything dead, confirm the name really does not
    # appear outside its own definition spans.
    survivors = []
    for d in list(dead):
        spans = sorted(defs[d], key=lambda x: (x[0], x[1]))
        outside = []
        for f in pair:
            code = srcs[f][1]
            cuts, prev = [], 0
            for sf, s, e in spans:
                if sf != f:
                    continue
                cuts.append(code[prev:s])
                prev = e
            cuts.append(code[prev:])
            outside.append(''.join(cuts))
        if re.search(r'\b%s\b' % re.escape(d), '\n'.join(outside)):
            survivors.append(d)
            dead.remove(d)
    if survivors and not want_names:
        print('  kept %d referenced-but-unreachable: %s'
              % (len(survivors), ' '.join(survivors)))
    if want_names:
        print(','.join(dead))
        return 0

    print('%s: %d definitions, %d reachable, %d dead'
          % (mod, len(defs), len(keep), len(dead)))
    print('  roots (referenced elsewhere in the tree): %s'
          % (' '.join(sorted(roots)) or '(none -- whole module is dead)'))
    total = 0
    for d in dead:
        n = sum(srcs[f][0].count('\n', s, e) for f, s, e in defs[d])
        total += n
        where = ' '.join(sorted(set(f for f, _, _ in defs[d])))
        print('    %-44s %-26s %4d lines' % (d, where, n))
    print('  %d dead lines' % total)
    return 0


if __name__ == '__main__':
    sys.exit(main())
