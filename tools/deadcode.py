#!/usr/bin/env python3
"""Accurate dead-function elimination for a single translation unit.

Earlier attempts failed in two opposite ways:

  * matching only call syntax `name(` missed function *pointers* -- thread entry
    points passed to SetThreadFuncAndData() looked dead and got deleted;
  * matching bare identifiers `name` kept everything alive, because a mention
    in a comment counted as a reference.

The fix is to blank comments and string/char literals first (preserving byte
offsets so line numbers stay valid), then match bare identifiers on what's left.
That catches function pointers without being fooled by prose.

Usage:
  deadcode.py report <file.cc> <root,root,...>
  deadcode.py strip  <file.cc> <root,root,...>   [--header <file.h>]
"""
import re
import sys

sys.path.insert(0, __file__.rsplit('/', 1)[0])
import tools_strip as T

KEYWORDS = {
    'if', 'for', 'while', 'switch', 'return', 'sizeof', 'else', 'do', 'case',
    'break', 'continue', 'goto', 'typedef', 'struct', 'union', 'enum', 'static',
    'const', 'inline', 'unlikely', 'likely', 'assert', 'defined', 'alignof',
    'and', 'or', 'not', 'new', 'delete', 'throw', 'catch', 'try', 'namespace',
}


def blank_noncode(src):
    """Replace comments and literals with spaces, keeping length identical."""
    out = list(src)
    i = 0
    n = len(src)
    while i < n:
        c = src[i]
        if c == '/' and i + 1 < n and src[i + 1] == '/':
            j = src.find('\n', i)
            j = n if j < 0 else j
            for k in range(i, j):
                out[k] = ' '
            i = j
        elif c == '/' and i + 1 < n and src[i + 1] == '*':
            j = src.find('*/', i + 2)
            j = n if j < 0 else j + 2
            for k in range(i, j):
                if src[k] != '\n':
                    out[k] = ' '
            i = j
        elif c in '"\'':
            j = i + 1
            while j < n:
                if src[j] == '\\':
                    j += 2
                    continue
                if src[j] == c:
                    j += 1
                    break
                j += 1
            for k in range(i, min(j, n)):
                if src[k] != '\n':
                    out[k] = ' '
            i = j
        else:
            i += 1
    return ''.join(out)


def index_functions(src, code):
    """name -> (start, end) for every top-level function definition."""
    defs = {}
    pat = re.compile(
        r'^(?:static\s+)?(?:THREAD_FUNC_DECL\s+)?[A-Za-z_][\w:*&<>,\s]*?\b([A-Za-z_]\w*)\s*\(',
        re.M)
    for m in pat.finditer(code):
        name = m.group(1)
        if name in KEYWORDS:
            continue
        try:
            op = code.index('(', m.start())
            cp = T.match_paren(code, op)
        except ValueError:
            continue
        after = code[cp + 1:cp + 200].lstrip()
        if after.startswith(';'):
            continue                      # declaration, not definition
        br = code.find('{', cp)
        if br < 0 or code[cp + 1:br].strip() not in ('', 'const'):
            continue
        try:
            end = T.match_brace(code, br) + 1
        except ValueError:
            continue
        defs.setdefault(name, (m.start(), end))
    return defs


def reachable(code, defs, roots):
    names = set(defs)
    keep = set()
    stack = [r for r in roots if r in defs]
    while stack:
        n = stack.pop()
        if n in keep:
            continue
        keep.add(n)
        s, e = defs[n]
        for tok in set(re.findall(r'\b([A-Za-z_]\w*)\b', code[s:e])):
            if tok in names and tok not in keep:
                stack.append(tok)
    return keep


def main():
    # An unrecognised mode used to fall through to the strip path below, which
    # rewrites the file in place with no backup and no confirmation:
    # `deadcode.py list adjutant_pgs_eval.cc ""` reported "removed 5/5 dead
    # functions" and took the file from 3054 lines to 1280.  Validate first.
    if len(sys.argv) < 4 or sys.argv[1] not in ('report', 'strip'):
        sys.exit(__doc__)
    mode, path = sys.argv[1], sys.argv[2]
    roots = [r for r in sys.argv[3].split(',') if r]
    src = open(path).read()
    code = blank_noncode(src)
    assert len(code) == len(src)
    defs = index_functions(src, code)

    # The matching header often carries inline/static functions that call into
    # this .cc.  Anything they reference is a root, or we would delete symbols
    # the header still needs.
    hdr_path = None
    if '--header' in sys.argv:
        hdr_path = sys.argv[sys.argv.index('--header') + 1]
    elif path.endswith('.cc') and __import__('os').path.exists(path[:-3] + '.h'):
        hdr_path = path[:-3] + '.h'
    if hdr_path:
        hsrc = open(hdr_path).read()
        hcode = blank_noncode(hsrc)
        # Only the header's *inline bodies* and macros count.  A bare prototype
        # says the symbol exists, not that anything calls it -- treating
        # prototypes as roots would keep every function alive.
        hbodies = []
        for _n, (hs, he) in index_functions(hsrc, hcode).items():
            hbodies.append(hcode[hs:he])
        for line in hcode.split('\n'):
            if line.lstrip().startswith('#define'):
                hbodies.append(line)
        for tok in set(re.findall(r'\b([A-Za-z_]\w*)\b', '\n'.join(hbodies))):
            if tok in defs:
                roots.append(tok)

    keep = reachable(code, defs, roots)
    dead = sorted(set(defs) - keep)

    # Safety net.  index_functions() does not recognise every signature style
    # (macro-prefixed returns, some multi-line forms), so a caller can go
    # unindexed and its body never scanned.  Before deleting anything, confirm
    # the name genuinely does not occur outside its own definition.
    survivors = []
    for d in list(dead):
        s0, e0 = defs[d]
        outside = code[:s0] + code[e0:]
        if re.search(r'\b%s\b' % re.escape(d), outside):
            survivors.append(d)
            dead.remove(d)
    if survivors and mode == 'report':
        print('    (kept %d referenced-but-unreachable: %s)'
              % (len(survivors), ' '.join(survivors[:6])))

    if mode == 'report':
        print('%s: %d functions, %d reachable, %d dead'
              % (path, len(defs), len(keep), len(dead)))
        for d in dead:
            s, e = defs[d]
            print('    %-44s %d lines' % (d, src.count('\n', s, e)))
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
    print('%s: removed %d/%d dead functions' % (path, removed, len(dead)))

    if '--header' in sys.argv:
        hdr = sys.argv[sys.argv.index('--header') + 1]
        h = open(hdr).read()
        n = 0
        for d in dead:
            h2 = re.sub(r'^[A-Za-z][^\n;]*\b%s\([^;]*;\n' % re.escape(d), '', h, flags=re.M)
            if h2 != h:
                h = h2
                n += 1
        open(hdr, 'w').write(h)
        print('%s: removed %d declarations' % (hdr, n))
    return 0


if __name__ == '__main__':
    sys.exit(main())
