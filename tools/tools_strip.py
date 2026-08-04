#!/usr/bin/env python3
"""Surgical block removal for the plink2 --score-only reduction.

adjutant.cc's flag parser and command dispatcher are built from long chains of
brace-delimited if/else-if branches.  Editing those by hand across ~30 commands
is error-prone, so this does it mechanically:

  * find a branch by a substring of its condition
  * brace-match to its true end (string/char/comment aware)
  * delete it, repairing the chain so the remaining branches stay well-formed

Usage:
  tools_strip.py flag  <file> <cond-substr> [<cond-substr>...]
  tools_strip.py func  <file> <function-name> [<function-name>...]
  tools_strip.py check <file>
"""
import re
import sys


def scan(src, i):
    """Advance past a string literal, char literal or comment starting at i.

    Returns the index just past it, or None if src[i] starts none of those.
    """
    c = src[i]
    if c == '"' or c == "'":
        j = i + 1
        while j < len(src):
            if src[j] == '\\':
                j += 2
                continue
            if src[j] == c:
                return j + 1
            j += 1
        return j
    if src.startswith('//', i):
        j = src.find('\n', i)
        return len(src) if j < 0 else j
    if src.startswith('/*', i):
        j = src.find('*/', i + 2)
        return len(src) if j < 0 else j + 2
    return None


def match_brace(src, open_idx):
    """Given the index of a '{', return the index of its matching '}'."""
    depth = 0
    i = open_idx
    while i < len(src):
        nxt = scan(src, i)
        if nxt is not None:
            i = nxt
            continue
        if src[i] == '{':
            depth += 1
        elif src[i] == '}':
            depth -= 1
            if depth == 0:
                return i
        i += 1
    raise ValueError('unbalanced braces from offset %d' % open_idx)


def find_branch(src, needle):
    """Locate the if/else-if branch whose header is/contains `needle`.

    `needle` must start at the branch's own 'if (' keyword -- e.g.
    'if (strequal_k_unsafe(flagname_p2, "xport"))'.  Anchoring on the header
    itself (rather than scanning backwards for an enclosing 'if') is what keeps
    this from silently selecting the *previous* branch.

    Returns (kind, start, open_paren) where kind is 'if' or 'elif' and start is
    the index of the owning 'if'/'else' keyword.
    """
    pos = src.find(needle)
    if pos < 0:
        return None
    if src.count(needle) > 1:
        raise ValueError('ambiguous needle (%d matches): %.60s'
                         % (src.count(needle), needle))
    if not src.startswith('if', pos):
        raise ValueError('needle must begin at its own "if (": %.60s' % needle)
    open_paren = src.index('(', pos)
    head = src[:pos].rstrip()
    if head.endswith('else'):
        return ('elif', len(head) - 4, open_paren)
    return ('if', pos, open_paren)


def cut_branch(src, needle):
    found = find_branch(src, needle)
    if not found:
        return src, False
    kind, start, open_paren = found
    close_paren = match_paren(src, open_paren)
    brace = src.find('{', close_paren)
    if brace < 0:
        return src, False
    end = match_brace(src, brace) + 1

    tail = src[end:]
    stripped = tail.lstrip()
    lead = len(tail) - len(stripped)

    if kind == 'elif':
        # "} else if (X) {...}" -> drop it; preceding '}' still closes the
        # previous branch and the following 'else if' reattaches to it.
        return src[:start] + tail[lead:], True

    # First branch of a chain: promote the next 'else if' to 'if',
    # or drop a trailing bare 'else'.
    if stripped.startswith('else if'):
        return src[:start] + 'if' + stripped[len('else if'):], True
    if stripped.startswith('else'):
        rest = stripped[len('else'):].lstrip()
        if rest.startswith('{'):
            inner_end = match_brace(rest, 0)
            return src[:start] + rest[1:inner_end].strip() + rest[inner_end + 1:], True
    return src[:start] + tail[lead:], True


def match_paren(src, open_idx):
    depth = 0
    i = open_idx
    while i < len(src):
        nxt = scan(src, i)
        if nxt is not None:
            i = nxt
            continue
        if src[i] == '(':
            depth += 1
        elif src[i] == ')':
            depth -= 1
            if depth == 0:
                return i
        i += 1
    raise ValueError('unbalanced parens from offset %d' % open_idx)


def cut_function(src, name):
    """Delete a whole top-level function definition by name."""
    pat = re.compile(
        r'^(?:[A-Za-z_][\w:<>,*&\s]*?\s)?' + re.escape(name) + r'\s*\(',
        re.M)
    for m in pat.finditer(src):
        open_paren = src.index('(', m.start() + len(m.group(0)) - 1)
        close_paren = match_paren(src, open_paren)
        after = src[close_paren + 1:]
        head = after[:200].lstrip()
        if head.startswith(';'):
            continue  # forward declaration, not the definition
        brace = src.find('{', close_paren)
        if brace < 0 or src[close_paren + 1:brace].strip() not in ('', 'const'):
            continue
        end = match_brace(src, brace) + 1
        start = m.start()
        # Absorb an immediately preceding comment block / THREAD_FUNC_DECL line.
        return src[:start].rstrip('\n') + '\n' + src[end:].lstrip('\n'), True
    return src, False


def balanced(src):
    depth = 0
    i = 0
    while i < len(src):
        nxt = scan(src, i)
        if nxt is not None:
            i = nxt
            continue
        if src[i] == '{':
            depth += 1
        elif src[i] == '}':
            depth -= 1
        i += 1
    return depth


def main():
    # `cut_branch if mode == 'flag' else cut_function` below treats every
    # unrecognised mode as `func` and writes the result, so an unchecked typo
    # silently deletes functions.  Validate against the documented set first.
    if len(sys.argv) < 3 or sys.argv[1] not in ('flag', 'func', 'check'):
        sys.exit(__doc__)
    mode, path = sys.argv[1], sys.argv[2]
    src = open(path).read()

    if mode == 'check':
        # Some files carry a nonzero baseline because mutually exclusive
        # #if/#else branches each open a brace (e.g. RealMain vs main in
        # adjutant.cc).  The compiler sees one; a plain scanner sees both.
        # What matters is that the depth does not *change*, so the expected
        # baseline is passed in.
        expected = int(sys.argv[3]) if len(sys.argv) > 3 else 0
        d = balanced(src)
        ok = d == expected
        print('%s: brace depth %d (expected %d) %s'
              % (path, d, expected, 'OK' if ok else 'MISMATCH'))
        return 0 if ok else 1

    before_depth = balanced(src)
    cut = cut_branch if mode == 'flag' else cut_function
    done, missing = [], []
    for needle in sys.argv[3:]:
        src, ok = cut(src, needle)
        (done if ok else missing).append(needle)

    after_depth = balanced(src)
    if after_depth != before_depth:
        sys.exit('ABORT: brace depth changed %d -> %d' % (before_depth, after_depth))

    open(path, 'w').write(src)
    print('%s: removed %d, not found %d' % (path, len(done), len(missing)))
    if missing:
        print('  not found: ' + ', '.join(missing))
    return 0


if __name__ == '__main__':
    sys.exit(main())
