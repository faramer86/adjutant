#!/usr/bin/env python3
"""Iteratively delete code that references now-removed commands.

For each compiler error, delete the smallest enclosing construct:
  * if the error sits in an if/while condition -> delete that whole block
  * otherwise -> delete the single statement containing it

Brace balance is re-checked after every deletion; anything that would unbalance
the file is skipped and reported instead of being applied blindly.

Usage:
  autofix.py <file.cc>          # rewrites <file.cc> in place
"""
import re
import subprocess
import sys

import tools_strip as T

CXX = ['c++', '-c', '-O2', '-std=c++17', '-g', '-DZSTD_MULTITHREAD',
       '-Ilibdeflate', '-DACCELERATE_NEW_LAPACK', '-ferror-limit=0',
       '-fsyntax-only']

# Never delete a span touching these: they are the --score load/compute path.
# Without this guard the fixer happily removes e.g. the LoadPvar() call when a
# neighbouring statement is what actually failed to compile.
PROTECTED = (
    'LoadPvar', 'LoadPsam', 'LoadPhenos', 'ScoreReport',
    'LoadAlleleAndGenoCounts', 'ComputeAlleleFreqs', 'ReadAlleleFreqs',
    'ApplyVariantBpFilters', 'AllocAndPopulateIdHtableMt', 'KeepOrRemove',
    'ExtractExcludeFlagNorange', 'SortedXidboxInitAlloc',
    'PgfiInitPhase1', 'PgfiInitPhase2', 'PgrInit', 'CleanupPgfi', 'CleanupPgr',
    'InitScore', 'CleanupScore', 'Plink2Core', 'CmdlineParsePhase',
    'SetChrCodes', 'InitChrInfo', 'FinalizeChrset', 'CleanupChrInfo',
    'SampleSortFileMap', 'UpdateSampleSubsets', 'KeepRemoveIf',
    'EnforceFreqConstraints', 'EnforceGenoThresh', 'EnforceHweThresh',
    'EnforceMinBpSpace', 'RandomThinProb', 'CleanupLogfile', 'InitLogfile',
)


def protected(chunk):
    return any(name in chunk for name in PROTECTED)


def err_lines(path):
    out = subprocess.run(CXX + [path], capture_output=True, text=True).stderr
    seen = set()
    for l in out.split('\n'):
        m = re.match(r'^%s:(\d+):' % re.escape(path), l)
        if m and ': error:' in l:
            seen.add(int(m.group(1)))
    return sorted(seen)


def stmt_span(lines, i):
    """Line range [a,b] of the statement containing line index i (0-based)."""
    a = i
    # walk back to the start of the statement
    while a > 0:
        prev = lines[a - 1].rstrip()
        if prev.endswith((';', '{', '}', ':')) or not prev:
            break
        a -= 1
    b = i
    while b < len(lines) - 1 and not lines[b].rstrip().endswith((';', '{', '}')):
        b += 1
    return a, b


def block_span(lines, i):
    """If line i is inside an if/while header, return the whole block's span."""
    a, b = stmt_span(lines, i)
    head = '\n'.join(lines[a:b + 1]).strip()
    if not re.match(r'^\}?\s*(else\s+)?(if|while)\s*\(', head):
        return None
    text = '\n'.join(lines)
    off = sum(len(x) + 1 for x in lines[:a])
    try:
        op = text.index('(', off)
        cp = T.match_paren(text, op)
        br = text.index('{', cp)
        if text[cp + 1:br].strip():
            return None
        end = T.match_brace(text, br)
    except ValueError:
        return None
    endline = text[:end].count('\n')
    # absorb a trailing 'else' chain? no - keep it simple, bail if present
    after = text[end + 1:end + 40].lstrip()
    if after.startswith('else'):
        return None
    return a, endline


def main():
    # The target used to default to adjutant.cc, so a bare `autofix.py` was
    # already the destructive invocation: 30 rounds of compiler-guided deletion
    # written back in place, with no backup and no confirmation.  Name the file.
    # A leading `-` is rejected too, so the sibling script's `--apply` lands on
    # the usage string rather than on open('--apply').
    if len(sys.argv) != 2 or sys.argv[1].startswith('-'):
        sys.exit(__doc__)
    path = sys.argv[1]
    baseline = T.balanced(open(path).read())
    total = 0
    for rnd in range(1, 30):
        errs = err_lines(path)
        if not errs:
            print('round %d: clean' % rnd)
            break
        src = open(path).read()
        lines = src.split('\n')
        # Work bottom-up so earlier line numbers stay valid, and apply one
        # deletion at a time so a single bad span can be skipped rather than
        # poisoning the whole round.
        removed = skipped = 0
        for e in sorted(errs, reverse=True):
            i = e - 1
            if i >= len(lines):
                continue
            span = block_span(lines, i) or stmt_span(lines, i)
            if not span:
                skipped += 1
                continue
            if protected('\n'.join(lines[span[0]:span[1] + 1])):
                skipped += 1
                continue
            cand = lines[:span[0]] + lines[span[1] + 1:]
            if T.balanced('\n'.join(cand)) != baseline:
                skipped += 1
                continue
            removed += span[1] - span[0] + 1
            lines = cand
        if not removed:
            print('round %d: nothing removable, %d error lines remain'
                  % (rnd, len(errs)))
            break
        open(path, 'w').write('\n'.join(lines))
        total += removed
        print('round %d: removed %d lines (%d spans skipped), %d error lines'
              % (rnd, removed, skipped, len(errs)))
    print('total lines removed: %d' % total)
    print('remaining errors: %d' % len(err_lines(path)))


if __name__ == '__main__':
    main()
