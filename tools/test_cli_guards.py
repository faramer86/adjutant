#!/usr/bin/env python3
"""An unrecognised mode argument must never rewrite the source file.

`deadcode.py list adjutant_pgs_eval.cc ""` -- a typo, the modes are `report`
and `strip` -- printed "removed 5/5 dead functions" and took a 3054-line file
down to 1280.  There is no backup and no confirmation; recovery was by hand.

Every one of these scripts dispatched on mode with the destructive path as the
*fallback* (`if mode == 'report': ... return` then strip below, or an `else:`,
or `cut_branch if mode == 'flag' else cut_function`), so any string that was
not the read-only mode meant "write".  linkmap.py alone got it right, with
`else: sys.exit(__doc__)`.

autofix.py and reduce_plink2.py have no mode to mistype, and were left out of
that fix for exactly that reason.  They wrote `adjutant.cc` anyway -- one via
`sys.argv[1] if len(sys.argv) > 1 else 'adjutant.cc'`, the other by hardcoding
the name -- so `python3 tools/autofix.py` with no arguments at all *was* the
destructive invocation.  Both now want it spelled out: a path for autofix,
`--apply` for reduce_plink2.

Run:  python3 tools/test_cli_guards.py
"""
import os
import shutil
import subprocess
import sys
import tempfile

TOOLS = os.path.dirname(os.path.abspath(__file__))

FIXTURE_CC = '''\
#include "adjutant_fixture.h"

typedef struct FixtureLonelyStructStruct {
  int a;
} FixtureLonelyStruct;

static int FixtureDeadHelper(int x) {
  return x + 1;
}

int FixtureLiveEntry(int x) {
  return x * 2;
}
'''

FIXTURE_H = '''\
#ifndef ADJUTANT_FIXTURE_H
#define ADJUTANT_FIXTURE_H

int FixtureLiveEntry(int x);

void FixtureVanishedDecl(int x);

static inline int FixtureDeadInline(int x) {
  return x - 1;
}

#endif
'''

# autofix.py is driven by compiler errors: it deletes the smallest enclosing
# statement of every error line.  No #include, so the fixture stands alone --
# an unresolved header is itself an error, on line 1, and the deletion that
# follows is not what we mean to measure.
FIXTURE_AUTOFIX = '''\
int FixtureLiveEntry(int x) {
  FixtureVanishedCall(x);
  return x * 2;
}
'''

# reduce_plink2.py's first pass drops includes of deleted headers and the
# Plink2Cmdline fields of removed commands.  It indexes on the struct markers
# rather than searching for them, so both have to be here or it raises before
# reaching the write.
FIXTURE_REDUCE = '''\
#include "plink2_glm.h"

typedef struct Plink2CmdlineStruct {
  GlmInfo glm_info;
  double score_thresh;
} Plink2Cmdline;
'''

BOGUS = 'list'   # the actual typo that caused the loss

# (script, target, bogus-mode argv, read-only argv, destructive argv).
# tools_strip.py's read-only mode takes a different tail from its destructive
# one -- `check <file> [depth]` wants an integer, not a function name.
CC = 'adjutant_fixture.cc'
H = 'adjutant_fixture.h'
CASES = [
    ('deadcode.py',    CC, [BOGUS, CC, ''],  ['report', CC, ''],  ['strip', CC, '']),
    ('deadheader.py',  H,  [BOGUS, H],       ['report', H],       ['strip', H]),
    ('deadtype.py',    CC, [BOGUS, CC],      ['report', CC],      ['strip', CC]),
    ('deaddecl.py',    H,  [BOGUS, H],       ['report', H],       ['strip', H]),
    ('tools_strip.py', CC, [BOGUS, CC, 'FixtureDeadHelper'],
     ['check', CC], ['func', CC, 'FixtureDeadHelper']),
]

# The two scripts with no mode to mistype, where the bare invocation was
# itself the destructive one.  (script, fixture, correct argv, stray argvs) --
# the strays are the plausible near-misses: nothing at all, and the other
# script's calling convention.
MAIN_CASES = [
    ('autofix.py', FIXTURE_AUTOFIX, ['adjutant.cc'],
     [[], ['--apply']]),
    ('reduce_plink2.py', FIXTURE_REDUCE, ['--apply'],
     [[], ['adjutant.cc'], ['--dry-run']]),
]

# Bare invocation must print usage and exit non-zero, not raise IndexError.
NO_ARGV = ['deadcode.py', 'deadheader.py', 'deadtype.py', 'tools_strip.py',
           'closure.py', 'prune_named.py', 'extroots.py']

# deaddecl.py documents `report` as its default when no mode is given, and
# report is read-only.  That default stays; it just must not be reachable by
# *mistyping* a mode, which is what the CASES table above covers.
DEFAULTS_TO_REPORT = ['deaddecl.py']

failures = []


def check(cond, label):
    print('%s  %s' % ('PASS' if cond else 'FAIL', label))
    if not cond:
        failures.append(label)


def run(script, argv, cwd):
    return subprocess.run([sys.executable, os.path.join(TOOLS, script)] + argv,
                          cwd=cwd, capture_output=True, text=True)


def workdir(tmp, name):
    """A scratch tree holding just the fixture pair.

    deadheader/deadtype/deaddecl glob `adjutant*.cc` and `adjutant*.h` out of
    the cwd, so the fixture has to be named to match and live alone.
    """
    d = os.path.join(tmp, name)
    os.makedirs(d)
    open(os.path.join(d, 'adjutant_fixture.cc'), 'w').write(FIXTURE_CC)
    open(os.path.join(d, 'adjutant_fixture.h'), 'w').write(FIXTURE_H)
    return d


def maindir(tmp, name, body):
    """A scratch tree whose one source file is named `adjutant.cc`.

    Deliberately not workdir(): autofix.py and reduce_plink2.py write to that
    exact name, and adding it to workdir()'s tree would hand the `adjutant*.cc`
    globbers above a second file to chew on.
    """
    d = os.path.join(tmp, name)
    os.makedirs(d)
    open(os.path.join(d, 'adjutant.cc'), 'w').write(body)
    return d


def main():
    tmp = tempfile.mkdtemp(prefix='cli_guards_')
    try:
        for script, target, bogus_argv, ro_argv, rw_argv in CASES:
            # Control: the destructive mode really does rewrite this fixture.
            # Without this, "file unchanged" below could pass on an inert
            # fixture and prove nothing.
            d = workdir(tmp, script + '.control')
            before = open(os.path.join(d, target)).read()
            run(script, rw_argv, d)
            after = open(os.path.join(d, target)).read()
            check(after != before,
                  '%-15s %-6s rewrites the fixture (control)' % (script, rw_argv[0]))

            # The bug: an unrecognised mode must refuse, not fall through.
            d = workdir(tmp, script + '.bogus')
            before = open(os.path.join(d, target)).read()
            r = run(script, bogus_argv, d)
            after = open(os.path.join(d, target)).read()
            check(after == before,
                  '%-15s %-6s leaves the file untouched' % (script, BOGUS))
            check(r.returncode != 0,
                  '%-15s %-6s exits non-zero' % (script, BOGUS))
            check('usage' in (r.stdout + r.stderr).lower(),
                  '%-15s %-6s prints usage' % (script, BOGUS))

            # Guard the documented read-only mode against the fix.
            d = workdir(tmp, script + '.ro')
            before = open(os.path.join(d, target)).read()
            r = run(script, ro_argv, d)
            after = open(os.path.join(d, target)).read()
            check(after == before and r.returncode == 0,
                  '%-15s %-6s still works, writes nothing' % (script, ro_argv[0]))

        for script, fixture, ok_argv, stray_argvs in MAIN_CASES:
            # Control: invoked properly, it really does rewrite adjutant.cc.
            # Same reason as above -- otherwise "untouched" proves nothing.
            d = maindir(tmp, script + '.control', fixture)
            before = open(os.path.join(d, 'adjutant.cc')).read()
            run(script, ok_argv, d)
            after = open(os.path.join(d, 'adjutant.cc')).read()
            check(after != before, '%-16s %-13s rewrites adjutant.cc (control)'
                  % (script, ' '.join(ok_argv)))

            for n, argv in enumerate(stray_argvs):
                label = ' '.join(argv) or '(none)'
                d = maindir(tmp, '%s.stray%d' % (script, n), fixture)
                before = open(os.path.join(d, 'adjutant.cc')).read()
                r = run(script, argv, d)
                after = open(os.path.join(d, 'adjutant.cc')).read()
                out = r.stdout + r.stderr
                check(after == before, '%-16s %-13s leaves adjutant.cc untouched'
                      % (script, label))
                check(r.returncode != 0 and 'usage' in out.lower(),
                      '%-16s %-13s usage, non-zero' % (script, label))
                check('Traceback' not in out,
                      '%-16s %-13s no traceback' % (script, label))

        for script in NO_ARGV:
            d = workdir(tmp, script + '.noargv')
            r = run(script, [], d)
            out = r.stdout + r.stderr
            check('Traceback' not in out,
                  '%-15s bare invocation: no traceback' % script)
            check(r.returncode != 0 and 'usage' in out.lower(),
                  '%-15s bare invocation: usage, non-zero' % script)

        for script in DEFAULTS_TO_REPORT:
            d = workdir(tmp, script + '.default')
            before = open(os.path.join(d, H)).read()
            r = run(script, [], d)
            after = open(os.path.join(d, H)).read()
            check('Traceback' not in (r.stdout + r.stderr) and after == before
                  and r.returncode == 0,
                  '%-15s bare invocation: reports, writes nothing' % script)
    finally:
        shutil.rmtree(tmp, ignore_errors=True)

    print()
    if failures:
        print('%d FAILED' % len(failures))
        for f in failures:
            print('    %s' % f)
        return 1
    print('all guards hold')
    return 0


if __name__ == '__main__':
    sys.exit(main())
