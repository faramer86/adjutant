#!/usr/bin/env python3
"""Ask the linker what is dead, instead of guessing from the source.

deadcode.py reasons about one translation unit at a time.  Two kinds of dead
code slip through it:

  1. a *ring* of dead functions that call each other.  adjutant_matrix.cc's
     matrix-inversion family is thirteen functions, every one of them
     referenced -- by another member of the ring.  Nothing outside calls any of
     them.  Per-TU reachability says "0 dead"; the linker strips all thirteen.

  2. a function whose only call site is in a branch that cannot be taken.
     `if (randmem)` in adjutant.cc, where randmem is a literal 0 that is never
     assigned, keeps RandomizeBigstack looking alive to any grep.

The fix is to link twice and diff the two `-Wl,-dead_strip -Wl,-map` reports:

    -O2                dead = (never called) OR (inlined away)
    -O0 -fno-inline    dead = (never called)                       <- no inlining

    dead in BOTH   -> TIER 1, genuinely uncalled.  Safe to delete.
    dead at -O2 only -> TIER 2, either inlined (keep!) or only reachable
                        through a folded-away branch (delete).  Needs a human.

Tier 2 is a *question*, not an answer -- Count12Vec6 and RandomizeBigstack land
in the same bucket and only one of them is dead.  Read the call site.

Usage:
  linkmap.py build            # produce both maps under /tmp  (~1 min)
  linkmap.py report           # classify, using the maps already built
"""
import collections
import glob
import os
import re
import subprocess
import sys

O2_MAP = '/tmp/adjutant_o2.map'
NI_MAP = '/tmp/adjutant_noinline.map'
NI_OBJ = '/tmp/adjutant_noinline_obj'

VENDOR = ('libdeflate/lib/*.o libdeflate/lib/arm/*.o libdeflate/lib/x86/*.o '
          'zstd/lib/common/*.o zstd/lib/compress/*.o zstd/lib/decompress/*.o')
# Our own code: the application plus the bundled plink2/pgenlib units.  The
# vendored zstd/libdeflate trees are out of scope -- they are unmodified
# upstream and we do not prune them.
OURS = re.compile(r'^(adjutant[\w]*|plink2_\w+|pgenlib_\w+|SFMT)\.o$')

CXX = ('-std=c++17 -g -ffp-contract=off -DZSTD_MULTITHREAD -DIGNORE_BUNDLED_ZSTD '
       '-Ilibdeflate -I/opt/homebrew/include -I/usr/local/include '
       '-DACCELERATE_NEW_LAPACK -w')


def sh(cmd):
    return subprocess.run(cmd, shell=True, capture_output=True, text=True)


def link(objs, out, mapfile):
    r = sh('c++ %s -o %s -framework Accelerate -ldl -lpthread -lz '
           '-Wl,-dead_strip -Wl,-map,%s' % (' '.join(objs), out, mapfile))
    if not os.path.exists(mapfile):
        sys.exit('link failed:\n' + r.stderr[:2000])


def build():
    vendor = sum((glob.glob(p) for p in VENDOR.split()), [])
    if not vendor:
        sys.exit('vendored objects missing -- run `make` first')

    # -O2: whatever the normal build produced.
    normal = sorted(glob.glob('adjutant*.o') + glob.glob('include/*.o'))
    if not normal:
        sys.exit('no objects -- run `make` first')
    link(normal + vendor, '/tmp/adjutant_o2.bin', O2_MAP)
    print('  %s' % O2_MAP)

    # -O0 -fno-inline: recompile only our own units, leave the vendored trees.
    os.makedirs(NI_OBJ, exist_ok=True)
    srcs = sorted(glob.glob('adjutant*.cc') + glob.glob('include/plink2_*.cc')
                  + glob.glob('include/pgenlib_*.cc'))
    procs = []
    for f in srcs:
        o = os.path.join(NI_OBJ, f.replace('/', '_')[:-3] + '.o')
        procs.append((f, subprocess.Popen(
            'c++ -c -O0 -fno-inline -fno-inline-functions %s -o %s %s' % (CXX, o, f),
            shell=True, stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)))
    for f, p in procs:
        if p.wait():
            sys.exit('compile failed: %s\n%s' % (f, p.stderr.read().decode()[:1500]))
    sh('cc -c -O0 -fno-inline -std=gnu99 -g -DZSTD_MULTITHREAD -w -Ilibdeflate '
       '-o %s/SFMT.o include/SFMT.c' % NI_OBJ)
    link(sorted(glob.glob(NI_OBJ + '/*.o')) + vendor,
         '/tmp/adjutant_noinline.bin', NI_MAP)
    print('  %s' % NI_MAP)


def parse(path):
    files, mode, dead, live = {}, None, [], []
    for ln in open(path).read().split('\n'):
        if ln.startswith('# Object files:'):
            mode = 'obj'; continue
        if ln.startswith('# Sections:'):
            mode = None; continue
        if ln.startswith('# Symbols:'):
            mode = 'sym'; continue
        if ln.startswith('# Dead Stripped Symbols:'):
            mode = 'dead'; continue
        if not ln.strip() or ln.startswith('#'):
            continue
        if mode == 'obj':
            m = re.match(r'\[\s*(\d+)\]\s+(.*)', ln)
            if m:
                files[int(m.group(1))] = m.group(2).strip().rsplit('/', 1)[-1]
        elif mode == 'dead':
            m = re.match(r'<<dead>>\s+(0x[0-9a-fA-F]+)\s+\[\s*(\d+)\]\s+(.*)', ln)
            if m:
                dead.append((int(m.group(1), 16),
                             files.get(int(m.group(2)), '?'), m.group(3).strip()))
        elif mode == 'sym':
            m = re.match(r'0x[0-9a-fA-F]+\s+(0x[0-9a-fA-F]+)\s+\[\s*(\d+)\]\s+(.*)', ln)
            if m:
                live.append((int(m.group(1), 16),
                             files.get(int(m.group(2)), '?'), m.group(3).strip()))
    return dead, live


def junk(n):
    return ('literal' in n or n.startswith('l_') or 'GCC_except' in n
            or n.startswith('_OBJC') or 'unwind' in n or n.startswith('ltmp'))


def classify(path):
    dead, live = parse(path)
    def ours(f):
        return bool(OURS.match(f.rsplit('/', 1)[-1]))
    # A symbol that survives in ANY object is live: the entries stripped from
    # the others are duplicate weak/COMDAT copies, not dead code.
    liveset = set(n for _, f, n in live if ours(f) and not junk(n))
    deadset = {}
    for s, f, n in dead:
        if ours(f) and not junk(n) and n not in liveset:
            deadset.setdefault(n, (s, f.rsplit('/', 1)[-1]))
    return deadset


def demangle(names):
    out = subprocess.run(['c++filt'], input='\n'.join(names),
                         capture_output=True, text=True).stdout.split('\n')
    return {n: re.sub(r'\(.*', '', re.sub(r'^plink2::', '', d))
            for n, d in zip(names, out)}


def report():
    for p in (O2_MAP, NI_MAP):
        if not os.path.exists(p):
            sys.exit('%s missing -- run `linkmap.py build` first' % p)
    o2, ni = classify(O2_MAP), classify(NI_MAP)
    both = sorted(set(o2) & set(ni), key=lambda n: -o2[n][0])
    only = sorted(set(o2) - set(ni), key=lambda n: -o2[n][0])
    dem = demangle(both + only)

    for title, group, note in (
        ('TIER 1 -- never called in either build.  Safe to delete.', both, ''),
        ('TIER 2 -- dead at -O2 only: INLINED (keep) or unreachable branch '
         '(delete).', only, '  read each call site before acting')):
        print('=' * 76)
        print(title)
        if note:
            print(note)
        print('=' * 76)
        by = collections.defaultdict(list)
        for n in group:
            size, f = o2[n]
            by[f].append((size, dem.get(n, n)))
        for f in sorted(by):
            print('\n  --- %s ---' % f)
            for size, name in sorted(by[f], key=lambda x: -x[0]):
                print('    %6dB  %s' % (size, name))
        print('\n  (%d symbols)\n' % len(group))


if __name__ == '__main__':
    cmd = sys.argv[1] if len(sys.argv) > 1 else 'report'
    if cmd == 'build':
        build()
    elif cmd == 'report':
        report()
    else:
        sys.exit(__doc__)
