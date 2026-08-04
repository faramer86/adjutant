"""Roots for file X = symbols in X.o that some *other* .o references.

Usage:
  extroots.py <module> <file.o> [<file.o>...]   e.g. extroots.py plink2_common *.o
"""
import subprocess, sys, re
if len(sys.argv) < 3:
    sys.exit(__doc__)
target = sys.argv[1]          # e.g. plink2_common
objs = sys.argv[2:]
def syms(o, kind):
    out = subprocess.run(['nm','-g',o], capture_output=True, text=True).stdout
    r=set()
    for l in out.split('\n'):
        p=l.split()
        if len(p)>=2 and p[-2]==kind:
            r.add(p[-1])
    return r
mine = syms(target+'.o','T')
needed=set()
for o in objs:
    if o == target+'.o': continue
    needed |= syms(o,'U')
roots = sorted(mine & needed)
# demangle C++ names to bare function names
names=set()
for r in roots:
    # Itanium mangling: _ZN6plink2<len><name>E...  -- honour the length prefix,
    # otherwise a greedy match swallows the parameter encoding too.
    m=re.match(r'_*_ZN6plink2(\d+)', r)
    if m:
        n=int(m.group(1)); start=m.end()
        names.add(r[start:start+n])
    else:
        m2=re.match(r'^_*([A-Za-z_]\w*)$', r)
        if m2: names.add(m2.group(1))
print(','.join(sorted(names)))
