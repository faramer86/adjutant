#!/usr/bin/env python3
"""Reduce adjutant.cc to the --score-only command surface.

Compiler-guided: strip the Plink2Cmdline fields belonging to removed commands,
then repeatedly ask clang which flag-parse branches and command dispatches
still reference things that no longer exist, and delete those whole blocks.

Run from the repo root.  Idempotent-ish: always start from a pristine adjutant.cc.

Usage:
  reduce_plink2.py --apply      # rewrites ./adjutant.cc in place
"""
import json
import re
import subprocess
import sys

import tools_strip as T

CXX = ['c++', '-c', '-O2', '-std=c++17', '-g', '-DZSTD_MULTITHREAD',
       '-Ilibdeflate', '-DACCELERATE_NEW_LAPACK', '-ferror-limit=0',
       '-fsyntax-only']

# Flags that survive: they are filters/loaders the --score path genuinely needs.
KEEP_FLAGS = {
    '--maf', '--max-maf', '--mac', '--max-mac', '--hwe',
    '--mach-r2-filter', '--minimac3-r2-filter',
}

DEAD_HEADERS = [
    'plink2_adjust', 'plink2_export', 'plink2_export_legacy', 'plink2_fasta',
    'plink2_glm', 'plink2_glm_linear', 'plink2_glm_logistic',
    'plink2_glm_shared', 'plink2_import', 'plink2_import_legacy', 'plink2_ld',
    'plink2_merge', 'plink2_misc', 'plink2_perm', 'plink2_set',
]

DEAD_FIELDS = [
    'SortMode sample_sort_mode;', 'SortMode sort_vars_mode;',
    'GrmFlags grm_flags;', 'double grm_sparse_cutoff;', 'PcaFlags pca_flags;',
    'WriteCovarFlags write_covar_flags;',
    'PhenoTransformFlags pheno_transform_flags;',
    'FaFlags fa_flags;', 'RangeList vscore_col_idx_range_list;',
    'UpdateAllelesInfo update_alleles_info;', 'UpdateSexInfo update_sex_info;',
    'LdInfo ld_info;', 'SdiffInfo sdiff_info;',
    'KingFlags king_flags;', 'double king_cutoff;', 'double king_table_filter;',
    'double king_table_subset_thresh;',
    'FreqRptFlags freq_rpt_flags;', 'MissingRptFlags missing_rpt_flags;',
    'GenoCountsFlags geno_counts_flags;', 'HardyFlags hardy_flags;',
    'HetFlags het_flags;', 'SampleCountsFlags sample_counts_flags;',
    'RecoverVarIdsFlags recover_var_ids_flags;', 'VscoreFlags vscore_flags;',
    'AlleleAlphanumFlags allele_alphanum_flags;', 'RmDupMode rmdup_mode;',
    'GlmInfo glm_info;', 'AdjustInfo adjust_info;', 'FstInfo fst_info;',
    'PgenDiffInfo pgen_diff_info;', 'PermConfig perm_config;',
    'ExportfInfo exportf_info;', 'GwasSsfInfo gwas_ssf_info;',
    'ClumpInfo clump_info;', 'VcorInfo vcor_info;',
    'PhenoSvdInfo pheno_svd_info;', 'CheckSexInfo check_sex_info;',
    'MendelInfo mendel_info;', 'FlipInfo flip_info;',
]

# Cross-checks against commands that no longer exist.
OBSOLETE_GUARD = re.compile(
    r'[ \t]*if \(unlikely\('
    r'(?:pc\.flip_info\.subset_fname \|\| '
    r'\(make_plink2_flags & kfMakePlink2FillMissingWithRef\)'
    r'|pc\.freq_rpt_flags & kfAlleleFreqCol(?:Minimac3R2|MachR2))\)\) \{\n'
    r'(?:[^\n]*\n)*?[ \t]*goto main_ret_INVALID_CMDLINE_A;\n[ \t]*\}\n')

CASE_RE = re.compile(r"^\s{6}case '(.)':")
FLAG_RE = re.compile(
    r'^\s{8}(?:\} else )?if \((?:likely\()?strequal_k_unsafe\(flagname_p2, "([^"]*)"\)')
ALT_RE = re.compile(
    r'^\s{8}(?:\} else )?if \((?:likely\()?StrStartsWith\w*\(flagname_p2, "([^"]*)"')
DISP_RE = re.compile(r'^\s{2,6}if \((?:pcp->)?command_flags1 & (kfCommand1\w+)\)')


def errors(path='adjutant.cc'):
    out = subprocess.run(CXX + [path], capture_output=True, text=True).stderr
    return [l for l in out.split('\n') if ': error:' in l]


def err_lines(path='adjutant.cc'):
    seen = set()
    for l in errors(path):
        m = re.match(r'^%s:(\d+):' % re.escape(path), l)
        if m:
            seen.add(int(m.group(1)))
    return sorted(seen)


def owners(src):
    """line number -> ('flag', full-name) | ('disp', kfCommand1X) | None."""
    out = {}
    cur_case = cur = None
    for i, l in enumerate(src.split('\n'), start=1):
        m = CASE_RE.match(l)
        if m:
            cur_case = m.group(1)
            cur = None
        m = FLAG_RE.match(l) or ALT_RE.match(l)
        if m:
            cur = ('flag', '--' + (cur_case or '?') + m.group(1), l.strip())
        out[i] = cur
    return out


def headers(src):
    """full flag name -> exact 'if (...)' header substring."""
    out = {}
    cur_case = None
    for l in src.split('\n'):
        m = CASE_RE.match(l)
        if m:
            cur_case = m.group(1)
        m = FLAG_RE.match(l) or ALT_RE.match(l)
        if m:
            name = '--' + (cur_case or '?') + m.group(1)
            stripped = l.strip()
            out.setdefault(name, stripped[stripped.index('if ('):])
    return out


def main():
    # There is no path argument and no mode here -- the target is hardcoded
    # below -- so a bare `reduce_plink2.py` was itself the destructive
    # invocation.  Nothing is written unless the write was asked for.
    if sys.argv[1:] != ['--apply']:
        sys.exit(__doc__)
    src = open('adjutant.cc').read()
    base = T.balanced(src)

    # 1. includes of deleted headers
    for d in DEAD_HEADERS:
        src = src.replace('#include "%s.h"\n' % d, '')

    # 2. Plink2Cmdline fields for removed commands
    a = src.index('typedef struct Plink2CmdlineStruct {')
    b = src.index('} Plink2Cmdline;')
    body = src[a:b]
    for f in DEAD_FIELDS:
        body2 = re.sub(r'^[ \t]*' + re.escape(f) + r'[ \t]*\n', '', body, flags=re.M)
        if body2 == body:
            print('  warn: field not found: %s' % f)
        body = body2
    src = src[:a] + body + src[b:]

    # 3. obsolete cross-check guards inside flags we are keeping
    src, n_guard = OBSOLETE_GUARD.subn('', src)
    open('adjutant.cc', 'w').write(src)
    print('includes+fields+%d guards done' % n_guard)

    # 4. compiler-guided removal of dead flag branches / dispatches
    total_flags = total_disp = 0
    for round_no in range(1, 8):
        src = open('adjutant.cc').read()
        base = T.balanced(src)
        own = owners(src)
        hdrs = headers(src)
        bad_flags, unattributed = set(), []
        for ln in err_lines():
            o = own.get(ln)
            if o and o[0] == 'flag' and o[1] not in KEEP_FLAGS:
                bad_flags.add(o[1])
            elif o is None:
                unattributed.append(ln)
        if not bad_flags:
            print('round %d: no more removable flag branches '
                  '(%d unattributed errors remain)' % (round_no, len(unattributed)))
            break
        removed = 0
        for f in sorted(bad_flags):
            h = hdrs.get(f)
            if not h:
                continue
            try:
                src2, ok = T.cut_branch(src, h)
            except ValueError as e:
                print('  skip %s: %s' % (f, e))
                continue
            if ok and T.balanced(src2) == base:
                src = src2
                removed += 1
        open('adjutant.cc', 'w').write(src)
        total_flags += removed
        print('round %d: removed %d flag branches (%d errors left)'
              % (round_no, removed, len(errors())))

    print('TOTAL flag branches removed: %d' % total_flags)
    return 0


if __name__ == '__main__':
    sys.exit(main())
