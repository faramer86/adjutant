// This file is part of PLINK 2.0, copyright (C) 2005-2026 Shaun Purcell,
// Christopher Chang.
//
// This program is free software: you can redistribute it and/or modify it
// under the terms of the GNU General Public License as published by the Free
// Software Foundation, either version 3 of the License, or (at your option)
// any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT
// ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
// FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
// more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

#include "adjutant_help.h"

#include <string.h>

#include "include/plink2_bits.h"
#include "adjutant_cmdline.h"

#ifdef __cplusplus
namespace plink2 {
#endif

const char kCmdlineFormatStr[] = "\n"
"USAGE\n"
"  adjutant <input> --score <weights> [flags...]\n"
"  adjutant --sscore <scores> --give-best-pgs --pheno <file> [flags...]\n"
"  adjutant --bim <file> --pgs-model <weights> <chr> <pos> <eff> <oth> <w>\n"
"  adjutant --help [flag...]\n\n";

PglErr DispHelp(const char* const* argvk, uint32_t param_ct) {
  // yes, this is overkill.  But it should be a good template for other
  // command-line programs to use.
  uint32_t param_ctl = BitCtToWordCt(param_ct);
  PglErr reterr = kPglRetSuccess;
  const char** new_argv = nullptr;
  HelpCtrl help_ctrl;
  uint32_t arg_uidx;
  uint32_t arg_idx;
  uint32_t net_unmatched_ct;
  int32_t col_num;
  int32_t leading_dashes;
  help_ctrl.iters_left = param_ct? 2 : 0;
  help_ctrl.param_ct = param_ct;
  help_ctrl.unmatched_ct = param_ct;
  help_ctrl.param_slens = nullptr;
  help_ctrl.all_match_arr = nullptr;
  help_ctrl.argv = nullptr;
  if (param_ct) {
    if (unlikely(pgl_malloc(param_ct * sizeof(int32_t), &help_ctrl.param_slens) ||
                 pgl_malloc(param_ctl * 3 * sizeof(intptr_t), &help_ctrl.all_match_arr))) {
      goto DispHelp_ret_NOMEM;
    }
    leading_dashes = 0;
    for (arg_uidx = 0; arg_uidx != param_ct; ++arg_uidx) {
      if (argvk[arg_uidx][0] == '-') {
        leading_dashes = 1;
        break;
      }
    }
    if (leading_dashes) {
      if (unlikely(pgl_malloc(param_ct * sizeof(intptr_t), &new_argv))) {
        goto DispHelp_ret_NOMEM;
      }
      for (arg_uidx = 0; arg_uidx != param_ct; ++arg_uidx) {
        if (argvk[arg_uidx][0] == '-') {
          if (argvk[arg_uidx][1] == '-') {
            new_argv[arg_uidx] = &(argvk[arg_uidx][2]);
          } else {
            new_argv[arg_uidx] = &(argvk[arg_uidx][1]);
          }
        } else {
          new_argv[arg_uidx] = argvk[arg_uidx];
        }
      }
      help_ctrl.argv = new_argv;
    } else {
      help_ctrl.argv = argvk;
    }
    for (arg_idx = 0; arg_idx != param_ct; ++arg_idx) {
      help_ctrl.param_slens[arg_idx] = strlen(help_ctrl.argv[arg_idx]);
    }
    ZeroWArr(param_ctl * 3, help_ctrl.all_match_arr);
    help_ctrl.prefix_match_arr = &(help_ctrl.all_match_arr[param_ctl]);
    help_ctrl.perfect_match_arr = &(help_ctrl.all_match_arr[param_ctl * 2]);
    help_ctrl.preprint_newline = 1;
  } else {
    help_ctrl.argv = nullptr;
    fputs(kCmdlineFormatStr, stdout);
    fputs(
"adjutant is PLINK 2 cut down to polygenic scoring.\n"
"\n"
"  1  SCORE      genotypes + weights  ->  .sscore\n"
"  2  EVALUATE   scores + a phenotype ->  .best_pgs, .pgs_corr\n"
"  3  HARMONISE  weights + a .bim     ->  .matched\n"
"\n"
"1 and 2 run in one command.  2 also takes a ready-made score matrix via\n"
"--sscore, needing no genotypes.  3 is standalone: two text files in, two out.\n"
"\n"
"EXAMPLES\n"
"  Score one model:\n"
"    adjutant --pfile data --score weights.txt 1 2 header-read --out result\n"
"\n"
"  Score ten models and rank them against a trait:\n"
"    adjutant --pfile data --score weights.txt 1 2 header-read \\\n"
"             --score-col-nums 3-12 --pheno pheno.txt --pheno-name TRAIT \\\n"
"             --covar covar.txt --give-best-pgs --out result\n"
"\n"
"  Rank a score matrix you already have:\n"
"    adjutant --sscore scores.sscore --pheno pheno.txt --covar covar.txt \\\n"
"             --give-best-pgs --out result\n"
"\n"
"  Align a published PGS to your dataset, then score with it:\n"
"    adjutant --bim target.bim --pgs-model weights.txt 1 2 3 4 5 --out h\n"
"    adjutant --pfile target --score h.matched 3 4 6 header-read --out result\n"
"\n"
"FLAGS\n"
"  input       --pfile --bfile --sscore\n"
"  scoring     --score --score-col-nums --score-partitions --read-freq\n"
"  filters     --keep --remove --extract --exclude\n"
"              --exclude-mhc-hg19/hg38 --exclude-xmhc-hg19/hg38\n"
"  phenotype   --pheno --pheno-name --pheno-col-nums --covar --covar-name\n"
"  evaluation  --give-best-pgs --give-corr-matrix --pgs-pheno-pairs\n"
"              --prevalence\n"
"  harmonise   --bim --pgs-model --pgs-models-list\n"
"              --allow-ambiguous --resolve-ambiguous\n"
"  runtime     --out --threads --memory --seed --silent --help --version\n"
"\n"
"OUTPUT   (<prefix> from --out, default 'adjutant')\n"
"  .sscore        one row per sample, one column per scored model\n"
"  .best_pgs      one row per (phenotype, model), ranked\n"
"  .pgs_corr      correlation matrix across models\n"
"  .matched       harmonised weights, ready for --score\n"
"  .match_report  per-variant harmonisation status\n"
"  .log           full record of the run\n"
"\n"
"Detail on any flag:  adjutant --help <flag>\n"
, stdout);
  }
  do {
    if (!param_ct) {
      fputs("\nINPUT   One of --pfile / --bfile (genotypes) or --sscore (ready-made scores).\n\n", stdout);
    }
    HelpPrint("pfile\0bfile\0", &help_ctrl, 1,
"  --pfile <prefix> ['vzs']  .pgen + .pvar + .psam.  'vzs' = .pvar is Zstd-\n"
"                            compressed.\n"
"  --bfile <prefix> ['vzs']  .bed + .bim + .fam.\n"
               );
    HelpPrint("sscore\0", &help_ctrl, 1,
"  --sscore <file>           Load a score matrix instead of scoring genotypes:\n"
"                            a .sscore from plink2 or adjutant, joined to\n"
"                            --pheno/--covar on IID.  Column names become PGS\n"
"                            IDs, minus any _AVG / _SUM suffix -- kept where\n"
"                            dropping it would collide.\n"
"                            Needs --give-best-pgs or --give-corr-matrix.\n"
               );

    if (!param_ct) {
      fputs("\nSCORING\n\n", stdout);
    }
    HelpPrint("score\0", &help_ctrl, 1,
"  --score <file> [i] [j] [k] [modifiers...]\n"
"                            Apply a linear scoring system.  One line per\n"
"                            (variant, allele) pair: variant ID in column i\n"
"                            (default 1), allele in j (default i+1), weight in\n"
"                            k (default j+1).\n"
"\n"
"    Modifiers:\n"
"      header-read           first line holds score names (else SCORE1, ...)\n"
"      header                skip the first line without naming scores\n"
"      no-mean-imputation    drop missing genotypes instead of imputing them\n"
"      center                shift genotypes to mean zero\n"
"      variance-standardize  ...and to variance one\n"
"      dominant / recessive  dosage>1 -> 1 / max(dosage-1, 0).  Not on chrX\n"
"      se                    treat weights as standard errors; report the\n"
"                            standard error of the score\n"
"      zs                    Zstd-compress the output\n"
"      list-variants[-zs]    also write <prefix>.sscore.vars\n"
"      ignore-dup-ids        keep the first of duplicated variant IDs\n"
"      cols=[+/-]<set>       add or drop output columns; cols=+scoresums adds\n"
"                            the unaveraged sums\n"
"\n"
               );
    HelpPrint("score-col-nums\0score\0", &help_ctrl, 1,
"  --score-col-nums <n...>   Score several weight columns in one pass, e.g.\n"
"                            3-12.  Each becomes its own output column, named\n"
"                            from the header.\n"
               );
    HelpPrint("score-partitions\0", &help_ctrl, 1,
"  --score-partitions <file> Split each model into variant subsets and score\n"
"                            every subset separately, in one pass.\n"
"\n"
"    <file> is a variant ID column then one 0/1 column per partition, named by\n"
"    the header.  Values must be exactly 0 or 1.\n"
"      SNP     PART_A  PART_B\n"
"      rs1231  0         1\n"
"\n"
"    Output is the ordinary .sscore, one column per (model, partition),\n"
"    grouped by model with the unpartitioned total first:\n"
"      W1_AVG  W1_PART_A_AVG  W1_PART_B_AVG  W2_AVG ...\n"
"\n"
"    Each partition's _AVG divides by ITS OWN variants' allele count, so\n"
"    partitions of different sizes stay comparable.  They therefore do NOT sum\n"
"    to the total; use cols=+scoresums for that.  ALLELE_CT, DENOM and\n"
"    NAMED_ALLELE_DOSAGE_SUM describe the total.\n"
"\n"
               );
    HelpPrint("read-freq\0", &help_ctrl, 1,
"  --read-freq <file>        Allele frequencies for mean-imputation, instead\n"
"                            of computing them from the data.\n"
               );

    if (!param_ct) {
      fputs("\nFILTERS   Optional.  Each takes file(s) of IDs, one per line.\n\n", stdout);
    }
    HelpPrint("keep\0remove\0extract\0exclude\0", &help_ctrl, 1,
"  --keep <f...> / --remove <f...>\n"
"                            Keep / drop samples.  Two columns, FID then IID;\n"
"                            an IID-only list matches nothing.\n"
"  --extract <f...> / --exclude <f...>\n"
"                            Keep / drop variants, by ID.\n"
               );
    HelpPrint("exclude-mhc-hg19\0exclude-mhc-hg38\0exclude-xmhc-hg19\0exclude-xmhc-hg38\0", &help_ctrl, 1,
"  --exclude-mhc-hg19 / --exclude-mhc-hg38\n"
"                            Drop chr6 MHC.  hg19 28477797-33448354,\n"
"                            hg38 28510120-33480577.\n"
"  --exclude-xmhc-hg19 / --exclude-xmhc-hg38\n"
"                            Drop the extended MHC.  hg19 25726291-33368421,\n"
"                            hg38 25726063-33400644.\n"
"\n"
"  Bounds are inclusive.  One region flag per run; not valid with --sscore.\n"
"  Applies to --score and to --bim/--pgs-model alike.  Nothing auto-detects\n"
"  the build, and the wrong flag silently drops the wrong window.\n"
"\n"
"  xMHC is NOT a superset of MHC: it reaches further 5' but stops ~80kb short\n"
"  3', so it can drop fewer variants.  The two regions come from different\n"
"  published definitions.\n"
               );

    if (!param_ct) {
      fputs("\nPHENOTYPE AND COVARIATES   Required by --give-best-pgs.\n\n", stdout);
    }
    HelpPrint("pheno\0pheno-name\0pheno-col-nums\0", &help_ctrl, 1,
"  --pheno <file>            Load phenotypes.  ADDS to any in the .psam/.fam\n"
"                            rather than replacing them, so --pheno-name is\n"
"                            usually wanted too.  A categorical phenotype is\n"
"                            skipped with a warning, not an error.\n"
"  --pheno-name <name...>    Select phenotype column(s) by name.\n"
"  --pheno-col-nums <n...>   ...or by 1-based index.  Not with --pheno-name.\n"
"\n"
"    Several selections are allowed; each gets its own block of output rows.\n"
"\n"
               );
    HelpPrint("covar\0covar-name\0", &help_ctrl, 1,
"  --covar <file>            Load covariates.  Without them the null model is\n"
"                            intercept-only.\n"
"  --covar-name <name...>    Select covariate column(s) by name.\n"
               );

    if (!param_ct) {
      fputs("\nPGS EVALUATION\n\n", stdout);
    }
    HelpPrint("give-best-pgs\0", &help_ctrl, 1,
"  --give-best-pgs [no-na | cols=<names>]\n"
"                            Rank scores by how much each adds to a\n"
"                            covariates-only null model.  Writes\n"
"                            <prefix>.best_pgs, one row per (phenotype,\n"
"                            model), sorted by R2_INCR within each phenotype.\n"
"\n"
"    Needs a phenotype and either --score or --sscore.  Quantitative\n"
"    phenotypes are fitted by OLS, case/control by logistic regression.\n"
"\n"
"    20 columns; the README has the full reference.  The sort key R2_INCR is\n"
"    observed-scale for both trait types, not a pseudo-R2 of the logistic fit.\n"
"    R2L_* need --prevalence; OR and AUC_* are case/control only.\n"
"\n"
"    'no-na' omits the columns this run cannot fill, which is 8 of the 20 when\n"
"    every phenotype is quantitative.  N_CASE, OR and AUC_* go when no\n"
"    phenotype is case/control; R2L_* go when --prevalence is absent as well.\n"
"    The rule is the configuration, not the data, so the columns follow from\n"
"    the command line and a column left empty by a failed fit still appears.\n"
"\n"
"    'cols=' names the columns to emit, comma-separated, and emits them in the\n"
"    order named -- e.g. cols=PGS_ID,PHENO,R2_INCR,P.  That is the way to make\n"
"    a positional parser safe: you fix the positions, so nothing a later\n"
"    version adds or removes can move them.  An unknown or repeated name is a\n"
"    parse error.  The names are the header ones without the leading '#':\n"
"      PGS_ID PHENO TYPE N N_CASE BETA OR SE P R2_NULL R2_FULL R2_INCR\n"
"      R2_INCR_L95 R2_INCR_U95 R2L_INCR R2L_INCR_L95 R2L_INCR_U95\n"
"      AUC_NULL AUC_FULL AUC_INCR\n"
"\n"
"    Read the interval before the ranking.  Over many candidates the winner is\n"
"    partly the luckiest, and R2_INCR_L95/U95 describe a model chosen in\n"
"    advance, not one picked as the maximum of many.  To quote a number for the\n"
"    winner, select on one subset and report on another.\n"
"\n"
               );
    HelpPrint("give-corr-matrix\0", &help_ctrl, 1,
"  --give-corr-matrix        Pearson correlation matrix across the scored\n"
"                            models -> <prefix>.pgs_corr.\n"
"\n"
"    Needs only scores: no phenotype, no --give-best-pgs.  A full square\n"
"    matrix, 12 significant digits, computed over EVERY sample -- not the\n"
"    complete-case set behind any one R2.  A constant column is NA, not 0.\n"
"    The header line starts with '#', so read it as:\n"
"      R:      read.table(f, header=TRUE, comment.char=\"\", row.names=1)\n"
"      pandas: read_csv(f, sep='\\t', index_col=0)\n"
"\n"
               );
    HelpPrint("pgs-pheno-pairs\0", &help_ctrl, 1,
"  --pgs-pheno-pairs <file>  Evaluate only the listed (phenotype, PGS) pairs\n"
"                            instead of every combination.\n"
"\n"
"    Two whitespace-separated columns, '#' comments:\n"
"      QT1      PGS000123\n"
"      CC_DIAB  PGS000789\n"
"\n"
"    A filter, not a selector: all phenotypes still load, so this is mutually\n"
"    exclusive with --pheno-name / --pheno-col-nums.  An unknown name is an\n"
"    error, not a skipped row.\n"
"\n"
               );
    HelpPrint("prevalence\0", &help_ctrl, 1,
"  --prevalence <val>        Population prevalence in (0,1), for liability-\n"
"                            scale R^2 (Lee 2012).  Case/control only.\n"
               );
    if (!param_ct) {
      fputs("\nHARMONISATION   Standalone.  Aligns a weight file to a target .bim.\n\n", stdout);
    }
    HelpPrint("bim\0pgs-model\0pgs-models-list\0", &help_ctrl, 1,
"  --bim <file>              Target .bim, defining the wanted allele coding.\n"
"  --pgs-model <file> <chr> <pos> <effect> <other> <weight> [<freq>]\n"
"  --pgs-model <file> id <id> <effect> <other> <weight> [<freq>]\n"
"                            Weight file plus the 1-based column numbers of\n"
"                            each field.  <freq> is only for\n"
"                            --resolve-ambiguous.\n"
"  --pgs-models-list <list> <chr> <pos> <effect> <other> <weight> [<freq>]\n"
"  --pgs-models-list <list> id <id> <effect> <other> <weight> [<freq>]\n"
"                            The same for a file of weight-file paths, merged\n"
"                            into one matrix; the column numbers apply to all\n"
"                            of them.  Non-overlapping variants get weight 0;\n"
"                            each column is named from the file's basename.\n"
"\n"
"  --pgs-model and --pgs-models-list match variants on chr:pos, or on the\n"
"  .bim's variant ID when the column list starts with 'id'.  Alleles are then\n"
"  tried in four representations -- exact, reversed, flipped, flipped+reversed\n"
"  -- and the beta is negated when the effect and other alleles swap.\n"
"  Reimplements bigsnpr::snp_match().\n"
"\n"
"  Weights MUST be betas.  Negating an odds ratio is meaningless; a warning\n"
"  fires if every weight is positive with mean in [0.5, 2.0].\n"
"\n"
"  'id' survives a build mismatch between the weight file and the target, but\n"
"  trades it for rsID versioning: merged or retired rsIDs can fail to match,\n"
"  or match the wrong variant, with nothing in the coordinates to show it.\n"
"  Read the match rate the run prints; a low one means the two files disagree\n"
"  about what an rsID means.  Rows with a '.' ID are skipped on both sides;\n"
"  where a .bim repeats an ID the first row wins, the rest going to\n"
"  <prefix>.duplicate_ids.\n"
"\n"
"  <prefix>.matched is CHR POS ID A1 A2 then one weight column per model,\n"
"  every identifying field taken from the matched .bim row -- so it carries\n"
"  the target's coordinates, not the source's.  Feed it to --score with column\n"
"  3 as the ID, 4 as the allele, and --score-col-nums for the weights.\n"
"  <prefix>.match_report gives a per-variant STATUS with both sides of the\n"
"  join: SRC_ID/SRC_CHR/SRC_POS as supplied, BIM_CHR/BIM_POS as matched.\n"
"\n"
               );
    HelpPrint("allow-ambiguous\0resolve-ambiguous\0", &help_ctrl, 1,
"  --allow-ambiguous         Keep strand-ambiguous (A/T, C/G) variants on their\n"
"                            exact match instead of dropping them.  Assumes\n"
"                            both files are on the same strand -- if they are\n"
"                            not, those betas get the wrong sign.\n"
"  --resolve-ambiguous [maf] Decide them by allele frequency instead: needs\n"
"                            <freq> on --pgs-model and --read-freq for the\n"
"                            target, and the two must come from COMPARABLE\n"
"                            ANCESTRIES -- the comparison cannot tell a strand\n"
"                            difference from a population one.  Refuses when\n"
"                            EITHER the source or the target MAF reaches [maf]\n"
"                            (default 0.4), where the two strand hypotheses\n"
"                            are indistinguishable.\n"
"\n"
"  The two flags cannot be combined: mixing resolved and merely-kept\n"
"  palindromes in one output would hide the difference in evidence.\n"
               );

    if (!param_ct) {
      fputs("\nRUNTIME\n\n", stdout);
    }
    HelpPrint("out\0threads\0memory\0seed\0silent\0help\0version\0", &help_ctrl, 1,
"  --out <prefix>            Output prefix.  Default 'adjutant'.\n"
"  --threads <n>             Compute threads.  Default: detected cores.\n"
"  --memory <MiB>            Main workspace size.\n"
"  --seed <n>                Random seed.\n"
"  --silent                  Suppress stdout; the .log is still written.\n"
"  --help [flag...]          This text, or detail on the named flags.\n"
"  --version                 Print version and exit.\n"
               );
  } while (help_ctrl.iters_left--);
  if (help_ctrl.unmatched_ct) {
    net_unmatched_ct = help_ctrl.unmatched_ct;
    printf("\nNo help entr%s for", (help_ctrl.unmatched_ct == 1)? "y" : "ies");
    col_num = (help_ctrl.unmatched_ct == 1)? 17 : 19;
    arg_uidx = 0;
    // er, should replace the \n logic with a WordWrap() call
    while (help_ctrl.unmatched_ct) {
      arg_uidx = AdvTo0Bit(help_ctrl.all_match_arr, arg_uidx);
      help_ctrl.unmatched_ct--;
      if (help_ctrl.unmatched_ct) {
        if (net_unmatched_ct == 2) {
          if (help_ctrl.param_slens[arg_uidx] + col_num > 76) {
            putc_unlocked('\n', stdout);
            col_num = 2 + help_ctrl.param_slens[arg_uidx];
          } else {
            putc_unlocked(' ', stdout);
            col_num += 3 + help_ctrl.param_slens[arg_uidx];
          }
          putc_unlocked('\'', stdout);
          fputs(argvk[arg_uidx], stdout);
          putc_unlocked('\'', stdout);
        } else {
          if (help_ctrl.param_slens[arg_uidx] + col_num > 75) {
            putc_unlocked('\n', stdout);
            col_num = 3 + help_ctrl.param_slens[arg_uidx];
          } else {
            putc_unlocked(' ', stdout);
            col_num += 4 + help_ctrl.param_slens[arg_uidx];
          }
          putc_unlocked('\'', stdout);
          fputs(argvk[arg_uidx], stdout);
          fputs("',", stdout);
        }
        if (help_ctrl.unmatched_ct == 1) {
          if (col_num > 76) {
            fputs("\nor", stdout);
            col_num = 2;
          } else {
            fputs(" or", stdout);
            col_num += 3;
          }
        }
      } else {
        putc_unlocked((help_ctrl.param_slens[arg_uidx] + col_num > 75)? '\n' : ' ', stdout);
        putc_unlocked('\'', stdout);
        fputs(argvk[arg_uidx], stdout);
        fputs("\'.\n", stdout);
      }
      arg_uidx++;
    }
  }
  if (param_ct) {
    while (0) {
    DispHelp_ret_NOMEM:
      reterr = kPglRetNomem;
    }
    free_cond(help_ctrl.param_slens);
    free_cond(help_ctrl.all_match_arr);
    free_cond(new_argv);
  }
  return reterr;
}

#ifdef __cplusplus
}
#endif
