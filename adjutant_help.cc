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
"  3  HARMONISE  weights + a .bim     ->  .matched, .match_report\n"
"\n"
"1 and 2 run in one command, and 2 can start from a ready-made .sscore\n"
"instead of genotypes.  3 runs on its own.\n"
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
"  harmonise   --bim --pgs-model --pgs-models-list ['keep-allele']\n"
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
"Each flag is summarised below.  For its notes and formats as well:\n"
"  adjutant --help <flag>\n"
, stdout);
  }
  // Plain --help prints each flag's summary.  --help <flag> also prints the
  // notes after it: the HelpPrint calls under `if (param_ct)` repeat the
  // flag names of the summary they follow, so they print only when asked.
  do {
    if (!param_ct) {
      fputs("\nINPUT   One of --pfile / --bfile (genotypes) or --sscore (ready-made scores).\n\n", stdout);
    }
    HelpPrint("pfile\0bfile\0", &help_ctrl, 1,
"  --pfile <prefix> ['vzs']  .pgen + .pvar + .psam.  'vzs': the .pvar is\n"
"                            Zstd-compressed.\n"
"  --bfile <prefix> ['vzs']  .bed + .bim + .fam.\n"
               );
    HelpPrint("sscore\0", &help_ctrl, 1,
"  --sscore <file>           A score matrix in place of genotypes: a .sscore\n"
"                            from plink2 or adjutant, joined to --pheno and\n"
"                            --covar on IID.  Needs --give-best-pgs or\n"
"                            --give-corr-matrix.\n"
               );
    if (param_ct) {
      HelpPrint("sscore\0", &help_ctrl, 1,
"    Column names, minus any _AVG or _SUM suffix, become the PGS IDs; a\n"
"    suffix stays where dropping it would make two names alike.\n"
                 );
    }

    if (!param_ct) {
      fputs("\nSCORING\n\n", stdout);
    }
    HelpPrint("score\0", &help_ctrl, 1,
"  --score <file> [i] [j] [k] [modifiers...]\n"
"                            Apply a linear score.  One line per (variant,\n"
"                            allele): variant ID in column i (default 1),\n"
"                            allele in j (default i+1), weight in k (default\n"
"                            j+1).\n"
"\n"
"    Modifiers:\n"
"      header-read           first line names the scores (else SCORE1, ...)\n"
"      header                skip the first line\n"
"      no-mean-imputation    drop missing genotypes instead of imputing them\n"
"      center                shift genotypes to mean 0\n"
"      variance-standardize  ...and to variance 1\n"
"      dominant / recessive  dosage>1 -> 1 / max(dosage-1, 0); not on chrX\n"
"      se                    weights are standard errors; report the score's\n"
"                            standard error\n"
"      zs                    Zstd-compress the output\n"
"      list-variants[-zs]    also write <prefix>.sscore.vars\n"
"      ignore-dup-ids        keep the first of duplicated variant IDs\n"
"      cols=[+/-]<set>       add or drop output columns; +scoresums adds the\n"
"                            unaveraged sums\n"
"\n"
               );
    HelpPrint("score-col-nums\0score\0", &help_ctrl, 1,
"  --score-col-nums <n...>   Score several weight columns in one pass, e.g.\n"
"                            3-12, each to an output column named from the\n"
"                            header.\n"
               );
    HelpPrint("score-partitions\0", &help_ctrl, 1,
"  --score-partitions <file> Also score subsets of each model's variants, in\n"
"                            the same pass.\n"
               );
    if (param_ct) {
      HelpPrint("score-partitions\0", &help_ctrl, 1,
"    <file> has a variant ID column, then one column of 0s and 1s per\n"
"    partition, named by the header:\n"
"      SNP     PART_A  PART_B\n"
"      rs1231  0       1\n"
"\n"
"    The .sscore gets a column per (model, partition), by model, total first:\n"
"      W1_AVG  W1_PART_A_AVG  W1_PART_B_AVG  W2_AVG ...\n"
"\n"
"    A partition's _AVG divides by its own allele count, so partitions of\n"
"    different sizes compare but do not add up to the total; cols=+scoresums\n"
"    gives sums.  ALLELE_CT, DENOM and NAMED_ALLELE_DOSAGE_SUM describe the\n"
"    total.\n"
                 );
    }
    HelpPrint("read-freq\0", &help_ctrl, 1,
"  --read-freq <file>        Allele frequencies for mean-imputation, instead\n"
"                            of computing them from the data.\n"
               );

    if (!param_ct) {
      fputs("\nFILTERS   Optional.\n\n", stdout);
    }
    HelpPrint("keep\0remove\0extract\0exclude\0", &help_ctrl, 1,
"  --keep <f...> / --remove <f...>\n"
"                            Keep / drop the samples listed, FID then IID on\n"
"                            each line; an IID-only list matches nothing.\n"
"  --extract <f...> / --exclude <f...>\n"
"                            Keep / drop the variants listed, one ID per line.\n"
               );
    HelpPrint("exclude-mhc-hg19\0exclude-mhc-hg38\0exclude-xmhc-hg19\0exclude-xmhc-hg38\0", &help_ctrl, 1,
"  --exclude-mhc-hg19 / --exclude-mhc-hg38\n"
"                            Drop the chr6 MHC, on the genome build named.\n"
"  --exclude-xmhc-hg19 / --exclude-xmhc-hg38\n"
"                            Drop the extended MHC (not a superset of MHC).\n"
               );
    if (param_ct) {
      HelpPrint("exclude-mhc-hg19\0exclude-mhc-hg38\0exclude-xmhc-hg19\0exclude-xmhc-hg38\0", &help_ctrl, 1,
"    Regions, bounds inclusive:\n"
"      MHC   hg19 28477797-33448354   hg38 28510120-33480577\n"
"      xMHC  hg19 25726291-33368421   hg38 25726063-33400644\n"
"\n"
"    One region flag per run, on --score or on --bim/--pgs-model, not on\n"
"    --sscore.  The build is not detected, so the wrong flag drops the wrong\n"
"    window.  xMHC starts further 5' than MHC but ends ~80kb sooner 3', so it\n"
"    can drop fewer variants.\n"
                 );
    }

    if (!param_ct) {
      fputs("\nPHENOTYPE AND COVARIATES   Required by --give-best-pgs.\n\n", stdout);
    }
    HelpPrint("pheno\0pheno-name\0pheno-col-nums\0", &help_ctrl, 1,
"  --pheno <file>            Load phenotypes, added to any in the .psam/.fam,\n"
"                            so --pheno-name is usually wanted too.  A\n"
"                            categorical phenotype is skipped with a warning.\n"
"  --pheno-name <name...>    Select phenotype columns by name...\n"
"  --pheno-col-nums <n...>   ...or by 1-based index, not both.  Each selected\n"
"                            phenotype gets its own block of output rows.\n"
               );
    HelpPrint("covar\0covar-name\0", &help_ctrl, 1,
"  --covar <file>            Load covariates.  Without them the null model is\n"
"                            intercept-only.\n"
"  --covar-name <name...>    Select covariate columns by name.\n"
               );

    if (!param_ct) {
      fputs("\nPGS EVALUATION\n\n", stdout);
    }
    HelpPrint("give-best-pgs\0", &help_ctrl, 1,
"  --give-best-pgs [no-na | cols=<names>]\n"
"                            Rank scores by what each adds to a covariates-\n"
"                            only null model.  Writes <prefix>.best_pgs, one\n"
"                            row per (phenotype, model), sorted by R2_INCR.\n"
"                            Needs a phenotype and --score or --sscore.\n"
"                            'no-na' drops the columns a run cannot fill;\n"
"                            'cols=' picks the columns and their order.\n"
               );
    if (param_ct) {
      HelpPrint("give-best-pgs\0", &help_ctrl, 1,
"    Quantitative phenotypes are fitted by OLS, case/control by logistic\n"
"    regression.  The columns (the README has the full reference):\n"
"      PGS_ID PHENO TYPE N N_CASE BETA OR SE P R2_NULL R2_FULL R2_INCR\n"
"      R2_INCR_L95 R2_INCR_U95 R2L_INCR R2L_INCR_L95 R2L_INCR_U95\n"
"      AUC_NULL AUC_FULL AUC_INCR\n"
"    R2_INCR is observed-scale for both trait types, not a pseudo-R2.  R2L_*\n"
"    need --prevalence; OR and AUC_* are case/control only.\n"
"\n"
"    no-na         drop N_CASE, OR and AUC_* when no phenotype is\n"
"                  case/control, and R2L_* then too or without --prevalence.\n"
"                  The command line decides, not the data: a column a failed\n"
"                  fit left empty still shows.\n"
"    cols=<names>  emit exactly these columns, in this order, e.g.\n"
"                  cols=PGS_ID,PHENO,R2_INCR,P.  Names as above; an unknown\n"
"                  or repeated one is an error.  Use it to parse by\n"
"                  position: no later version can move these columns.\n"
"\n"
"    Over many candidates the winner is partly the luckiest, and\n"
"    R2_INCR_L95/U95 hold for a model chosen in advance.  To quote a number\n"
"    for the winner, select on one subset and report on another.\n"
                 );
    }
    HelpPrint("give-corr-matrix\0", &help_ctrl, 1,
"  --give-corr-matrix        Pearson correlations between the scored models\n"
"                            -> <prefix>.pgs_corr.  Needs only scores.\n"
               );
    if (param_ct) {
      HelpPrint("give-corr-matrix\0", &help_ctrl, 1,
"    A square matrix, 12 significant digits, over every sample (not the\n"
"    complete cases behind any one R2).  A constant column is NA, not 0.\n"
"    The header line starts with '#', so read it as:\n"
"      R:      read.table(f, header=TRUE, comment.char=\"\", row.names=1)\n"
"      pandas: read_csv(f, sep='\\t', index_col=0)\n"
                 );
    }
    HelpPrint("pgs-pheno-pairs\0", &help_ctrl, 1,
"  --pgs-pheno-pairs <file>  Evaluate only the listed (phenotype, PGS) pairs.\n"
               );
    if (param_ct) {
      HelpPrint("pgs-pheno-pairs\0", &help_ctrl, 1,
"    Two whitespace-separated columns, '#' comments:\n"
"      QT1      PGS000123\n"
"      CC_DIAB  PGS000789\n"
"    Every phenotype still loads, so this excludes --pheno-name and\n"
"    --pheno-col-nums.  An unknown name is an error.\n"
                 );
    }
    HelpPrint("prevalence\0", &help_ctrl, 1,
"  --prevalence <val>        Population prevalence in (0,1), for liability-\n"
"                            scale R^2 (Lee 2012).  Case/control only.\n"
               );

    if (!param_ct) {
      fputs("\nHARMONISATION   Standalone.  Aligns weight files to a target .bim.\n\n", stdout);
    }
    HelpPrint("bim\0pgs-model\0pgs-models-list\0", &help_ctrl, 1,
"  --bim <file>              Target .bim; its alleles set the output coding.\n"
"  --pgs-model <file> <chr> <pos> <effect> <other> <weight> [<freq>]\n"
"  --pgs-model <file> id <id> <effect> <other> <weight> [<freq>]\n"
"                            A weight file and the 1-based column of each\n"
"                            field, matched on chr:pos or, after 'id', on the\n"
"                            .bim's variant ID.  Weights must be betas.  A\n"
"                            trailing 'keep-allele' keeps each weight on its\n"
"                            own effect allele.\n"
"  --pgs-models-list <list> <chr> <pos> <effect> <other> <weight> [<freq>]\n"
"  --pgs-models-list <list> id <id> <effect> <other> <weight> [<freq>]\n"
"                            The same for a file of weight-file paths, merged\n"
"                            into one matrix with a column per file.\n"
               );
    if (param_ct) {
      HelpPrint("bim\0pgs-model\0pgs-models-list\0", &help_ctrl, 1,
"    <freq> is only for --resolve-ambiguous.  With --pgs-models-list the\n"
"    column numbers apply to every file; each column is named from the\n"
"    file's basename, with 0 where a model lacks the variant.\n"
"\n"
"    Matching.  Once the position or ID matches, alleles are tried exact,\n"
"    reversed, flipped, then flipped+reversed, and the beta is negated when\n"
"    effect and other swap (as in bigsnpr::snp_match()).  Where the .bim has\n"
"    several records at one chr:pos, as for a split multi-allelic site, a\n"
"    variant goes to the record its alleles fit.  A source row matching a\n"
"    record an earlier row already took is DUPLICATE.\n"
"\n"
"    'keep-allele' after the column numbers writes each weight on the\n"
"    source's own effect allele at its own sign (the complement after a\n"
"    strand flip), so raw score sums equal the source model's.  By default a\n"
"    weight moves to the .bim's A1, negated, and the sums shift by 2*beta per\n"
"    reversed variant.  With --pgs-models-list, models that weight a variant\n"
"    on different alleles get one row per allele.  .match_report STATUS\n"
"    still says REVERSED or FLIPPED; OUT_A1 and OUT_BETA are what was written.\n"
"\n"
"    Weights must be betas: negating an odds ratio is meaningless.  A warning\n"
"    fires if every weight is positive with mean in [0.5, 2.0].\n"
"\n"
"    'id' survives a build mismatch but trusts rsIDs: a merged or retired\n"
"    rsID can miss, or match the wrong variant with nothing to show it.\n"
"    Check the match rate the run prints.  '.' IDs are skipped on both\n"
"    sides; where the .bim repeats an ID the first row wins and the rest go\n"
"    to <prefix>.duplicate_ids.\n"
"\n"
"    <prefix>.matched is CHR POS ID A1 A2, then a weight column per model,\n"
"    every field from the matched .bim row.  Feed it to --score with column\n"
"    3 as the ID, 4 as the allele, and --score-col-nums for the weights.\n"
"    <prefix>.match_report gives each variant's STATUS and both sides of the\n"
"    join: SRC_ID/SRC_CHR/SRC_POS as supplied, BIM_CHR/BIM_POS as matched.\n"
                 );
    }
    HelpPrint("allow-ambiguous\0resolve-ambiguous\0", &help_ctrl, 1,
"  --allow-ambiguous         Keep strand-ambiguous (A/T, C/G) variants on\n"
"                            their exact match; safe only if both files share\n"
"                            a strand.\n"
"  --resolve-ambiguous [maf] Decide them by allele frequency instead.  Not\n"
"                            with --allow-ambiguous.\n"
               );
    if (param_ct) {
      HelpPrint("allow-ambiguous\0resolve-ambiguous\0", &help_ctrl, 1,
"    With --allow-ambiguous, variants on opposite strands keep a beta of the\n"
"    wrong sign.  --resolve-ambiguous needs <freq> on --pgs-model, --read-freq\n"
"    for the target, and comparable ancestries: frequencies cannot tell a\n"
"    strand difference from a population one.  A variant stays unresolved\n"
"    when its source or target MAF reaches [maf] (default 0.4).\n"
                 );
    }

    if (!param_ct) {
      fputs("\nRUNTIME\n\n", stdout);
    }
    HelpPrint("out\0threads\0memory\0seed\0silent\0help\0version\0", &help_ctrl, 1,
"  --out <prefix>            Output prefix.  Default 'adjutant'.\n"
"  --threads <n>             Compute threads.  Default: detected cores.\n"
"  --memory <MiB>            Main workspace size.\n"
"  --seed <n>                Random seed.\n"
"  --silent                  Suppress stdout; the .log is still written.\n"
"  --help [flag...]          This text, or a flag's summary and notes.\n"
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
