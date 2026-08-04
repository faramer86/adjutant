# Protected core — what must never be pruned

Derived from the 16-flag surface.  Anything not reachable from these entry
points is a candidate for removal; anything reachable is protected.

## Call chain

    main()                                   adjutant.cc
      CmdlineParsePhase1/2                   adjutant_cmdline.cc
      Plink2Core()                           adjutant.cc
        LoadPsam()                           adjutant_psam.cc     <- --pfile/--bfile
        LoadPvar()                           adjutant_pvar.cc     <- --pfile/--bfile
        PgfiInitPhase1/2, PgrInit            include/pgenlib_read.cc
        KeepOrRemove()                       adjutant_filter.cc   <- --keep/--remove
        ExtractExcludeFlagNorange()          adjutant_filter.cc   <- --extract/--exclude
        LoadAlleleAndGenoCounts()            adjutant_data.cc     \ allele freqs for
        ComputeAlleleFreqs()                 adjutant_filter.cc   / mean-imputation
        ReadAlleleFreqs()                    adjutant_filter.cc   <- --read-freq
        ScoreReport()                        adjutant_matrix_calc.cc  <- --score
          CalcScoreThread()                  adjutant_matrix_calc.cc
          RowMajorMatrixMultiplyStridedIncr  adjutant_matrix.h    (BLAS) -- an
              inline over ColMajorMatrixMultiplyStridedAddassign, which is now
              the only definition left in adjutant_matrix.cc
          InitCstream/Cswrite                adjutant_compress_stream.cc
      DispHelp()                             adjutant_help.cc     <- --help

## Support layers (reachable from the above, not independently prunable)

    adjutant_common.cc      shared genotype/sample/variant helpers
    adjutant_cmdline.cc     flag tokenising, ID hash tables, range lists
    adjutant_decompress.cc  text-stream setup over plink2_text
    include/                pgenlib (.pgen/.bed decode), plink2_text (line reader),
                            plink2_bits/base/string/htable/memory/simd/thread,
                            plink2_bgzf + plink2_zstfile (compressed I/O)

## Invariants the test suites enforce

    score_check.sh    --score output byte-identical to stock plink2 (89 files)
    input_matrix.sh   36 input/filter/modifier combinations
    verify.sh         1000G panel, --pfile == --bfile
