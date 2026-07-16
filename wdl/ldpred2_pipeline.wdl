## ldpred2_pipeline.wdl
##
## End-to-end LDpred2 pipeline in one submission:
##   summary statistics ─▶ (LD reference) ─▶ ldpred2_infer (weights)
##                                            └─▶ ldpred2_score (per-individual PGS)
##
## Chains the two modular sub-workflows. Each is also runnable standalone:
##   - wdl/ldpred2_infer.wdl   (weights from sumstats + LD reference)
##   - wdl/ldpred2_score.wdl   (per-individual PGS from weights + bgen)
##
## Cromwell needs the imported WDLs zipped as dependencies at submission, e.g.
##   ( cd wdl && zip ldpred2_deps.zip ldpred2_infer.wdl ldpred2_score.wdl )
## then submit ldpred2_pipeline.wdl with --deps ldpred2_deps.zip.
##
## WDL syntax: draft-2, matching the sub-workflows.

import "ldpred2_infer.wdl" as infer
import "ldpred2_score.wdl" as score

workflow ldpred2_pipeline {
    String docker

    # --- inference inputs ---
    Array[File] sumstats_files
    Array[File]? precomputed_ld_parts
    Array[String] ld_bgen_by_chr = []
    File? variants
    Int n_ld = 10000

    # --- scoring inputs ---
    File bgen_chunk_list

    call infer.ldpred2_infer as infer_step {
        input:
            docker = docker,
            sumstats_files = sumstats_files,
            precomputed_ld_parts = precomputed_ld_parts,
            ld_bgen_by_chr = ld_bgen_by_chr,
            variants = variants,
            n_ld = n_ld
    }

    # Score all inferred weight files (one .sscore per trait x weight file).
    call score.ldpred2_score as score_step {
        input:
            docker = docker,
            bgen_chunk_list = bgen_chunk_list,
            score_files = infer_step.weights
    }

    output {
        Array[File] weights = infer_step.weights   # LDpred2-auto weights per trait
        Array[File] reports = infer_step.reports    # per-trait h2/p QC
        Array[File] sscores = score_step.sscores    # per-individual PGS
        File pgs_qc         = score_step.qc
    }
}
