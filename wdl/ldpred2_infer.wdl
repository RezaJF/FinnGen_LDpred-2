## ldpred2_infer.wdl
##
## LDpred2-auto weight INFERENCE. Two paths for the LD reference:
##   (a) IN-SAMPLE (default): scatter over per-chromosome bgens to compute an
##       in-cohort LD reference (build_ld_chr), then infer per trait.
##   (b) PRECOMPUTED: pass `precomputed_ld_parts` (per-chr .rds parts, or a
##       one-file bundle as a single-element array) and skip the build.
## Then scatter over `sumstats_files` -> one harmonised weight file per trait,
## ready for ldpred2_score.wdl.
##
## WDL syntax: draft-2 (no `version` line; ${var} interpolation; optional args
## via ${"--flag " + optvar}), matching ldpred2_score.wdl.

task build_ld_chr {
    File bgen
    File bgi = bgen + ".bgi"
    File sample = bgen + ".sample"
    File variants            # variant list to include (chr,pos,a0,a1), e.g. HapMap3+
    Int n_ld = 10000
    Int window_kb = 3000
    String docker
    Int cpu = 4
    String mem = "32 GB"
    Int disk_gb = ceil(size(bgen, "G")) + 30
    Int preemptible = 1
    String chr_tag = basename(bgen, ".bgen")

    command <<<
        set -euxo pipefail
        export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
        Rscript /scripts/build_ld_reference.R \
            --bgen ${bgen} --sample ${sample} --variants ${variants} \
            --n_ld ${n_ld} --window_kb ${window_kb} \
            --out ${chr_tag}.ldpart.rds --ncores ${cpu}
    >>>
    output { File part = "${chr_tag}.ldpart.rds" }
    runtime {
        docker: "${docker}"
        cpu: cpu
        memory: mem
        disks: "local-disk " + disk_gb + " HDD"
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        preemptible: preemptible
        noAddress: true
    }
}

task infer {
    File sumstats
    Array[File] ld_parts
    String docker
    # sumstats column mapping (defaults match ldpred2_infer.R)
    String chr_col = "chr"
    String pos_col = "pos"
    String a1_col = "a1"
    String a0_col = "a0"
    String beta_col = "beta"
    String se_col = "beta_se"
    String n_col = "n_eff"
    String? freq_col
    String? info_col
    Float? n_eff_const
    Int n_chains = 30
    Int cpu = 8
    String mem = "32 GB"
    Int disk_gb = ceil(size(ld_parts, "G")) + ceil(size(sumstats, "G")) + 20
    Int preemptible = 1
    String trait_name = basename(basename(basename(sumstats, ".gz"), ".tsv"), ".txt")

    command <<<
        set -euxo pipefail
        export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
        # Output PLAIN .tsv weights (not .gz): the scoring image's data.table is
        # built without zlib, so score_chunk.R must read uncompressed weight files.
        Rscript /scripts/ldpred2_infer.R \
            --sumstats ${sumstats} \
            --ld_parts ${sep="," ld_parts} \
            --out ${trait_name}.weights.tsv \
            --trait_name ${trait_name} \
            --chr_col ${chr_col} --pos_col ${pos_col} \
            --a1_col ${a1_col} --a0_col ${a0_col} \
            --beta_col ${beta_col} --se_col ${se_col} --n_col ${n_col} \
            ${"--freq_col " + freq_col} ${"--info_col " + info_col} \
            ${"--n_eff " + n_eff_const} \
            --n_chains ${n_chains} --ncores ${cpu} \
            --report ${trait_name}.report.tsv
    >>>
    output {
        File weights = "${trait_name}.weights.tsv"
        File report  = "${trait_name}.report.tsv"
    }
    runtime {
        docker: "${docker}"
        cpu: cpu
        memory: mem
        disks: "local-disk " + disk_gb + " HDD"
        zones: "europe-west1-b europe-west1-c europe-west1-d"
        preemptible: preemptible
        noAddress: true
    }
}

workflow ldpred2_infer {
    String docker
    Array[File] sumstats_files

    # LD reference — provide EITHER precomputed parts OR the in-sample build inputs
    Array[File]? precomputed_ld_parts
    Array[String] ld_bgen_by_chr = []   # per-chromosome bgens (with .bgi + .sample) for in-sample LD
    File? variants                      # HapMap3+ variant list (chr,pos,a0,a1) for in-sample build
    Int n_ld = 10000

    # build in-sample LD only when no precomputed parts were given
    if (!defined(precomputed_ld_parts)) {
        scatter (b in ld_bgen_by_chr) {
            call build_ld_chr {
                input: docker = docker, bgen = b, variants = select_first([variants]), n_ld = n_ld
            }
        }
    }

    Array[File] ld_parts = select_first([precomputed_ld_parts, build_ld_chr.part])

    scatter (ss in sumstats_files) {
        call infer { input: docker = docker, sumstats = ss, ld_parts = ld_parts }
    }

    output {
        Array[File] weights = infer.weights
        Array[File] reports = infer.report
    }
}
